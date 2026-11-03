// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef PAPER2_TRANSPORT_DIAGNOSTICS_HH
#define PAPER2_TRANSPORT_DIAGNOSTICS_HH

#include <array>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

// Serial, rectangular YaspGrid/BOX diagnostic. A region is a union of complete
// vertex-centred control volumes; the cut is a midpoint between vertex layers.
// Fluxes use the residual's kernel, never the cell-centred VTK velocity.
template<class TypeTag>
class TransportDiagnostics
{
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using GV = GetPropType<TypeTag, Properties::GridVariables>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Solution = GetPropType<TypeTag, Properties::SolutionVector>;
    using Residual = GetPropType<TypeTag, Properties::LocalResidual>;
    struct Region {
        std::string name;
        double requested, cut;
        double initialSalt = 0, initialWater = 0, previousSalt = 0, previousWater = 0;
        double cumulativeSalt = 0, cumulativeWater = 0, cumulativeAdv = 0, cumulativeDiff = 0;
    };
    struct Budget {
        double salt = 0, water = 0;
        double total = 0, adv = 0, diff = 0, waterOut = 0;
        double otherSaltOut = 0, otherWaterOut = 0;
    };
    const Problem& problem_;
    const GG& gg_;
    GV& gv_;
    std::vector<double> dofZ_, dofX_;
    std::vector<std::size_t> dofLevel_;
    std::size_t levelCount_ = 0;
    std::vector<Region> regions_;
    std::ofstream budget_, surface_, weighted_, localFlux_;
    double bottom_, top_, width_, height_, tolerance_, threshold_;
    double lastTime_ = 0, initialAllSalt_ = 0, previousDirSalt_ = 0;
    double initialAllWater_ = 0, previousDirWater_ = 0;
    double cumulativeBottomSaltIn_ = 0, cumulativeBottomWaterIn_ = 0;
    double cumulativeExternalSaltOut_ = 0, cumulativeExternalWaterOut_ = 0;
    double surfaceLow_ = -1, surfaceHigh_ = -1, globalLow_ = -1, globalHigh_ = -1;
    std::size_t step_ = 0;
    bool initialized_ = false, finalized_ = false, stopAtThreshold_, writeLocalFluxes_;

    static void checkStream(const std::ofstream& stream) {
        if (!stream) throw std::runtime_error("Paper 2 diagnostic output write failed");
    }
    static void numberOrNull(std::ostream& out, double value) {
        if (value < 0) out << "null"; else out << value;
    }
    void writeStop(const std::string& reason) const {
        // Fixed filename in a new, dedicated run directory. No user text in JSON.
        std::ofstream out("stop_record.json.tmp");
        out << std::setprecision(17)
            << "{\n  \"schema_version\": 1,\n  \"reason\": \"" << reason
            << "\",\n  \"last_accepted_time_s\": " << lastTime_
            << ",\n  \"accepted_steps\": " << step_
            << ",\n  \"initial_diagnostics_written\": " << (initialized_ ? "true" : "false")
            << ",\n  \"requested_end_s\": " << getParam<double>("TimeLoop.TEnd")
            << ",\n  \"surface_threshold_low_s\": "; numberOrNull(out, surfaceLow_);
        out << ",\n  \"surface_threshold_high_s\": "; numberOrNull(out, surfaceHigh_);
        out << ",\n  \"global_threshold_low_s\": "; numberOrNull(out, globalLow_);
        out << ",\n  \"global_threshold_high_s\": "; numberOrNull(out, globalHigh_);
        out << "\n}\n";
        out.flush(); checkStream(out); out.close();
        std::filesystem::rename("stop_record.json.tmp", "stop_record.json");
    }

public:
    TransportDiagnostics(const Problem& problem, const GG& gg, GV& gv)
        : problem_(problem), gg_(gg), gv_(gv)
        , bottom_(gg.bBoxMin()[1]), top_(gg.bBoxMax()[1])
        , width_(gg.bBoxMax()[0]-gg.bBoxMin()[0]), height_(top_-bottom_)
        , tolerance_(1e-10*std::max(1.0, height_))
        , threshold_(getParam<double>("Transport.SolubilityMoleFraction", 0.1))
        , stopAtThreshold_(getParam<bool>("Transport.StopAtSolubility", true))
        , writeLocalFluxes_(getParam<bool>("Transport.WriteLocalFluxes", false))
    {
        static_assert(GG::GridView::dimension == 2, "Paper 2 diagnostics currently support 2D BOX only");
        if (gg.gridView().comm().size() != 1)
            throw std::runtime_error("Transport diagnostics require one MPI rank (no ghost reduction implemented)");
        if (!(width_ > 0 && height_ > 0 && threshold_ > 0 && threshold_ < 1))
            throw std::runtime_error("Invalid geometry or solubility threshold");
        if (!getParam<bool>("Problem.KrwMualemExponent", false))
            throw std::runtime_error("Paper 2 requires exponent-form krw");
        for (const auto* name : {"transport_budget.csv", "surface_events.csv", "stop_record.json", "transport_metadata.json", "weighted_online.csv", "local_vertical_flux.csv"})
            if (std::filesystem::exists(name))
                throw std::runtime_error(std::string("Refusing to overwrite diagnostic: ")+name);

        std::vector<double> levels;
        dofZ_.resize(gg.numDofs()); dofX_.resize(gg.numDofs()); dofLevel_.resize(gg.numDofs());
        for (const auto& v : vertices(gg.gridView())) {
            const auto z = v.geometry().center()[1];
            levels.push_back(z); dofZ_[gg.vertexMapper().index(v)] = z;
            dofX_[gg.vertexMapper().index(v)] = v.geometry().center()[0];
        }
        std::sort(levels.begin(), levels.end());
        levels.erase(std::unique(levels.begin(), levels.end(), [this](double a, double b) {
            return std::abs(a-b) < tolerance_;
        }), levels.end());
        if (levels.size() < 3) throw std::runtime_error("Need at least three vertical vertex layers");
        levelCount_ = levels.size();
        for (std::size_t d=0; d<dofZ_.size(); ++d)
            dofLevel_[d] = std::lower_bound(levels.begin(), levels.end(), dofZ_[d])-levels.begin();
        const auto lockedCuts = getParam<std::vector<double>>("Transport.CutCoordinates");
        if (lockedCuts.size()!=7) throw std::runtime_error("WP6 requires seven locked cut coordinates");
        std::size_t cutIndex=0;
        for (double fraction : {0.9, 0.75, 0.5, 0.995, 0.99, 0.975, 0.95}) {
            const double target = lockedCuts[cutIndex++];
            double best = 0.5*(levels[0]+levels[1]);
            for (std::size_t i = 1; i+1 < levels.size(); ++i) {
                const double candidate = 0.5*(levels[i]+levels[i+1]);
                if (std::abs(candidate-target) < std::abs(best-target)) best = candidate;
            }
            if (std::abs(best-target)>1e-12) throw std::runtime_error("Locked WP6 cut is not a conservative BOX face");
            regions_.push_back({"upper_"+std::to_string(fraction), fraction, best});
        }
        regions_.push_back({"free_domain", 0.0, 0.5*(levels[0]+levels[1])});
        budget_.open("transport_budget.csv"); surface_.open("surface_events.csv");
        weighted_.open("weighted_online.csv"); weighted_ << std::setprecision(17);
        weighted_ << "step,time_s,weighted_top_mean,weighted_top_sigma,weighted_fourier,weighted_projection,volume_perturbation_rms,harmonic2_top,harmonic3_top,projection_defined\n";
        budget_ << std::setprecision(17);
        surface_ << std::setprecision(17);
        budget_ << "step,time_s,dt_s,region,requested_z_over_H,actual_z_over_H,inventory_salt_m,inventory_initial_salt_m,inventory_water_m,"
                   "cut_salt_adv_down_m_s,cut_salt_diff_down_m_s,cut_salt_total_down_m_s,cut_water_down_m_s,"
                   "other_salt_out_m_s,other_water_out_m_s,cumulative_cut_salt_m,cumulative_cut_adv_m,cumulative_cut_diff_m,"
                   "cumulative_net_water_out_m,salt_residual_m,water_residual_m,step_salt_residual_m,step_water_residual_m,"
                   "whole_salt_inventory_m,whole_water_inventory_m,bottom_salt_reaction_in_m_s,bottom_water_reaction_in_m_s,"
                   "whole_salt_residual_m,whole_water_residual_m,bottom_reaction_defined,threshold_not_yet_observed\n";
        surface_ << "step,time_s,dt_s,x_top_area_mean,x_top_max,x_domain_max,x_domain_min,x_top_node_mean,sigma_top_node,"
                    "top_area_per_unit_thickness_m,surface_threshold_reached,global_threshold_reached\n";
        if (writeLocalFluxes_) {
            localFlux_.open("local_vertical_flux.csv"); localFlux_ << std::setprecision(17);
            localFlux_ << "step,time_s,dt_s,region,actual_z_over_H,face_x_m,face_area_m,"
                          "water_down_integrated_m2_s,water_down_density_m_s,"
                          "salt_adv_down_integrated_m2_s,salt_adv_down_density_m_s,"
                          "salt_diff_down_integrated_m2_s,salt_diff_down_density_m_s,"
                          "salt_total_down_integrated_m2_s,salt_total_down_density_m_s,upwind_mole_fraction\n";
            checkStream(localFlux_);
        }
        std::ofstream meta("transport_metadata.json");
        meta << std::setprecision(17)
             << "{\"schema_version\":2,\"inventory\":\"phi*Sw*x, per horizontal area\","
                "\"salt_inventory_unit\":\"m\",\"salt_flux_unit\":\"m/s\","
                "\"geometry\":\"2D rectangular BOX, unit out-of-plane thickness\","
                "\"time_integration\":\"accepted-step backward Euler endpoint flux\","
                "\"bottom_flux\":\"reaction inferred from prescribed bottom control volumes; not an independent boundary-flux check\","
                "\"whole_budget\":\"equivalent to free-domain balance after including prescribed-node reaction\","
                "\"vtk_velocity\":\"cell-centre Darcy reconstruction including krw; not conservative face flux\","
                "\"local_flux\":\"optional conservative SCVF contributions; sums of integrated columns divided by width reproduce seven-plane budget fluxes\","
                "\"local_flux_enabled\":" << (writeLocalFluxes_ ? "true" : "false") << ','
             << "\"solubility_threshold\":" << threshold_ << ",\"stop_at_solubility\":" << (stopAtThreshold_ ? "true" : "false")
             << ",\"cuts\":[";
        for (std::size_t i=0;i<regions_.size();++i) {
            if (i) meta << ',';
            meta << "{\"name\":\"" << regions_[i].name << "\",\"requested_z_over_H\":" << regions_[i].requested
                 << ",\"actual_z_over_H\":" << (regions_[i].cut-bottom_)/height_ << '}';
        }
        meta << "]}\n"; meta.flush(); checkStream(meta);
        checkStream(budget_); checkStream(surface_);
        writeStop("initialized_not_completed");
    }

    // If stack unwinding reaches this object, retain the last accepted step.
    // SIGKILL/power loss cannot run a destructor; the per-step checkpoint remains.
    ~TransportDiagnostics() {
        if (!finalized_) { try { writeStop("exception_or_incomplete"); } catch (...) {} }
    }

    bool thresholdReached() const { return globalHigh_ >= 0 || surfaceHigh_ >= 0; }
    bool shouldStop() const { return stopAtThreshold_ && thresholdReached(); }
    std::string thresholdReason() const {
        return surfaceHigh_ >= 0 ? "surface_solubility_threshold" : "interior_solubility_threshold";
    }
    void finish(const std::string& reason) { writeStop(reason); finalized_ = true; }

    void write(double time, double dt, const Solution& x)
    {
        if (!std::isfinite(time) || !std::isfinite(dt) ||
            (initialized_ && (!(dt > 0) || std::abs(time-lastTime_-dt) > 1e-9*std::max(1.0,time))) ||
            (!initialized_ && (time != 0 || dt != 0)))
            throw std::runtime_error("Invalid accepted-step time for transport integration");
        const std::size_t outputStep = initialized_ ? step_+1 : 0;
        std::vector<Budget> b(regions_.size());
        double allSalt=0, allWater=0, dirSalt=0, dirWater=0;
        double topWeighted=0, topArea=0;
        std::vector<double> topWeight(x.size(),0.0), volumeWeight(x.size(),0.0);
        double externalSalt=0, externalWater=0;
        double topMax=-std::numeric_limits<double>::infinity(), domainMax=topMax;
        double domainMin=std::numeric_limits<double>::infinity(), nodeSum=0;
        std::size_t nodeCount=0;
        for (std::size_t d=0;d<x.size();++d) {
            const double value=x[d][1];
            if (!std::isfinite(value) || !std::isfinite(x[d][0]))
                throw std::runtime_error("Non-finite accepted state");
            domainMin=std::min(domainMin,value); domainMax=std::max(domainMax,value);
            if (std::abs(dofZ_[d]-top_) < tolerance_) {
                topMax=std::max(topMax,value); nodeSum+=value; ++nodeCount;
            }
        }
        for (const auto& element : elements(gg_.gridView())) {
            auto geometry=localView(gg_); geometry.bind(element);
            auto vv=localView(gv_.curGridVolVars()); vv.bind(element,geometry,x);
            auto cache=localView(gv_.gridFluxVarsCache()); cache.bind(element,geometry,vv);
            for (const auto& scv : scvs(geometry)) {
                const auto storage=Residual::normalizedStorage(problem_,vv[scv]);
                const double water=storage[0]*scv.volume()/width_;
                const double salt=storage[1]*scv.volume()/width_;
                if (!std::isfinite(salt) || !std::isfinite(water))
                    throw std::runtime_error("Non-finite control-volume inventory");
                volumeWeight[scv.dofIndex()]+=scv.volume()/width_;
                const double z=dofZ_[scv.dofIndex()];
                allSalt+=salt; allWater+=water;
                if (std::abs(z-bottom_) < tolerance_) { dirSalt+=salt; dirWater+=water; }
                for (std::size_t i=0;i<b.size();++i) if (z>regions_[i].cut) {
                    b[i].salt+=salt; b[i].water+=water;
                }
            }
            for (const auto& face : scvfs(geometry)) {
                const auto& inside=geometry.scv(face.insideScvIdx());
                const double zi=dofZ_[inside.dofIndex()];
                if (face.boundary()) {
                    if (std::abs(face.ipGlobal()[1]-bottom_) < tolerance_) continue; // Dirichlet, handled by reaction
                    const auto flux=problem_.neumann(element,geometry,vv,cache,face);
                    const double salt=flux[1]*face.area()/width_;
                    const double water=flux[0]*face.area()/width_;
                    externalSalt+=salt; externalWater+=water;
                    for (std::size_t i=0;i<b.size();++i) if (zi>regions_[i].cut) {
                        b[i].otherSaltOut+=salt; b[i].otherWaterOut+=water;
                    }
                    if (std::abs(face.ipGlobal()[1]-top_) < tolerance_) {
                        topArea+=face.area(); topWeighted+=face.area()*x[inside.dofIndex()][1];
                        topWeight[inside.dofIndex()]+=face.area();
                    }
                    continue;
                }
                const double zo=dofZ_[geometry.scv(face.outsideScvIdx()).dofIndex()];
                const auto flux=Residual::transportFlux(problem_,element,geometry,vv,face,cache);
                for (std::size_t i=0;i<b.size();++i) {
                    const int sign=int(zi>regions_[i].cut)-int(zo>regions_[i].cut);
                    b[i].adv+=sign*flux.advective/width_; b[i].diff+=sign*flux.diffusive/width_;
                    b[i].total+=sign*flux.total/width_; b[i].waterOut+=sign*flux.water/width_;
                    if (writeLocalFluxes_ && sign) {
                        const double area=face.area();
                        localFlux_ << outputStep << ',' << time << ',' << dt << ',' << regions_[i].name << ','
                                   << (regions_[i].cut-bottom_)/height_ << ',' << face.ipGlobal()[0] << ',' << area
                                   << ',' << sign*flux.water << ',' << sign*flux.water/area
                                   << ',' << sign*flux.advective << ',' << sign*flux.advective/area
                                   << ',' << sign*flux.diffusive << ',' << sign*flux.diffusive/area
                                   << ',' << sign*flux.total << ',' << sign*flux.total/area
                                   << ',' << flux.upwindMoleFraction << '\n';
                    }
                }
            }
        }
        if (!nodeCount || std::abs(topArea-width_) > 1e-9*width_)
            throw std::runtime_error("Top boundary quadrature does not cover the horizontal width");
        double variance=0; const double nodeMean=nodeSum/nodeCount;
        for (std::size_t d=0;d<x.size();++d) if (std::abs(dofZ_[d]-top_) < tolerance_)
            variance+=(x[d][1]-nodeMean)*(x[d][1]-nodeMean);
        if (topMax>=threshold_ && surfaceHigh_<0) { surfaceLow_=initialized_ ? lastTime_ : 0; surfaceHigh_=time; }
        if (domainMax>=threshold_ && globalHigh_<0) { globalLow_=initialized_ ? lastTime_ : 0; globalHigh_=time; }
        if (!initialized_) {
            initialAllSalt_=allSalt; initialAllWater_=allWater;
            for (std::size_t i=0;i<b.size();++i) {
                regions_[i].initialSalt=regions_[i].previousSalt=b[i].salt;
                regions_[i].initialWater=regions_[i].previousWater=b[i].water;
            }
        } else ++step_;
        // Outflow from the prescribed bottom node volumes into free domain is
        // minus the free-domain downward flux. Add their storage derivative.
        const double bottomSaltIn=initialized_ ? (dirSalt-previousDirSalt_)/dt-b.back().total : 0;
        const double bottomWaterIn=initialized_ ? (dirWater-previousDirWater_)/dt-b.back().waterOut : 0;
        cumulativeBottomSaltIn_+=dt*bottomSaltIn; cumulativeBottomWaterIn_+=dt*bottomWaterIn;
        cumulativeExternalSaltOut_+=dt*externalSalt; cumulativeExternalWaterOut_+=dt*externalWater;
        const double wholeSaltResidual=allSalt-initialAllSalt_+cumulativeExternalSaltOut_-cumulativeBottomSaltIn_;
        const double wholeWaterResidual=allWater-initialAllWater_+cumulativeExternalWaterOut_-cumulativeBottomWaterIn_;
        for (std::size_t i=0;i<b.size();++i) {
            auto& r=regions_[i]; const auto& v=b[i];
            r.cumulativeSalt+=dt*(v.total+v.otherSaltOut);
            r.cumulativeWater+=dt*(v.waterOut+v.otherWaterOut);
            r.cumulativeAdv+=dt*v.adv; r.cumulativeDiff+=dt*v.diff;
            // Current boundary law has zero salt on top/sides; fail rather than
            // silently mislabel cumulativeSalt as just cut flux in a new model.
            if (std::abs(v.otherSaltOut)>1e-30)
                throw std::runtime_error("Diagnostic schema requires zero top/side salt flux");
            budget_ << step_ << ',' << time << ',' << dt << ',' << r.name << ',' << r.requested << ',' << (r.cut-bottom_)/height_
                    << ',' << v.salt << ',' << r.initialSalt << ',' << v.water
                    << ',' << v.adv << ',' << v.diff << ',' << v.total << ',' << v.waterOut
                    << ',' << v.otherSaltOut << ',' << v.otherWaterOut << ',' << r.cumulativeSalt
                    << ',' << r.cumulativeAdv << ',' << r.cumulativeDiff << ',' << r.cumulativeWater
                    << ',' << v.salt-r.initialSalt+r.cumulativeSalt << ',' << v.water-r.initialWater+r.cumulativeWater
                    << ',' << v.salt-r.previousSalt+dt*(v.total+v.otherSaltOut)
                    << ',' << v.water-r.previousWater+dt*(v.waterOut+v.otherWaterOut)
                    << ',' << allSalt << ',' << allWater << ',' << bottomSaltIn << ',' << bottomWaterIn
                    << ',' << wholeSaltResidual << ',' << wholeWaterResidual << ',' << int(initialized_) << ',' << int(!thresholdReached()) << '\n';
            r.previousSalt=v.salt; r.previousWater=v.water;
        }
        surface_ << step_ << ',' << time << ',' << dt << ',' << topWeighted/topArea << ',' << topMax << ',' << domainMax
                 << ',' << domainMin << ',' << nodeMean << ',' << std::sqrt(variance/nodeCount) << ',' << topArea
                 << ',' << int(surfaceHigh_>=0) << ',' << int(globalHigh_>=0) << '\n';
        // Read-only weighted statistics; never enter the physical residual.
        const double mean=topWeighted/topArea;
        std::vector<double> rowSum(levelCount_,0.0), rowWeight(levelCount_,0.0);
        for (std::size_t d=0;d<x.size();++d) {
            rowSum[dofLevel_[d]]+=volumeWeight[d]*x[d][1];
            rowWeight[dofLevel_[d]]+=volumeWeight[d];
        }
        for (std::size_t l=0;l<levelCount_;++l) rowSum[l]/=rowWeight[l];
        double wvar=0, re=0, im=0, proj=0, den=0, rms=0, volume=0;
        std::array<double,2> hre={0,0}, him={0,0};
        const double k=2*std::acos(-1.0)*problem_.waveNumber()/width_;
        for (std::size_t d=0;d<x.size();++d) {
            const double phase=k*(dofX_[d]-gg_.bBoxMin()[0]);
            const double v=x[d][1]-mean;
            wvar+=topWeight[d]*v*v; re+=topWeight[d]*v*std::cos(phase); im-=topWeight[d]*v*std::sin(phase);
            for (std::size_t j=0;j<2;++j) {
                hre[j]+=topWeight[d]*v*std::cos((j+2)*phase);
                him[j]-=topWeight[d]*v*std::sin((j+2)*phase);
            }
            const double anomaly=x[d][1]-rowSum[dofLevel_[d]];
            rms+=volumeWeight[d]*anomaly*anomaly; volume+=volumeWeight[d];
            if (problem_.hasEigenprofile()) {
                const double basis=problem_.eigenAmplitudeAt(dofZ_[d])*std::cos(phase);
                proj+=volumeWeight[d]*anomaly*basis; den+=volumeWeight[d]*basis*basis;
            }
        }
        weighted_ << step_ << ',' << time << ',' << mean << ',' << std::sqrt(wvar/topArea)
                  << ',' << 2*std::hypot(re,im)/width_ << ',';
        if (den>0) weighted_ << proj/den; // Background projection is undefined, not zero.
        weighted_ << ',' << std::sqrt(rms/volume)
                  << ',' << 2*std::hypot(hre[0],him[0])/width_ << ',' << 2*std::hypot(hre[1],him[1])/width_
                  << ',' << int(den>0) << '\n';
        weighted_.flush(); checkStream(weighted_);
        if (writeLocalFluxes_) { localFlux_.flush(); checkStream(localFlux_); }
        budget_.flush(); surface_.flush(); checkStream(budget_); checkStream(surface_);
        previousDirSalt_=dirSalt; previousDirWater_=dirWater; lastTime_=time; initialized_=true;
        writeStop("in_progress_or_interrupted");
    }
};
} // namespace Dumux
#endif
