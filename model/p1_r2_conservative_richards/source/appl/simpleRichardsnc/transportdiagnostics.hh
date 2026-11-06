// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef P1_R2_CONSERVATIVE_TRANSPORT_DIAGNOSTICS_HH
#define P1_R2_CONSERVATIVE_TRANSPORT_DIAGNOSTICS_HH

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

// Serial 2D BOX diagnostics for total-liquid, NaCl, and reconstructed H2O
// mass. A unit out-of-plane thickness is used: inventories are kg and
// integrated face fluxes are kg/s.
template<class TypeTag>
class TransportDiagnostics
{
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using GV = GetPropType<TypeTag, Properties::GridVariables>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Solution = GetPropType<TypeTag, Properties::SolutionVector>;
    using Residual = GetPropType<TypeTag, Properties::LocalResidual>;

    struct Component {
        double total = 0.0;
        double salt = 0.0;
        double water() const { return total - salt; }
    };

    struct Cut {
        std::string name;
        double z = 0.0;
        Component flux;
        double saltAdvective = 0.0;
        double saltDiffusive = 0.0;
    };

    const Problem& problem_;
    const GG& gg_;
    GV& gv_;
    std::vector<double> dofZ_;
    std::vector<Cut> cuts_;
    std::ofstream budget_;
    std::ofstream cutFlux_;
    std::ofstream faceFlux_;
    std::ofstream surface_;
    double bottom_;
    double top_;
    double width_;
    double height_;
    double tolerance_;
    double threshold_;
    double lastTime_ = 0.0;
    Component initial_;
    Component previous_;
    Component previousBottom_;
    Component cumulativeExternalOut_;
    Component cumulativeBottomIn_;
    double surfaceLow_ = -1.0;
    double surfaceHigh_ = -1.0;
    double globalLow_ = -1.0;
    double globalHigh_ = -1.0;
    std::size_t step_ = 0;
    bool initialized_ = false;
    bool finalized_ = false;
    bool stopAtThreshold_;
    bool writeFaceFluxes_;

    static void checkStream(const std::ofstream& stream)
    {
        if (!stream)
            throw std::runtime_error("P1-R2 diagnostic output write failed");
    }

    static void numberOrNull(std::ostream& out, double value)
    {
        if (value < 0.0)
            out << "null";
        else
            out << value;
    }

    void writeStop(const std::string& reason) const
    {
        std::ofstream out("stop_record.json.tmp");
        out << std::setprecision(17)
            << "{\n  \"schema_version\": 2,\n"
            << "  \"reason\": \"" << reason << "\",\n"
            << "  \"last_accepted_time_s\": " << lastTime_ << ",\n"
            << "  \"accepted_steps\": " << step_ << ",\n"
            << "  \"initial_diagnostics_written\": "
            << (initialized_ ? "true" : "false") << ",\n"
            << "  \"requested_end_s\": " << getParam<double>("TimeLoop.TEnd")
            << ",\n  \"surface_threshold_low_s\": ";
        numberOrNull(out, surfaceLow_);
        out << ",\n  \"surface_threshold_high_s\": ";
        numberOrNull(out, surfaceHigh_);
        out << ",\n  \"global_threshold_low_s\": ";
        numberOrNull(out, globalLow_);
        out << ",\n  \"global_threshold_high_s\": ";
        numberOrNull(out, globalHigh_);
        out << "\n}\n";
        out.flush();
        checkStream(out);
        out.close();
        std::filesystem::rename("stop_record.json.tmp", "stop_record.json");
    }

public:
    TransportDiagnostics(const Problem& problem, const GG& gg, GV& gv)
        : problem_(problem)
        , gg_(gg)
        , gv_(gv)
        , bottom_(gg.bBoxMin()[1])
        , top_(gg.bBoxMax()[1])
        , width_(gg.bBoxMax()[0] - gg.bBoxMin()[0])
        , height_(top_ - bottom_)
        , tolerance_(1e-10 * std::max(1.0, height_))
        , threshold_(getParam<double>("Transport.SolubilityMoleFraction", 0.1))
        , stopAtThreshold_(getParam<bool>("Transport.StopAtSolubility", true))
        , writeFaceFluxes_(getParam<bool>("Transport.WriteFaceFluxes", false))
    {
        static_assert(
            GG::GridView::dimension == 2,
            "P1-R2 diagnostics currently support 2D BOX only");
        if (gg.gridView().comm().size() != 1)
            throw std::runtime_error("P1-R2 diagnostics require one MPI rank");
        if (!(width_ > 0.0 && height_ > 0.0
              && threshold_ > 0.0 && threshold_ < 1.0))
            throw std::runtime_error("Invalid geometry or threshold");

        for (const auto* name : {
                 "component_budget.csv",
                 "cut_fluxes.csv",
                 "independent_face_fluxes.csv",
                 "surface_events.csv",
                 "transport_metadata.json",
                 "stop_record.json",
             })
            if (std::filesystem::exists(name))
                throw std::runtime_error(
                    std::string("Refusing to overwrite diagnostic: ") + name);

        std::vector<double> levels;
        dofZ_.resize(gg.numDofs());
        for (const auto& vertex : vertices(gg.gridView()))
        {
            const auto dof = gg.vertexMapper().index(vertex);
            const auto z = vertex.geometry().center()[1];
            dofZ_[dof] = z;
            levels.push_back(z);
        }
        std::sort(levels.begin(), levels.end());
        levels.erase(
            std::unique(
                levels.begin(),
                levels.end(),
                [this](double a, double b) {
                    return std::abs(a - b) < tolerance_;
                }),
            levels.end());
        if (levels.size() < 3)
            throw std::runtime_error("Need at least three vertical vertex layers");

        // The first cut isolates the strongly prescribed bottom BOX control
        // volumes. Its local storage change and internal SCVF flux reconstruct
        // the Dirichlet reaction without using a whole-domain residual remainder.
        cuts_.push_back(
            {"bottom_cv_interface", 0.5 * (levels[0] + levels[1])});
        const auto requestedCuts = getParam<std::vector<double>>(
            "Transport.CutCoordinates", std::vector<double>{});
        for (std::size_t i = 0; i < requestedCuts.size(); ++i)
        {
            const auto target = requestedCuts[i];
            double best = 0.5 * (levels[0] + levels[1]);
            for (std::size_t level = 1; level + 1 < levels.size(); ++level)
            {
                const auto candidate = 0.5 * (levels[level] + levels[level + 1]);
                if (std::abs(candidate - target) < std::abs(best - target))
                    best = candidate;
            }
            if (std::abs(best - target) > 1e-12)
                throw std::runtime_error(
                    "Transport cut is not a conservative BOX face");
            cuts_.push_back({"requested_" + std::to_string(i), best});
        }

        budget_.open("component_budget.csv");
        budget_ << std::setprecision(17)
                << "step,time_s,dt_s,total_mass_kg,salt_mass_kg,water_mass_kg,"
                << "external_total_out_kg_s,external_salt_out_kg_s,"
                << "external_water_out_kg_s,bottom_total_in_kg_s,"
                << "bottom_salt_in_kg_s,bottom_water_in_kg_s,"
                << "cumulative_external_total_out_kg,"
                << "cumulative_external_salt_out_kg,"
                << "cumulative_bottom_total_in_kg,"
                << "cumulative_bottom_salt_in_kg,"
                << "total_residual_kg,salt_residual_kg,water_residual_kg,"
                << "step_total_residual_kg,step_salt_residual_kg,"
                << "step_water_residual_kg,bottom_reaction_defined\n";

        cutFlux_.open("cut_fluxes.csv");
        cutFlux_ << std::setprecision(17)
                 << "step,time_s,dt_s,cut,cut_z_m,total_down_kg_s,"
                 << "salt_down_kg_s,water_down_kg_s,"
                 << "salt_advective_down_kg_s,salt_diffusive_down_kg_s\n";

        surface_.open("surface_events.csv");
        surface_ << std::setprecision(17)
                 << "step,time_s,dt_s,x_top_area_mean,x_top_max,"
                 << "x_domain_max,x_domain_min,top_area_m,"
                 << "top_total_out_kg_s,top_salt_out_kg_s,"
                 << "top_water_out_kg_s,surface_threshold_reached,"
                 << "global_threshold_reached\n";

        if (writeFaceFluxes_)
        {
            faceFlux_.open("independent_face_fluxes.csv");
            faceFlux_ << std::setprecision(17)
                      << "step,time_s,kind,cut,x_m,z_m,area_m2,"
                      << "total_flux_kg_s,salt_flux_kg_s,water_flux_kg_s,"
                      << "salt_advective_kg_s,salt_diffusive_kg_s,"
                      << "upwind_x,upwind_w,face_density_kg_m3\n";
            checkStream(faceFlux_);
        }

        std::ofstream meta("transport_metadata.json");
        meta << std::setprecision(17)
             << "{\n"
             << "  \"schema_version\": 3,\n"
             << "  \"equation_0\": \"total liquid mass\",\n"
             << "  \"equation_1\": \"dissolved NaCl mass\",\n"
             << "  \"inventory_unit\": \"kg per unit out-of-plane thickness\",\n"
             << "  \"integrated_flux_unit\": \"kg/s per unit out-of-plane thickness\",\n"
             << "  \"diffusion_reference\": \"mass averaged\",\n"
             << "  \"darcy_velocity_assumption\": \"mixture mass-average velocity approximation\",\n"
             << "  \"advection_face_rule\": \"upwind rho and mass fraction\",\n"
             << "  \"diffusion_face_rule\": \"arithmetic rho and D, BOX gradient of mass fraction\",\n"
             << "  \"bottom_reaction\": \"bottom control-volume storage change plus independently summed internal SCVF flux; not a whole-domain residual remainder\",\n"
             << "  \"time_integration\": \"accepted-step backward Euler endpoint flux\",\n"
             << "  \"evaporation_boundary_mode\": \""
             << problem_.evaporationBoundaryMode() << "\",\n"
             << "  \"bottom_boundary_mode\": \""
             << (problem_.bottomIsDirichlet() ? "reservoir" : "no_flux")
             << "\",\n"
             << "  \"face_flux_output\": "
             << (writeFaceFluxes_ ? "true" : "false") << "\n"
             << "}\n";
        meta.flush();
        checkStream(meta);
        checkStream(budget_);
        checkStream(cutFlux_);
        checkStream(surface_);
        writeStop("initialized_not_completed");
    }

    ~TransportDiagnostics()
    {
        if (!finalized_)
        {
            try
            {
                writeStop("exception_or_incomplete");
            }
            catch (...)
            {}
        }
    }

    bool thresholdReached() const
    {
        return globalHigh_ >= 0.0 || surfaceHigh_ >= 0.0;
    }

    bool shouldStop() const
    {
        return stopAtThreshold_ && thresholdReached();
    }

    std::string thresholdReason() const
    {
        return surfaceHigh_ >= 0.0
                   ? "surface_solubility_threshold"
                   : "interior_solubility_threshold";
    }

    void finish(const std::string& reason)
    {
        writeStop(reason);
        finalized_ = true;
    }

    void write(double time, double dt, const Solution& x)
    {
        if (!std::isfinite(time) || !std::isfinite(dt)
            || (initialized_
                && (!(dt > 0.0)
                    || std::abs(time - lastTime_ - dt)
                           > 1e-9 * std::max(1.0, time)))
            || (!initialized_ && (time != 0.0 || dt != 0.0)))
            throw std::runtime_error(
                "Invalid accepted-step time for P1-R2 budget");

        if (initialized_)
            ++step_;
        const auto outputStep = step_;
        Component inventory;
        Component bottomInventory;
        Component externalOut;
        Component topOut;
        std::vector<Cut> currentCuts = cuts_;
        for (auto& cut : currentCuts)
        {
            cut.flux = {};
            cut.saltAdvective = 0.0;
            cut.saltDiffusive = 0.0;
        }

        double topWeighted = 0.0;
        double topArea = 0.0;
        double topMax = -std::numeric_limits<double>::infinity();
        double domainMax = topMax;
        double domainMin = std::numeric_limits<double>::infinity();

        for (std::size_t dof = 0; dof < x.size(); ++dof)
        {
            if (!std::isfinite(x[dof][0]) || !std::isfinite(x[dof][1]))
                throw std::runtime_error("Non-finite accepted state");
            domainMax = std::max(domainMax, static_cast<double>(x[dof][1]));
            domainMin = std::min(domainMin, static_cast<double>(x[dof][1]));
            if (std::abs(dofZ_[dof] - top_) < tolerance_)
                topMax = std::max(topMax, static_cast<double>(x[dof][1]));
        }

        for (const auto& element : elements(gg_.gridView()))
        {
            auto geometry = localView(gg_);
            geometry.bind(element);
            auto vv = localView(gv_.curGridVolVars());
            vv.bind(element, geometry, x);
            auto cache = localView(gv_.gridFluxVarsCache());
            cache.bind(element, geometry, vv);

            for (const auto& scv : scvs(geometry))
            {
                const auto storage = Residual::massStorage(problem_, vv[scv]);
                const Component contribution{
                    storage[0] * scv.volume(),
                    storage[1] * scv.volume(),
                };
                if (!std::isfinite(contribution.total)
                    || !std::isfinite(contribution.salt))
                    throw std::runtime_error("Non-finite mass inventory");
                inventory.total += contribution.total;
                inventory.salt += contribution.salt;
                if (std::abs(dofZ_[scv.dofIndex()] - bottom_) < tolerance_)
                {
                    bottomInventory.total += contribution.total;
                    bottomInventory.salt += contribution.salt;
                }
            }

            for (const auto& face : scvfs(geometry))
            {
                const auto& inside = geometry.scv(face.insideScvIdx());
                const auto zi = dofZ_[inside.dofIndex()];
                if (face.boundary())
                {
                    const bool isBottom =
                        std::abs(face.ipGlobal()[1] - bottom_) < tolerance_;
                    if (isBottom && problem_.bottomIsDirichlet())
                        continue;

                    const auto flux =
                        problem_.neumann(element, geometry, vv, cache, face);
                    const Component faceOut{
                        flux[0] * face.area(),
                        flux[1] * face.area(),
                    };
                    externalOut.total += faceOut.total;
                    externalOut.salt += faceOut.salt;
                    const bool isTop =
                        std::abs(face.ipGlobal()[1] - top_) < tolerance_;
                    if (isTop)
                    {
                        topOut.total += faceOut.total;
                        topOut.salt += faceOut.salt;
                        topArea += face.area();
                        topWeighted +=
                            face.area() * x[inside.dofIndex()][1];
                    }
                    if (writeFaceFluxes_)
                        faceFlux_
                            << outputStep << ',' << time << ",boundary,,"
                            << face.ipGlobal()[0] << ','
                            << face.ipGlobal()[1] << ',' << face.area() << ','
                            << faceOut.total << ',' << faceOut.salt << ','
                            << faceOut.water() << ",nan,nan,nan,nan,nan\n";
                    continue;
                }

                const auto zo =
                    dofZ_[geometry.scv(face.outsideScvIdx()).dofIndex()];
                const auto flux = Residual::transportFlux(
                    problem_, element, geometry, vv, face, cache);
                for (std::size_t index = 0; index < currentCuts.size(); ++index)
                {
                    auto& cut = currentCuts[index];
                    const int sign =
                        int(zi > cut.z) - int(zo > cut.z);
                    if (!sign)
                        continue;
                    cut.flux.total += sign * flux.totalMass;
                    cut.flux.salt += sign * flux.saltTotal;
                    cut.saltAdvective += sign * flux.saltAdvective;
                    cut.saltDiffusive += sign * flux.saltDiffusive;
                    if (writeFaceFluxes_)
                        faceFlux_
                            << outputStep << ',' << time << ",internal_cut,"
                            << cut.name << ',' << face.ipGlobal()[0] << ','
                            << face.ipGlobal()[1] << ',' << face.area() << ','
                            << sign * flux.totalMass << ','
                            << sign * flux.saltTotal << ','
                            << sign * flux.waterMass << ','
                            << sign * flux.saltAdvective << ','
                            << sign * flux.saltDiffusive << ','
                            << flux.upwindMoleFraction << ','
                            << flux.upwindMassFraction << ','
                            << flux.faceDensity << '\n';
                }
            }
        }

        if (!(topArea > 0.0)
            || std::abs(topArea - width_) > 1e-9 * width_)
            throw std::runtime_error(
                "Top boundary quadrature does not cover the width");

        Component bottomIn;
        const bool bottomReactionDefined =
            initialized_ && problem_.bottomIsDirichlet();
        if (bottomReactionDefined)
        {
            bottomIn.total =
                (bottomInventory.total - previousBottom_.total) / dt
                - currentCuts.front().flux.total;
            bottomIn.salt =
                (bottomInventory.salt - previousBottom_.salt) / dt
                - currentCuts.front().flux.salt;
        }

        if (!initialized_)
        {
            initial_ = inventory;
            previous_ = inventory;
            previousBottom_ = bottomInventory;
        }
        else
        {
            cumulativeExternalOut_.total += dt * externalOut.total;
            cumulativeExternalOut_.salt += dt * externalOut.salt;
            cumulativeBottomIn_.total += dt * bottomIn.total;
            cumulativeBottomIn_.salt += dt * bottomIn.salt;
        }

        const Component cumulativeResidual{
            inventory.total - initial_.total + cumulativeExternalOut_.total
                - cumulativeBottomIn_.total,
            inventory.salt - initial_.salt + cumulativeExternalOut_.salt
                - cumulativeBottomIn_.salt,
        };
        const Component stepResidual{
            initialized_
                ? inventory.total - previous_.total
                      + dt * externalOut.total - dt * bottomIn.total
                : 0.0,
            initialized_
                ? inventory.salt - previous_.salt
                      + dt * externalOut.salt - dt * bottomIn.salt
                : 0.0,
        };

        budget_ << outputStep << ',' << time << ',' << dt << ','
                << inventory.total << ',' << inventory.salt << ','
                << inventory.water() << ',' << externalOut.total << ','
                << externalOut.salt << ',' << externalOut.water() << ','
                << bottomIn.total << ',' << bottomIn.salt << ','
                << bottomIn.water() << ','
                << cumulativeExternalOut_.total << ','
                << cumulativeExternalOut_.salt << ','
                << cumulativeBottomIn_.total << ','
                << cumulativeBottomIn_.salt << ','
                << cumulativeResidual.total << ','
                << cumulativeResidual.salt << ','
                << cumulativeResidual.water() << ','
                << stepResidual.total << ',' << stepResidual.salt << ','
                << stepResidual.water() << ','
                << int(bottomReactionDefined) << '\n';

        for (const auto& cut : currentCuts)
            cutFlux_ << outputStep << ',' << time << ',' << dt << ','
                     << cut.name << ',' << cut.z << ','
                     << cut.flux.total << ',' << cut.flux.salt << ','
                     << cut.flux.water() << ',' << cut.saltAdvective << ','
                     << cut.saltDiffusive << '\n';

        if (writeFaceFluxes_ && bottomReactionDefined)
            faceFlux_ << outputStep << ',' << time
                      << ",bottom_reaction_local_balance,bottom,,,,"
                      << bottomIn.total << ',' << bottomIn.salt << ','
                      << bottomIn.water()
                      << ",nan,nan,nan,nan,nan\n";

        if (topMax >= threshold_ && surfaceHigh_ < 0.0)
        {
            surfaceLow_ = initialized_ ? lastTime_ : 0.0;
            surfaceHigh_ = time;
        }
        if (domainMax >= threshold_ && globalHigh_ < 0.0)
        {
            globalLow_ = initialized_ ? lastTime_ : 0.0;
            globalHigh_ = time;
        }
        surface_ << outputStep << ',' << time << ',' << dt << ','
                 << topWeighted / topArea << ',' << topMax << ','
                 << domainMax << ',' << domainMin << ',' << topArea << ','
                 << topOut.total << ',' << topOut.salt << ','
                 << topOut.water() << ',' << int(surfaceHigh_ >= 0.0) << ','
                 << int(globalHigh_ >= 0.0) << '\n';

        budget_.flush();
        cutFlux_.flush();
        surface_.flush();
        checkStream(budget_);
        checkStream(cutFlux_);
        checkStream(surface_);
        if (writeFaceFluxes_)
        {
            faceFlux_.flush();
            checkStream(faceFlux_);
        }

        previous_ = inventory;
        previousBottom_ = bottomInventory;
        lastTime_ = time;
        initialized_ = true;
        writeStop("in_progress_or_interrupted");
    }
};

} // namespace Dumux

#endif
