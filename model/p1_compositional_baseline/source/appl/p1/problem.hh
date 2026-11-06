#ifndef P1_COMPOSITIONAL_PROBLEM_HH
#define P1_COMPOSITIONAL_PROBLEM_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<class TypeTag>
class P1Problem : public PorousMediumFlowProblem<TypeTag>
{
    using Parent = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int verticalIdx = dimWorld - 1;
    static constexpr int waterIdx = FluidSystem::H2OIdx;
    static constexpr int airIdx = FluidSystem::AirIdx;
    static constexpr int saltIdx = FluidSystem::NaClIdx;
    static constexpr int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr Scalar pi = 3.1415926535897932384626433832795;

public:
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    struct Audit
    {
        std::array<double, 3> inventory{};
        std::array<double, 3> boundaryOut{};
        double topSaltOut = 0.0;
        double topArea = 0.0;
        double topSaltWeighted = 0.0;
        double topDensityWeighted = 0.0;
        double topSwWeighted = 0.0;
        double topSaltCos = 0.0;
        double topSaltSin = 0.0;
        double topCos = 0.0;
        double topSin = 0.0;
        double minSw = 1.0;
        double maxSw = 0.0;
        double minSalt = std::numeric_limits<double>::max();
        double maxSalt = 0.0;
        double minPressure = std::numeric_limits<double>::max();
        double maxPressure = -std::numeric_limits<double>::max();
        std::array<double, 2> maxDarcyVelocity{};
        std::array<std::size_t, 4> phaseStateCount{};
    };

    explicit P1Problem(std::shared_ptr<const GridGeometry> gridGeometry)
    : Parent(gridGeometry)
    , name_(getParam<std::string>("Problem.Name", "p1"))
    , boundaryMode_(getParam<std::string>("Problem.BoundaryMode", "Closed"))
    , initialPhaseState_(getParam<std::string>("Problem.InitialPhaseState", "Both"))
    , saltProfile_(getParam<std::string>("Problem.InitialSaltProfile", "Constant"))
    , atmosphericPressure_(getParam<Scalar>("Problem.AtmosphericPressure", 1.0e5))
    , hydrostaticDensity_(getParam<Scalar>("Problem.HydrostaticDensity", 1006.0241044335427))
    , hydrostaticGasDensity_(getParam<Scalar>("Problem.HydrostaticGasDensity", 1.177716))
    , hydrostaticGasDensityGradient_(getParam<Scalar>("Problem.HydrostaticGasDensityGradient", 1.2e-4))
    , gravity_(getParam<Scalar>("Problem.GravitationalAcceleration", 9.8))
    , waterLevelDepth_(getParam<Scalar>("Problem.WaterLevelDepth", 0.6307756332876138))
    , initialSalt_(getParam<Scalar>("Problem.SaltMoleFraction", 0.0036))
    , initialAir_(getParam<Scalar>("Problem.InitialDissolvedAirMoleFraction", 1.0e-8))
    , initialGasSaturation_(getParam<Scalar>("Problem.InitialGasSaturation", -1.0))
    , perturbationAmplitude_(getParam<Scalar>("Problem.PerturbationAmplitude", 0.0))
    , waveNumber_(getParam<int>("Problem.WaveNumber", 1))
    , evaporationRate_(getParam<Scalar>("Problem.EvaporationRate", 0.0))
    , evaporationWaterDensity_(getParam<Scalar>("Problem.EvaporationWaterDensity", 998.203))
    , testComponent_(getParam<int>("Problem.TestComponent", waterIdx))
    , testMolarFlux_(getParam<Scalar>("Problem.TestMolarFlux", 0.0))
    {
        FluidSystem::init();
        if (boundaryMode_ != "Closed" && boundaryMode_ != "Flux" && boundaryMode_ != "Reservoir")
            DUNE_THROW(Dune::InvalidStateException, "Problem.BoundaryMode must be Closed, Flux, or Reservoir");
        if (initialPhaseState_ != "Both" && initialPhaseState_ != "Liquid")
            DUNE_THROW(Dune::InvalidStateException, "Problem.InitialPhaseState must be Both or Liquid");
        if (initialSalt_ <= 0.0 || initialSalt_ >= 0.1)
            DUNE_THROW(Dune::InvalidStateException, "P1 initial NaCl mole fraction must be in (0,0.1)");
        if (testComponent_ < 0 || testComponent_ >= 3)
            DUNE_THROW(Dune::InvalidStateException, "Problem.TestComponent must be 0, 1, or 2");
        if (saltProfile_ == "Eigen")
            readEigenprofile_(getParam<std::string>("Problem.EigenprofilePath"));

        projectedArea_ = 1.0;
        if constexpr (dimWorld > 1)
            for (int d = 0; d < verticalIdx; ++d)
                projectedArea_ *= this->gridGeometry().bBoxMax()[d] - this->gridGeometry().bBoxMin()[d];
        if (!(projectedArea_ > 0.0))
            DUNE_THROW(Dune::InvalidStateException, "Nonpositive projected area");
    }

    const std::string& name() const
    { return name_; }

    void setTimeState(Scalar time, Scalar dt)
    {
        time_ = time;
        dt_ = dt;
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes types;
        if (boundaryMode_ == "Reservoir" && onBottom_(globalPos))
            types.setAllDirichlet();
        else
            types.setAllNeumann();
        return types;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return stateAt_(globalPos, false); }

    NumEqVector neumann(const Element&,
                        const FVElementGeometry&,
                        const ElementVolumeVariables&,
                        const ElementFluxVariablesCache&,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        if (!onTop_(scvf.ipGlobal()))
            return flux;

        if (boundaryMode_ == "Flux")
            flux[testComponent_] = testMolarFlux_;
        else if (boundaryMode_ == "Reservoir")
            flux[waterIdx] = evaporationRate_*evaporationWaterDensity_/FluidSystem::molarMass(waterIdx);

        return flux;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return stateAt_(globalPos, true); }

    template<class GridVariables>
    Audit audit(const SolutionVector& x, const GridVariables& gridVariables) const
    {
        Audit audit;
        LocalResidual localResidual(this);

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bind(element);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bind(element, fvGeometry, x);
            auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                const double volume = scv.volume()*volVars.extrusionFactor();
                const int state = volVars.priVars().state();
                if (state >= 0 && state < static_cast<int>(audit.phaseStateCount.size()))
                    ++audit.phaseStateCount[state];

                for (int phaseIdx = 0; phaseIdx < 2; ++phaseIdx)
                {
                    double fractionSum = 0.0;
                    const bool activePhase = volVars.saturation(phaseIdx) > 1e-12;
                    for (int compIdx = 0; compIdx < 3; ++compIdx)
                    {
                        const double fraction = volVars.moleFraction(phaseIdx, compIdx);
                        if (activePhase && (!std::isfinite(fraction)
                            || fraction < -1e-10 || fraction > 1.0 + 1e-10))
                            DUNE_THROW(Dune::InvalidStateException, "Nonphysical component mole fraction");
                        fractionSum += fraction;
                        audit.inventory[compIdx] += volume*volVars.porosity()
                            *volVars.saturation(phaseIdx)*volVars.molarDensity(phaseIdx)*fraction;
                    }
                    if (activePhase && std::abs(fractionSum - 1.0) > 1e-8)
                        DUNE_THROW(Dune::InvalidStateException, "Active-phase mole fractions do not sum to one");
                }

                const double sw = volVars.saturation(liquidPhaseIdx);
                const double xs = volVars.moleFraction(liquidPhaseIdx, saltIdx);
                audit.minSw = std::min(audit.minSw, sw);
                audit.maxSw = std::max(audit.maxSw, sw);
                audit.minSalt = std::min(audit.minSalt, xs);
                audit.maxSalt = std::max(audit.maxSalt, xs);
                audit.minPressure = std::min(audit.minPressure, double(volVars.pressure(liquidPhaseIdx)));
                audit.maxPressure = std::max(audit.maxPressure, double(volVars.pressure(liquidPhaseIdx)));
                if (sw < -1e-10 || sw > 1.0 + 1e-10 || !std::isfinite(volVars.pressure(liquidPhaseIdx)))
                    DUNE_THROW(Dune::InvalidStateException, "Nonphysical saturation or pressure");
            }

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (!scvf.boundary())
                    continue;

                std::array<double, 3> faceFlux{};
                if (boundaryMode_ == "Reservoir" && onBottom_(scvf.ipGlobal()))
                {
                    const auto flux = localResidual.computeFlux(*this, element, fvGeometry,
                                                               elemVolVars, scvf, elemFluxVarsCache);
                    for (int compIdx = 0; compIdx < 3; ++compIdx)
                        faceFlux[compIdx] = flux[compIdx];
                }
                else
                {
                    const auto fluxDensity = neumann(element, fvGeometry, elemVolVars,
                                                     elemFluxVarsCache, scvf);
                    const double factor = scvf.area()*elemVolVars[scvf.insideScvIdx()].extrusionFactor();
                    for (int compIdx = 0; compIdx < 3; ++compIdx)
                        faceFlux[compIdx] = fluxDensity[compIdx]*factor;
                }

                for (int compIdx = 0; compIdx < 3; ++compIdx)
                    audit.boundaryOut[compIdx] += faceFlux[compIdx];

                if (onTop_(scvf.ipGlobal()))
                {
                    const auto& volVars = elemVolVars[scvf.insideScvIdx()];
                    const double area = scvf.area()*volVars.extrusionFactor();
                    const double xs = volVars.moleFraction(liquidPhaseIdx, saltIdx);
                    audit.topSaltOut += faceFlux[saltIdx];
                    audit.topArea += area;
                    audit.topSaltWeighted += area*xs;
                    audit.topDensityWeighted += area*volVars.density(liquidPhaseIdx);
                    audit.topSwWeighted += area*volVars.saturation(liquidPhaseIdx);
                    if constexpr (dimWorld > 1)
                    {
                        const double width = this->gridGeometry().bBoxMax()[0]
                                           - this->gridGeometry().bBoxMin()[0];
                        const double angle = 2.0*pi*waveNumber_
                                           *(scvf.ipGlobal()[0] - this->gridGeometry().bBoxMin()[0])/width;
                        audit.topSaltCos += area*xs*std::cos(angle);
                        audit.topSaltSin += area*xs*std::sin(angle);
                        audit.topCos += area*std::cos(angle);
                        audit.topSin += area*std::sin(angle);
                    }
                }
            }
        }

        using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
        VelocityOutput velocityOutput(gridVariables);
        if (velocityOutput.enableOutput())
        {
            typename VelocityOutput::VelocityVector velocity(this->gridGeometry().gridView().size(0));
            for (int phaseIdx = 0; phaseIdx < 2; ++phaseIdx)
            {
                for (auto& value : velocity) value = 0.0;
                for (const auto& element : elements(this->gridGeometry().gridView()))
                {
                    auto fvGeometry = localView(this->gridGeometry());
                    fvGeometry.bind(element);
                    auto elemVolVars = localView(gridVariables.curGridVolVars());
                    elemVolVars.bind(element, fvGeometry, x);
                    auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
                    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
                    velocityOutput.calculateVelocity(velocity, element, fvGeometry,
                                                     elemVolVars, elemFluxVarsCache, phaseIdx);
                }
                for (const auto& value : velocity)
                    audit.maxDarcyVelocity[phaseIdx] = std::max(audit.maxDarcyVelocity[phaseIdx],
                                                                double(value.two_norm()));
            }
        }
        return audit;
    }

    template<class GridVariables>
    void writeCellState(const std::string& path,
                        const SolutionVector& x,
                        const GridVariables& gridVariables) const
    {
        if (!getParam<bool>("Problem.WriteCellStates", false))
            return;

        std::ofstream output(path);
        output << std::setprecision(17)
               << "cell,coordinate_0_m,coordinate_1_m,bulk_volume_m3,porosity,phase_state,";
        for (int phaseIdx = 0; phaseIdx < 2; ++phaseIdx)
            output << "s" << phaseIdx << ",molar_density" << phaseIdx << "_mol_m3,"
                   << "x" << phaseIdx << "_water,x" << phaseIdx << "_air,x" << phaseIdx << "_salt,";
        output << "\n";

        std::size_t cell = 0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bind(element);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bind(element, fvGeometry, x);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                const auto center = scv.center();
                output << cell++ << ',' << center[0] << ',';
                if constexpr (dimWorld > 1) output << center[1];
                else output << center[0];
                output << ',' << scv.volume()*volVars.extrusionFactor()
                       << ',' << volVars.porosity() << ',' << volVars.priVars().state() << ',';
                for (int phaseIdx = 0; phaseIdx < 2; ++phaseIdx)
                {
                    output << volVars.saturation(phaseIdx) << ',' << volVars.molarDensity(phaseIdx) << ',';
                    for (int compIdx = 0; compIdx < 3; ++compIdx)
                        output << volVars.moleFraction(phaseIdx, compIdx) << ',';
                }
                output << '\n';
            }
        }
    }

    void initializeBudget(const Audit& audit)
    {
        initialInventory_ = audit.inventory;
        cumulativeOut_ = {};
        writeMetadata_(audit);
        record(audit, 0.0, 0.0);
    }

    void record(const Audit& audit, Scalar acceptedTime, Scalar acceptedDt)
    {
        for (int compIdx = 0; compIdx < 3; ++compIdx)
            cumulativeOut_[compIdx] += acceptedDt*audit.boundaryOut[compIdx];

        if (!budget_.is_open())
        {
            budget_.open("component_budget.csv");
            budget_ << std::setprecision(17)
                    << "time_s,dt_s,water_inventory_mol_m2,air_inventory_mol_m2,salt_inventory_mol_m2,"
                    << "water_out_mol_m2_s,air_out_mol_m2_s,salt_out_mol_m2_s,"
                    << "water_cumulative_out_mol_m2,air_cumulative_out_mol_m2,salt_cumulative_out_mol_m2,"
                    << "water_residual_mol_m2,air_residual_mol_m2,salt_residual_mol_m2,"
                    << "water_signal_relative_residual,air_signal_relative_residual,salt_signal_relative_residual,"
                    << "top_salt_out_mol_m2_s,top_salt_mean,top_liquid_density_kg_m3,top_sw,mode_amplitude,"
                    << "min_sw,max_sw,min_salt_x,max_salt_x,min_liquid_pressure_pa,max_liquid_pressure_pa,"
                    << "max_liquid_darcy_velocity_m_s,max_gas_darcy_velocity_m_s,"
                    << "state_first_only,state_both,state_second_only\n";
        }

        std::array<double, 3> residual{};
        std::array<double, 3> relative{};
        for (int compIdx = 0; compIdx < 3; ++compIdx)
        {
            residual[compIdx] = (audit.inventory[compIdx] - initialInventory_[compIdx]
                               + cumulativeOut_[compIdx])/projectedArea_;
            const double signal = std::max({std::abs(audit.inventory[compIdx] - initialInventory_[compIdx]),
                                            std::abs(cumulativeOut_[compIdx]), 1e-30})/projectedArea_;
            relative[compIdx] = std::abs(residual[compIdx])/signal;
        }

        const double topArea = std::max(audit.topArea, 1e-30);
        const double topMean = audit.topSaltWeighted/topArea;
        double modeAmplitude = 0.0;
        if constexpr (dimWorld > 1)
        {
            const double cosine = (audit.topSaltCos - topMean*audit.topCos)/topArea;
            const double sine = (audit.topSaltSin - topMean*audit.topSin)/topArea;
            modeAmplitude = 2.0*std::sqrt(cosine*cosine + sine*sine);
        }

        budget_ << acceptedTime << ',' << acceptedDt;
        for (const auto value : audit.inventory) budget_ << ',' << value/projectedArea_;
        for (const auto value : audit.boundaryOut) budget_ << ',' << value/projectedArea_;
        for (const auto value : cumulativeOut_) budget_ << ',' << value/projectedArea_;
        for (const auto value : residual) budget_ << ',' << value;
        for (const auto value : relative) budget_ << ',' << value;
        budget_ << ',' << audit.topSaltOut/projectedArea_
                << ',' << topMean
                << ',' << audit.topDensityWeighted/topArea
                << ',' << audit.topSwWeighted/topArea
                << ',' << modeAmplitude
                << ',' << audit.minSw << ',' << audit.maxSw
                << ',' << audit.minSalt << ',' << audit.maxSalt
                << ',' << audit.minPressure << ',' << audit.maxPressure
                << ',' << audit.maxDarcyVelocity[liquidPhaseIdx]
                << ',' << audit.maxDarcyVelocity[gasPhaseIdx]
                << ',' << audit.phaseStateCount[Indices::firstPhaseOnly]
                << ',' << audit.phaseStateCount[Indices::bothPhases]
                << ',' << audit.phaseStateCount[Indices::secondPhaseOnly] << '\n';
        budget_.flush();
    }

private:
    PrimaryVariables stateAt_(const GlobalPosition& globalPos, bool perturb) const
    {
        PrimaryVariables values(0.0);
        const Scalar z = globalPos[verticalIdx] - this->gridGeometry().bBoxMin()[verticalIdx];
        const Scalar liquidPressure = atmosphericPressure_
            - hydrostaticDensity_*gravity_*(waterLevelDepth_ + z);
        const Scalar height = this->gridGeometry().bBoxMax()[verticalIdx]
                            - this->gridGeometry().bBoxMin()[verticalIdx];
        const Scalar gasDepth = height - z;
        const Scalar gasPressure = atmosphericPressure_ + gravity_*
            (hydrostaticGasDensity_*gasDepth
             + 0.5*hydrostaticGasDensityGradient_*gasDepth*gasDepth);
        const Scalar pc = gasPressure - liquidPressure;
        const Scalar sw = this->spatialParams().fluidMatrixInteractionAtPos(globalPos).sw(pc);

        Scalar salt = initialSalt_;
        if (perturb)
        {
            if (saltProfile_ == "CosineVertical")
            {
                const Scalar height = this->gridGeometry().bBoxMax()[verticalIdx]
                                    - this->gridGeometry().bBoxMin()[verticalIdx];
                salt += perturbationAmplitude_*std::cos(pi*z/height);
            }
            else if (saltProfile_ == "Eigen")
            {
                Scalar horizontal = 1.0;
                if constexpr (dimWorld > 1)
                {
                    const Scalar width = this->gridGeometry().bBoxMax()[0]
                                       - this->gridGeometry().bBoxMin()[0];
                    horizontal = std::cos(2.0*pi*waveNumber_
                                         *(globalPos[0] - this->gridGeometry().bBoxMin()[0])/width);
                }
                salt += perturbationAmplitude_*interpolateEigen_(z)*horizontal;
            }
        }
        if (salt <= 0.0 || salt >= 0.1)
            DUNE_THROW(Dune::InvalidStateException, "Initial salt perturbation leaves P1 admissible range");

        values[Indices::pressureIdx] = liquidPressure;
        values[saltIdx] = salt;
        if (initialPhaseState_ == "Liquid")
        {
            values.setState(Indices::firstPhaseOnly);
            values[Indices::switchIdx] = initialAir_;
        }
        else
        {
            values.setState(Indices::bothPhases);
            values[Indices::switchIdx] = initialGasSaturation_ >= 0.0
                ? initialGasSaturation_ : 1.0 - sw;
        }
        return values;
    }

    bool onTop_(const GlobalPosition& globalPos) const
    {
        return globalPos[verticalIdx] > this->gridGeometry().bBoxMax()[verticalIdx] - 1e-9;
    }

    bool onBottom_(const GlobalPosition& globalPos) const
    {
        return globalPos[verticalIdx] < this->gridGeometry().bBoxMin()[verticalIdx] + 1e-9;
    }

    void readEigenprofile_(const std::string& path)
    {
        std::ifstream input(path);
        if (!input)
            DUNE_THROW(Dune::IOError, "Cannot open eigenprofile " + path);
        std::string line;
        while (std::getline(input, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream stream(line);
            Scalar z = 0.0;
            Scalar amplitude = 0.0;
            if (stream >> z >> amplitude)
            {
                eigenZ_.push_back(z);
                eigenValues_.push_back(amplitude);
            }
        }
        if (eigenZ_.size() < 2 || !std::is_sorted(eigenZ_.begin(), eigenZ_.end()))
            DUNE_THROW(Dune::InvalidStateException, "Eigenprofile must contain sorted z/amplitude pairs");
    }

    Scalar interpolateEigen_(Scalar z) const
    {
        if (z <= eigenZ_.front()) return eigenValues_.front();
        if (z >= eigenZ_.back()) return eigenValues_.back();
        const auto upper = std::lower_bound(eigenZ_.begin(), eigenZ_.end(), z);
        const std::size_t i = static_cast<std::size_t>(upper - eigenZ_.begin());
        const Scalar weight = (z - eigenZ_[i - 1])/(eigenZ_[i] - eigenZ_[i - 1]);
        return eigenValues_[i - 1] + weight*(eigenValues_[i] - eigenValues_[i - 1]);
    }

    void writeMetadata_(const Audit& initialAudit) const
    {
        std::ofstream metadata("model_metadata.json");
        metadata << std::setprecision(17)
                 << "{\n"
                 << "  \"model\": \"isothermal_two_phase_three_component_2pnc\",\n"
                 << "  \"use_moles\": true,\n"
                 << "  \"precipitation_reaction\": false,\n"
                 << "  \"diffusion_reference_system\": \"mass_averaged_fick_converted_to_molar_component_flux\",\n"
                 << "  \"inventory_unit\": \"mol/m2_projected\",\n"
                 << "  \"flux_unit\": \"mol/(m2_projected*s)\",\n"
                 << "  \"projected_area_m2_per_extrusion\": " << projectedArea_ << ",\n"
                 << "  \"boundary_mode\": \"" << boundaryMode_ << "\",\n"
                 << "  \"initial_water_mol_m2\": " << initialAudit.inventory[waterIdx]/projectedArea_ << ",\n"
                 << "  \"initial_air_mol_m2\": " << initialAudit.inventory[airIdx]/projectedArea_ << ",\n"
                 << "  \"initial_salt_mol_m2\": " << initialAudit.inventory[saltIdx]/projectedArea_ << "\n"
                 << "}\n";
    }

    std::string name_;
    std::string boundaryMode_;
    std::string initialPhaseState_;
    std::string saltProfile_;
    Scalar atmosphericPressure_;
    Scalar hydrostaticDensity_;
    Scalar hydrostaticGasDensity_;
    Scalar hydrostaticGasDensityGradient_;
    Scalar gravity_;
    Scalar waterLevelDepth_;
    Scalar initialSalt_;
    Scalar initialAir_;
    Scalar initialGasSaturation_;
    Scalar perturbationAmplitude_;
    int waveNumber_;
    Scalar evaporationRate_;
    Scalar evaporationWaterDensity_;
    int testComponent_;
    Scalar testMolarFlux_;
    Scalar projectedArea_ = 1.0;
    Scalar time_ = 0.0;
    Scalar dt_ = 0.0;
    std::vector<Scalar> eigenZ_;
    std::vector<Scalar> eigenValues_;
    std::array<double, 3> initialInventory_{};
    std::array<double, 3> cumulativeOut_{};
    std::ofstream budget_;
};

} // namespace Dumux

#endif
