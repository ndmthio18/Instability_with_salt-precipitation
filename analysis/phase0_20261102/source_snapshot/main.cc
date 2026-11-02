// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright 漏 DuMux-Course contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>
#include <type_traits>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <map>
#include <cmath>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

#include <dumux/io/rasterimagereader.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/discretization/box.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>


#include "model.hh"

/*!
 * \ingroup RichardsSalt
 * \brief Test problem for image denoising using a nonlinear diffusion model
 */
namespace Dumux
{
    template <class TypeTag>
    class RichardsSaltTestProblem : public FVProblem<TypeTag>
    {
        using ParentType = FVProblem<TypeTag>;
        using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
        using GridView = typename GridGeometry::GridView;
        using Element = typename GridView::template Codim<0>::Entity;
        using GlobalPosition = typename GridGeometry::LocalView::Element::Geometry::GlobalCoordinate;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
        using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
        using SubControlVolume = typename FVElementGeometry::SubControlVolume;
        using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
        using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
        using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
        using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
        using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
        using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
        using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
        using Vertex = typename GridView::template Codim<GridView::dimensionworld>::Entity;

        // copy some indices for convenience
        enum
        {
            pressureIdx = Indices::pressureIdx,
            moleFractionIdx = Indices::moleFractionIdx,

            conti0EqIdx = Indices::conti0EqIdx,
            conti1EqIdx = Indices::conti1EqIdx,
        };

        enum class PerturbationType
        {
            Periodic,
            Constant,
            Random,
            Eigenprofile
        };

        // Function to convert string to PerturbationType
        static PerturbationType stringToPerturbationType(const std::string &str)
        {
            if (str == "Periodic")
            {
                return PerturbationType::Periodic;
            }
            else if (str == "Constant")
            {
                return PerturbationType::Constant;
            }
            else if (str == "Random")
            {
                return PerturbationType::Random;
            }
            else if (str == "Eigenprofile")
            {
                return PerturbationType::Eigenprofile;
            }
            else
            {
                throw std::invalid_argument("Invalid perturbation type: " + str);
            }
        }

    public:
        RichardsSaltTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
            : ParentType(gridGeometry)
        {
            permeability_ = getParam<Scalar>("Problem.Permeability");
            porosity_ = getParam<Scalar>("Problem.Porosity");

            saltMoleFraction_ = getParam<Scalar>("Problem.SaltMoleFraction");
            stdSaltMoleFraction_ = getParam<Scalar>("Problem.StandardDeviationSaltMoleFraction");

            waterLevel_ = getParam<Scalar>("Problem.WaterLevel");

            amplitude_ = getParam<Scalar>("Problem.Amplitude");
            wavenumber_ = getParam<Scalar>("Problem.WaveNumber");

            perturbationType_ = stringToPerturbationType(getParam<std::string>("Problem.PerturbationType"));

            unsigned seed = getParam<unsigned>("Problem.RandomSeed", 1);
            std::default_random_engine generator(seed);
            std::normal_distribution<double> distribution(saltMoleFraction_, stdSaltMoleFraction_);

            auto numVertex = this->gridGeometry().vertexMapper().size();
            randomSaltMoleFraction_.resize(numVertex, 0.0);

            for (const auto &vertex : vertices(this->gridGeometry().gridView()))
            {
                int vIdx = this->gridGeometry().vertexMapper().index(vertex);
                randomSaltMoleFraction_[vIdx] = distribution(generator);
            }

            // Load the tabulated eigenprofile A(y) if the Eigenprofile perturbation is used.
            if (perturbationType_ == PerturbationType::Eigenprofile)
            {
                const auto path = getParam<std::string>("Problem.EigenprofilePath");
                std::ifstream in(path);
                if (!in)
                    throw std::runtime_error("Cannot open Problem.EigenprofilePath: " + path);
                std::string line;
                while (std::getline(in, line))
                {
                    if (line.empty() || line[0] == '#')
                        continue;
                    std::istringstream iss(line);
                    Scalar yy, aa;
                    if (iss >> yy >> aa)
                    {
                        eigenY_.push_back(yy);
                        eigenA_.push_back(aa);
                    }
                }
                if (eigenY_.size() < 2)
                    throw std::runtime_error("Eigenprofile table needs >=2 rows: " + path);
            }
        }

        // Linear interpolation of the tabulated eigenprofile A at depth coordinate y
        // (y assumed ascending in the table; clamps outside the tabulated range).
        Scalar eigenAmplitude_(Scalar y) const
        {
            if (y <= eigenY_.front()) return eigenA_.front();
            if (y >= eigenY_.back())  return eigenA_.back();
            auto it = std::lower_bound(eigenY_.begin(), eigenY_.end(), y);
            const std::size_t i = static_cast<std::size_t>(it - eigenY_.begin());
            const Scalar y1 = eigenY_[i - 1], y2 = eigenY_[i];
            const Scalar a1 = eigenA_[i - 1], a2 = eigenA_[i];
            return a1 + (a2 - a1) * (y - y1) / (y2 - y1);
        }

        BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
        {
            BoundaryTypes bcTypes;
            if (onLowerBoundary_(globalPos))
                bcTypes.setAllDirichlet();
            else
                bcTypes.setAllNeumann();
            return bcTypes;
        }
        PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
        {
            PrimaryVariables values(0.0);

            const Scalar densityW = rhoW(saltMoleFraction_);
            Scalar hydrostaticPressure = -(waterLevel_ + globalPos[1])*densityW*gravitationalAcceleration();
            values[pressureIdx] = nonWettingReferencePressure() + hydrostaticPressure;
            
            // We assume here that at the bottom boundary is constant
            values[moleFractionIdx] = saltMoleFraction_ ; 

            return values;
        }

        NumEqVector neumann(const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const ElementVolumeVariables &elemVolVars,
                            const ElementFluxVariablesCache &elemFluxVarsCache,
                            const SubControlVolumeFace &scvf) const
        {
            NumEqVector values(0.0);
           
           const auto globalPos = scvf.ipGlobal();
            if (onUpperBoundary_(globalPos))
            {
                static const Scalar evaporationRate =  getParam<Scalar>("Problem.EvaporationRate");  // m/s
                values[conti0EqIdx] = evaporationRate;
                values[conti1EqIdx] = 0;
                
                
            }
            
            return values;
        }
      

        PrimaryVariables initial(const Vertex &vertex) const
        {
            PrimaryVariables values(0.0);
            const auto globalPos = vertex.geometry().center();
            const Scalar densityW = rhoW(saltMoleFraction_);
            Scalar hydrostaticPressure = -(waterLevel_ + globalPos[1])*densityW*gravitationalAcceleration();
            values[pressureIdx] = nonWettingReferencePressure() + hydrostaticPressure;
            // Option A: Enable for analysis of the perturbation height 
            static const Scalar perturbationHeight = getParam<Scalar>("Problem.PerturbationHeight"); 

            switch (perturbationType_)
            {
            case PerturbationType::Periodic:
            {
                // Option A: Enable for analysis of the perturbation height 
                 if (globalPos[1] >= perturbationHeight)
                // Option B: Enable for analysis of the perturbing only the upper boundary 
                //if (globalPos[1] == this->gridGeometry().bBoxMax()[1])
                {
                    Scalar wavelength = (this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0]) / wavenumber_;
                    values[moleFractionIdx] = saltMoleFraction_ + amplitude_ * std::cos(2 * M_PI / wavelength * globalPos[0]);
                }
                else
                {
                   // No perturbation for cells below the perturbation height, set to default mole fraction
                    values[moleFractionIdx] = saltMoleFraction_;
                }
                // Option C: Enable for eigenprofile from linear stability analysis  
                // Note: Each case has its own eigenprofile from which we derive the corresponding A(z)-profile. For this, we need to first estimate, the eigenprofile curve. 
                // Scalar amplitude = 0.0; 
                // Adjust the region where we apply the A(z)-profile. Below we assume that A(z) = 0
                // if (globalPos[1]>= 0.125)
                // {
                //eigenprofile WL025 WN8     0.125
                // const auto coeff_a = 7.638942176065109e-13;
                // const auto coeff_b = 5.28770754848503e-05; 
                // const auto coeff_c = 0.008530870953167433; 
                // const auto coeff_d = 0.5796464075334609;
                // const auto coeff_e = 22.160923408812245;
                // const auto coeff_f = 534.1755833099919;
                // const auto coeff_g = 8578.394972202852;
                // const auto coeff_h = 94133.0822872317;                
                // const auto coeff_i = 708823.6280653266;
                // const auto coeff_j = 3602440.6779811466;
                // const auto coeff_k = 11805813.111630982;
                // const auto coeff_l = 22514896.328417305;
                // const auto coeff_m = 18974942.7858038;
                
                // eigenprofile WL030 WN8 0.075
                // const auto coeff_a = -1.7911296835278794e-12;
                // const auto coeff_b = 0.00014742699233301506;
                // const auto coeff_c = 0.022730824811586017;
                // const auto coeff_d = 1.4620596234361638;
                // const auto coeff_e = 52.54085963848574;
                // const auto coeff_f = 1184.3902150141814;
                // const auto coeff_g = 17725.1279233307;
                // const auto coeff_h = 180839.309667342;
                // const auto coeff_i = 1264306.5295694193;
                // const auto coeff_j = 5961877.929792657;
                // const auto coeff_k = 18126381.000921432;
                // const auto coeff_l = 32082935.12166746;
                // const auto coeff_m = 25113300.115454108;
                
                // eigenprofile WL035 WN18  0.125
                // const auto coeff_a = 6.616474872066744e-14;
                // const auto coeff_b = 0.0005899180182105107;
                // const auto coeff_c = 0.10395116383425879;
                // const auto coeff_d = 7.9179033397135505;
                // const auto coeff_e = 349.0268064212656;
                // const auto coeff_f = 10010.961339586056;
                // const auto coeff_g = 198423.36461733107;
                // const auto coeff_h = 2808275.962015743;
                // const auto coeff_i = 28834974.07325553;
                // const auto coeff_j = 215409449.6336354;
                // const auto coeff_k = 1159072330.0471587;
                // const auto coeff_l = 4375217567.562299;
                // const auto coeff_m = 10990672557.639498;
                // const auto coeff_n = 16492835368.956154;
                // const auto coeff_o = 11180991128.11837;

                // eigenprofile WL40 WN16  0.125
                // const auto coeff_a = -5.160018139221863e-14;
                // const auto coeff_b = 0.0026652014086449606;
                // const auto coeff_c = 0.46022343367175916;
                // const auto coeff_d = 34.22038402471369;
                // const auto coeff_e = 1468.1732115907091;
                // const auto coeff_f = 40894.66036972125;
                // const auto coeff_g = 785906.3906962734;
                // const auto coeff_h = 10774067.00292978;
                // const auto coeff_i = 107110433.15486877;
                // const auto coeff_j = 774770057.6719315;
                // const auto coeff_k = 4038591415.3715277;
                // const auto coeff_l = 14781489071.715094;
                // const auto coeff_m = 36048297196.70141;
                // const auto coeff_n = 52599322714.29569;
                // const auto coeff_o = 34737467870.16858;


                // amplitude = coeff_a + coeff_b*globalPos[1]- coeff_c*pow(globalPos[1],2)+ coeff_d*pow(globalPos[1],3)-coeff_e*pow(globalPos[1],4)+coeff_f*pow(globalPos[1],5)-coeff_g*pow(globalPos[1],6)+coeff_h*pow(globalPos[1],7)-coeff_i*pow(globalPos[1],8)+coeff_j*pow(globalPos[1],9)-coeff_k*pow(globalPos[1],10)+coeff_l*pow(globalPos[1],11)-coeff_m*pow(globalPos[1],12)+coeff_n*pow(globalPos[1],13)-coeff_o*pow(globalPos[1],14); 
                // }


                // Scalar wavelength = (this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0]) / wavenumber_;
                // values[moleFractionIdx] = saltMoleFraction_ + amplitude * std::cos(2 * M_PI / wavelength * globalPos[0]);

                break;
            }
            case PerturbationType::Constant:
            {
                values[moleFractionIdx] = saltMoleFraction_;
                break;
            }
            case PerturbationType::Random:
            {
                // Respect PerturbationHeight (as the Periodic branch does): set
                // PerturbationHeight just below the top to perturb only the top row of
                // vertices, matching the paper's "random perturbation in the top row".
                if (globalPos[1] >= perturbationHeight)
                {
                    int vertexIdx = this->gridGeometry().vertexMapper().index(vertex);
                    values[moleFractionIdx] = randomSaltMoleFraction_[vertexIdx];
                }
                else
                {
                    values[moleFractionIdx] = saltMoleFraction_;
                }
                break;
            }
            case PerturbationType::Eigenprofile:
            {
                // Most-unstable LSA eigenmode A(y) (tabulated) modulated laterally by the
                // fastest-growing wavelength: x(x0,y) = x0 + amplitude * A(y) * cos(2 pi x / lambda).
                Scalar wavelength = (this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0]) / wavenumber_;
                Scalar A = eigenAmplitude_(globalPos[1]);
                values[moleFractionIdx] = saltMoleFraction_ + amplitude_ * A * std::cos(2 * M_PI / wavelength * globalPos[0]);
                break;
            }
            }

            return values;
        }


        Scalar permeability() const
        {
            return permeability_;
        }

        // Accessors for the onset diagnostics assembled in main(). The seeded mode is
        // A(y) cos(2 pi x / lambda); reporting its amplitude directly is a far more
        // sensitive measure of that mode than the spread of the top row, which mixes it
        // with every other lateral structure present.
        Scalar eigenAmplitudeAt(Scalar y) const
        {
            return eigenAmplitude_(y);
        }

        Scalar waveNumber() const
        {
            return wavenumber_;
        }

        bool hasEigenprofile() const
        {
            return perturbationType_ == PerturbationType::Eigenprofile;
        }

        Scalar porosity() const
        {
            return porosity_;
        }

        Scalar nonWettingReferencePressure() const
        {
            return 1e5;
        }

        Scalar effectiveSw(Scalar sw) const
        {
            static const Scalar swr = getParam<Scalar>("Problem.Swr");
            return (sw - swr) / (1 - swr);
        }

        Scalar sw(Scalar pw) const
        {
            static const Scalar n = getParam<Scalar>("Problem.PcN");
            static const Scalar m = 1.0 - 1.0 / n;
            static const Scalar alpha = getParam<Scalar>("Problem.PcAlpha");
            static const Scalar swr = getParam<Scalar>("Problem.Swr");
            const auto pc = nonWettingReferencePressure() - pw;
            return swr + (1 - swr) / (std::pow(1 + std::pow(alpha * pc, n), m));
            if (pw>nonWettingReferencePressure())
                return 1.0; 
            else{
                return swr + (1-swr)/(std::pow(1+std::pow(alpha*pc,n),m));
            }
        }

        // Gravitational acceleration [m/s^2], used CONSISTENTLY by the hydrostatic
        // initial/boundary profiles and by the momentum balance. Default 9.81 reproduces
        // pristine upstream darus-4610, which is what produced the published numerical
        // results of Bringedal et al. (2025) (Tables 5.8, 5.9). Set 9.8 to match the
        // value stated in Bringedal et al. (2022) Table 2 / (2025) Table 5.1 -- note the
        // authors' own LSA driver implies 9.80015, so their two legs differ by 0.1 %.
        // Before this parameter existed our patch set the momentum term to 9.8 while
        // leaving the hydrostatic profiles at 9.81, which left the initial state out of
        // equilibrium and drove a spurious Darcy flux ~18x the evaporation rate.
        Scalar gravitationalAcceleration() const
        {
            static const Scalar g = getParam<Scalar>("Problem.GravitationalAcceleration", 9.81);
            return g;
        }

        Scalar krw(Scalar sw) const
        {
            static const Scalar n = getParam<Scalar>("Problem.PcN");
            static const Scalar m = 1.0 - 1.0 / n;
            static const Scalar sweL = getParam<Scalar>("Problem.SweL");

            // Problem.KrwMualemExponent selects how the Mualem parameter L enters.
            //  false (default): sweL as a constant FACTOR -> krw(Se=1) = sweL = 0.73.
            //                   This is the upstream darus-4610 behaviour and is what
            //                   the published numerical results of Bringedal et al.
            //                   (2025) were produced with. Kept as the default so no
            //                   existing result changes.
            //  true:            Se^L as an EXPONENT -> krw(Se=1) = 1. This is what the
            //                   paper's eq. (2.5)/(3.13) prescribe, and what the
            //                   authors' own LSA code (Fktilde.m, DaRUS darus-4711)
            //                   implements. See code/paper_comparison/UNSAT_NOTES.md.
            static const bool mualemExponent = getParam<bool>("Problem.KrwMualemExponent", false);

            const auto eSw = effectiveSw(sw);

            Scalar val = 1.0 - std::pow(1.0 - std::pow(eSw, 1.0 / m), m);
            return mualemExponent ? std::pow(eSw, sweL) * val * val
                                  : sweL * val * val;
        }

        Scalar viscosity() const
        {
            // Runtime-configurable so a dimensionless design point can be realized
            // without silently changing R_E.  The default preserves the published run.
            static const Scalar mu = getParam<Scalar>("Problem.ReferenceViscosity", 1.1e-3);
            return mu;
        }

        Scalar effectiveDiffusionCoefficient(Scalar sw) const
        {
            // Molecular salt diffusivity (before Millington-Quirk porosity/saturation
            // scaling).  Making it an input is required by the proposal_fixed_e
            // realization, where beta determines D for each design point.
            static const Scalar molecularD =
                getParam<Scalar>("Problem.SaltDiffusivity", 1.5e-9);
            return std::pow(porosity(), 1.5) * std::pow(sw, 3.5) * molecularD;
        }

        Scalar molarMassH2O() const
        {
            static const Scalar value = getParam<Scalar>("Problem.MolarMassH2O", 18e-3);
            return value; // [kg/mol]
        }

        Scalar molarMassNaCl() const
        {
            static const Scalar value = getParam<Scalar>("Problem.MolarMassNaCl", 58.4428e-3);
            return value; // [kg/mol]
        }

        Scalar rhoW(Scalar moleFraction) const
        {
            static const Scalar rhoH2O =
                getParam<Scalar>("Problem.PureWaterDensity", 1000.0); // kg/m3
            static const Scalar linear =
                getParam<Scalar>("Problem.DensityLinearCoefficient", 0.6708);
            static const Scalar quadratic =
                getParam<Scalar>("Problem.DensityQuadraticCoefficient", 0.3739);
            Scalar massFraction = moleFraction * molarMassNaCl() / ((1 - moleFraction) * molarMassH2O() + moleFraction * molarMassNaCl());

            return rhoH2O + 1000.0 * massFraction * (linear + quadratic * massFraction);
        }

        Scalar molarRhoW(Scalar moleFraction) const
        {
            Scalar molarMassBrine = (1 - moleFraction) * molarMassH2O() + moleFraction * molarMassNaCl();
            return rhoW(moleFraction) / molarMassBrine;
        }

        bool onUpperBoundary_(const GlobalPosition &globalPos) const
        {
            return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
        }

        bool onLowerBoundary_(const GlobalPosition &globalPos) const
        {
            return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
        }

    private:
        Scalar permeability_;
        Scalar porosity_;

        Scalar saltMoleFraction_;
        Scalar waterLevel_;

        PerturbationType perturbationType_;

        Scalar wavenumber_;
        Scalar amplitude_;

        std::vector<Scalar> randomSaltMoleFraction_;
        Scalar stdSaltMoleFraction_;

        // Eigenprofile reference perturbation: tabulated A(y) read from a file.
        std::vector<Scalar> eigenY_;
        std::vector<Scalar> eigenA_;

        static constexpr Scalar eps_ = 1.5e-7;
    };
} // end namespace Dumux

// Property tag and specializations
namespace Dumux::Properties::TTag
{

    struct RichardsSaltTest
    {
        using InheritsFrom = std::tuple<RichardsSaltModel, BoxModel>;

        using Scalar = double;
        using Grid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2>>;

        // Set the above implemented problem
        template <class TypeTag>
        using Problem = RichardsSaltTestProblem<TypeTag>;

        using EnableGridVolumeVariablesCache = std::true_type;
        using EnableGridFluxVariablesCache = std::true_type;
        using EnableGridGeometryCache = std::true_type;
    };

} // end namespace Dumux::Properties::TTag

// the main program
int main(int argc, char **argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::RichardsSaltTest;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto &leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

    using GlobalPosition = typename GridGeometry::LocalView::Element::Geometry::GlobalCoordinate;
    std::vector<GlobalPosition> cellVelocity(leafGridView.size(0)); // or appropriate size

    static const bool enableGravity = getParam<bool>("Problem.EnableGravity");
    GlobalPosition gravity(0.0);
    if (enableGravity)
        gravity[GridGeometry::GridView::dimension-1] = -getParam<Scalar>("Problem.GravitationalAcceleration", 9.81);

    // initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    
    vtkWriter.addVolumeVariable([](const auto &vv)
                                { return vv.pw(); }, "pw");
    vtkWriter.addVolumeVariable([](const auto &vv)
                                { return vv.sw(); }, "sw");
    vtkWriter.addVolumeVariable([](const auto &vv)
                                { return vv.moleFraction(); }, "moleFraction");
    vtkWriter.addVolumeVariable([](const auto &vv)
                                { return vv.rhoW(); }, "rhoW"); 
    vtkWriter.addField(cellVelocity, "velocity");                                                                       

    vtkWriter.write(0.0);

    // ---------------------------------------------------------------------------------
    // Onset diagnostic, written every time step independently of the VTU cadence.
    //
    // The onset time is the minimum over time of the spread of the salt mole fraction
    // along the top boundary (Bringedal et al. 2025, eqs. 5.1-5.2; same diagnostic as
    // the 2022 paper and as code/paper_comparison/extract_onset.py). Writing it here
    // decouples the quantity of interest from the VTU output, which at the paper's
    // protocol (TEnd = 604800 s, dt = 50 s) would otherwise cost ~56 GB per run.
    //
    // `sigma_top` below is the population standard deviation, i.e. what numpy.std gives.
    // The paper's sigma^top is the unnormalised line integral sqrt(int (x - xbar)^2 dx);
    // for the uniform lateral spacing used here the two differ by the constant
    // sqrt(domain width) and therefore share the same argmin.
    // ---------------------------------------------------------------------------------
    static const int outputInterval = getParam<int>("TimeLoop.OutputInterval", 1);

    using ModelIndices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    std::vector<std::size_t> topDofs;
    {
        const auto yTop = gridGeometry->bBoxMax()[GridGeometry::GridView::dimension-1];
        for (const auto &vertex : vertices(leafGridView))
        {
            const auto pos = vertex.geometry().center();
            if (std::abs(pos[GridGeometry::GridView::dimension-1] - yTop) < 1e-12)
                topDofs.push_back(gridGeometry->vertexMapper().index(vertex));
        }
    }

    // Two further amplitudes are written beside sigma_top, both projections onto the mode
    // that was actually seeded rather than onto whatever the top row happens to contain.
    //
    //   fourier_top : |sum_j (X_j - Xbar) exp(-i 2 pi n_c x_j / W)| * 2/N over the top row.
    //                 Isolates the seeded lateral wavenumber from the rest of the row.
    //   mode_proj   : <X(x,y) - Xbar(y), A(y) cos(2 pi x / lambda)> / <basis, basis> over
    //                 every vertex. The matched filter for the seeded mode: it uses the
    //                 whole depth instead of one row, which for a mode confined to a few
    //                 millimetres of a 200 mm column is where nearly all of its signal is.
    //
    // Both are diagnostics only; they are written alongside the existing columns and change
    // nothing about the solution or the frozen sigma_top definition.
    const Scalar xMinDomain = gridGeometry->bBoxMin()[0];
    const Scalar xMaxDomain = gridGeometry->bBoxMax()[0];
    const Scalar domainWidth = xMaxDomain - xMinDomain;
    const bool haveMode = problem->hasEigenprofile();
    const Scalar lambdaSeed = haveMode ? domainWidth / problem->waveNumber() : domainWidth;

    // Per-vertex basis values and the depth level each vertex belongs to, so the lateral
    // mean can be removed level by level before projecting.
    std::vector<Scalar> basis(gridGeometry->numDofs(), 0.0);
    std::vector<std::size_t> level(gridGeometry->numDofs(), 0);
    std::vector<Scalar> dofX(gridGeometry->numDofs(), 0.0);
    std::vector<Scalar> levelY;
    Scalar basisNormSq = 0.0;
    {
        constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
        std::map<long long, std::size_t> levelOf;
        for (const auto &vertex : vertices(leafGridView))
        {
            const auto pos = vertex.geometry().center();
            const auto dof = gridGeometry->vertexMapper().index(vertex);
            const Scalar yy = pos[dimWorld-1];
            const long long key = static_cast<long long>(std::llround(yy*1e9));
            auto it = levelOf.find(key);
            if (it == levelOf.end())
            {
                it = levelOf.emplace(key, levelY.size()).first;
                levelY.push_back(yy);
            }
            level[dof] = it->second;
            dofX[dof] = pos[0];
            basis[dof] = (haveMode ? problem->eigenAmplitudeAt(yy) : 1.0)
                         * std::cos(2*M_PI/lambdaSeed * pos[0]);
            basisNormSq += basis[dof]*basis[dof];
        }
    }

    std::ofstream sigmaFile(problem->name() + "_sigmatop.csv");
    sigmaFile << "# time[s] x_top_mean[mol/mol] sigma_top[mol/mol] n_top_dofs"
                 " fourier_top[mol/mol] mode_proj[mol/mol]\n";
    sigmaFile << std::scientific << std::setprecision(12);

    auto writeSigmaTop = [&](Scalar time)
    {
        if (topDofs.empty()) return;
        Scalar mean = 0.0;
        for (const auto dof : topDofs)
            mean += x[dof][ModelIndices::moleFractionIdx];
        mean /= static_cast<Scalar>(topDofs.size());
        Scalar var = 0.0;
        for (const auto dof : topDofs)
        {
            const Scalar d = x[dof][ModelIndices::moleFractionIdx] - mean;
            var += d*d;
        }
        var /= static_cast<Scalar>(topDofs.size());

        Scalar re = 0.0, im = 0.0;
        for (std::size_t k = 0; k < topDofs.size(); ++k)
        {
            const Scalar d = x[topDofs[k]][ModelIndices::moleFractionIdx] - mean;
            const Scalar phase = 2*M_PI/lambdaSeed * dofX[topDofs[k]];
            re += d*std::cos(phase);
            im -= d*std::sin(phase);
        }
        const Scalar fourierTop = 2.0*std::sqrt(re*re + im*im)
                                  / static_cast<Scalar>(topDofs.size());

        std::vector<Scalar> levelSum(levelY.size(), 0.0);
        std::vector<std::size_t> levelCount(levelY.size(), 0);
        for (std::size_t dof = 0; dof < basis.size(); ++dof)
        {
            levelSum[level[dof]] += x[dof][ModelIndices::moleFractionIdx];
            ++levelCount[level[dof]];
        }
        for (std::size_t l = 0; l < levelY.size(); ++l)
            if (levelCount[l] > 0) levelSum[l] /= static_cast<Scalar>(levelCount[l]);
        Scalar proj = 0.0;
        for (std::size_t dof = 0; dof < basis.size(); ++dof)
            proj += (x[dof][ModelIndices::moleFractionIdx] - levelSum[level[dof]])*basis[dof];
        proj = (basisNormSq > 0.0) ? std::abs(proj)/basisNormSq : 0.0;

        sigmaFile << time << " " << mean << " " << std::sqrt(var) << " "
                  << topDofs.size() << " " << fourierTop << " " << proj << "\n";
        sigmaFile.flush();
    };
    writeSigmaTop(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    //using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start();
    do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        for (const auto& element : elements(leafGridView)) {
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(gridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, x);
        const auto eIdx = gridGeometry->elementMapper().index(element);
        const auto geometry = element.geometry();

        const auto& localBasis = fvGeometry.feLocalBasis();
        using ShapeValue = typename Dune::FieldVector<double, 1>;
        std::vector<ShapeValue> shapeValues;

        // Evaluate gradP - rho*g at integration point
        Dune::FieldVector<Scalar, GridGeometry::GridView::dimension> gradP(0.0);
        const auto ipLocal = geometry.local(geometry.center());
        localBasis.evaluateFunction(ipLocal, shapeValues);
        using ShapeJacobian = typename std::decay_t<decltype(localBasis)>::Traits::JacobianType;
        std::vector<ShapeJacobian> shapeJacobian;
        const auto jacInvT = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, shapeJacobian);

        Scalar cellDensity(0.0);
        Scalar cellViscosity(0.0);
        for (auto&& scv : scvs(fvGeometry)) {
            const auto& volVars = elemVolVars[scv];

            // The global shape function gradient
            GlobalPosition gradN(0.0);
            jacInvT.mv(shapeJacobian[scv.localDofIndex()][0], gradN);
            gradP.axpy(volVars.pw(), gradN);
            cellDensity += shapeValues[scv.indexInElement()][0] * volVars.rhoW();
            cellViscosity += shapeValues[scv.indexInElement()][0] * problem->viscosity();
        }

        // Add gravity term
        gradP.axpy(-cellDensity, gravity);
        // Obtain element permeability from spatial parameters
        cellVelocity[eIdx] = -problem->permeability() * gradP;
        cellVelocity[eIdx] /= cellViscosity;
    }


        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance the time loop to the next step
        timeLoop->advanceTimeStep();

        // the onset diagnostic, always; the VTU only every OutputInterval-th step
        writeSigmaTop(timeLoop->time());
        if (outputInterval <= 1 || timeLoop->timeStepIndex() % outputInterval == 0
            || timeLoop->finished())
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());


    return 0;
} // end main
