#include <config.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
try
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (mpiHelper.size() != 1)
        DUNE_THROW(Dune::InvalidStateException, "P1 diagnostics currently require one MPI rank");
    DumuxMessage::print(true);
    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& gridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const Scalar tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const Scalar maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const Scalar initialDt = getParam<Scalar>("TimeLoop.DtInitial");
    const int outputEvery = getParam<int>("TimeLoop.OutputEverySteps", 0);
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, initialDt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter);
    if (outputEvery > 0)
        vtkWriter.write(0.0);

    const auto initialAudit = problem->audit(x, *gridVariables);
    problem->initializeBudget(initialAudit);
    problem->writeCellState("cell_state_initial.csv", x, *gridVariables);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridView, gridGeometry->dofMapper());
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    std::size_t acceptedSteps = 0;
    timeLoop->start();
    do
    {
        const Scalar acceptedDt = timeLoop->timeStepSize();
        problem->setTimeState(timeLoop->time(), acceptedDt);
        nonLinearSolver.solve(x, *timeLoop);

        const Scalar acceptedTime = timeLoop->time() + acceptedDt;
        const auto audit = problem->audit(x, *gridVariables);
        problem->record(audit, acceptedTime, acceptedDt);

        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
        ++acceptedSteps;

        if (outputEvery > 0 && (acceptedSteps % outputEvery == 0 || timeLoop->finished()))
            vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(std::min(maxDt,
            nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize())));
    }
    while (!timeLoop->finished());

    timeLoop->finalize(gridView.comm());
    problem->writeCellState("cell_state_final.csv", x, *gridVariables);
    std::ofstream status("run_status.json");
    status << std::setprecision(17)
           << "{\"status\":\"completed\",\"accepted_steps\":" << acceptedSteps
           << ",\"time_s\":" << timeLoop->time() << "}\n";
    Parameters::print();
    DumuxMessage::print(false);
    return 0;
}
catch (const std::exception& error)
{
    std::ofstream status("run_status.json");
    status << "{\"status\":\"failed\"}\n";
    std::cerr << "P1_FAILURE " << error.what() << '\n';
    return 1;
}
