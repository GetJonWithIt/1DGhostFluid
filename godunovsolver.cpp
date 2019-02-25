#include "godunovsolver.h"

// This class encapsulates the general first-order Godunov-type solver for the (stiffened gas) Euler equations, and its multimaterial extension, as detailed in Toro.
GodunovSolver::GodunovSolver()
{
}

// Evolves the computational domain by one timestep, using a Godunov-type scheme.
void GodunovSolver::computeGodunovTimeStep(vector<StateVector> & newCells, vector<StateVector> & currentCells, double cellSpacing, double timeStep)
{
    int cellCount = newCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> conservedVariableVector = newCells[i].computeConservedVariableVector();
        vector<double> leftFluxVector = HLLCSolver::computeHLLCFlux(currentCells[i], currentCells[i + 1]);
        vector<double> rightFluxVector = HLLCSolver::computeHLLCFlux(currentCells[i + 1], currentCells[i + 2]);

        newCells[i].setConservedVariableVector(LinearAlgebra::addVectors(conservedVariableVector, LinearAlgebra::multiplyVector(
                                                                             (timeStep / cellSpacing), LinearAlgebra::subtractVectors(leftFluxVector, rightFluxVector))));
    }
}

// Evolves the computational domain until finalTime, using a Godunov-type scheme.
vector<StateVector> GodunovSolver::solve(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<StateVector> currentCells = cells;

    while (currentTime < finalTime)
    {
        vector<StateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 1);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

        computeGodunovTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return currentCells;
}

// Evolves the computational domain for a given multimaterial system until finalTime, using a Godunov-type scheme and the original ghost fluid method.
MultimaterialSystem GodunovSolver::solveWithBasicGhostFluid(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<StateVector> currentLeftMaterialCells = leftMaterialCells;
    vector<StateVector> currentRightMaterialCells = rightMaterialCells;

    while (currentTime < finalTime)
    {
        MultimaterialSystem newMultimaterialSystem = GhostFluidMethod::applyBasicGhostFluidBoundaryConditions(
                    MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction));
        currentLeftMaterialCells = newMultimaterialSystem.getLeftMaterialCells();
        currentRightMaterialCells = newMultimaterialSystem.getRightMaterialCells();

        vector<StateVector> currentLeftMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentLeftMaterialCells, 1);
        vector<StateVector> currentRightMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentRightMaterialCells, 1);

        double timeStep = min(Solvers::computeStableTimeStep(currentLeftMaterialCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration),
                              Solvers::computeStableTimeStep(currentRightMaterialCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration));

        levelSetFunction = GhostFluidMethod::updateLevelSetFunction(levelSetFunction, currentLeftMaterialCells, currentRightMaterialCells, cellSpacing, timeStep);

        computeGodunovTimeStep(currentLeftMaterialCells, currentLeftMaterialCellsWithBoundary, cellSpacing, timeStep);
        computeGodunovTimeStep(currentRightMaterialCells, currentRightMaterialCellsWithBoundary, cellSpacing, timeStep);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction);
}

// Evolves the computational domain for a given multimaterial system until finalTime, using a Godunov-type sceme and the real ghost fluid method.
MultimaterialSystem GodunovSolver::solveWithRiemannGhostFluid(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<StateVector> currentLeftMaterialCells = leftMaterialCells;
    vector<StateVector> currentRightMaterialCells = rightMaterialCells;

    double timeStep = min(Solvers::computeStableTimeStep(currentLeftMaterialCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration),
                          Solvers::computeStableTimeStep(currentRightMaterialCells, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration));

    while (currentTime < finalTime)
    {
        MultimaterialSystem newMultimaterialSystem = GhostFluidMethod::applyRiemannGhostFluidBoundaryConditions(
                    MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction));
        //MultimaterialSystem newMultimaterialSystem = GhostFluidMethod::applyRiemannGhostFluidBoundaryConditions(
        //            MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction), cellSpacing, timeStep);
        currentLeftMaterialCells = newMultimaterialSystem.getLeftMaterialCells();
        currentRightMaterialCells = newMultimaterialSystem.getRightMaterialCells();

        vector<StateVector> currentLeftMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentLeftMaterialCells, 1);
        vector<StateVector> currentRightMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentRightMaterialCells, 1);

        timeStep = min(Solvers::computeStableTimeStep(currentLeftMaterialCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration),
                              Solvers::computeStableTimeStep(currentRightMaterialCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration));

        levelSetFunction = GhostFluidMethod::updateLevelSetFunction(levelSetFunction, currentLeftMaterialCells, currentRightMaterialCells, cellSpacing, timeStep);

        computeGodunovTimeStep(currentLeftMaterialCells, currentLeftMaterialCellsWithBoundary, cellSpacing, timeStep);
        computeGodunovTimeStep(currentRightMaterialCells, currentRightMaterialCellsWithBoundary, cellSpacing, timeStep);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction);
}
