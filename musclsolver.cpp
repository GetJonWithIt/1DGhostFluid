#include "musclsolver.h"

// This class encapsulates the second-order, slope-limited, TVD MUSCL-Hancock solver for the (stiffened gas) Euler equations, as detailed in Toro.
MUSCLSolver::MUSCLSolver()
{
}

// Computes the MUSCL-Hancock flux, using the interface between four neighbouring cells, with states given by leftLeftStateVector, leftStateVector, rightStateVector. and rightRightStateVector,
// respetively.
vector<double> MUSCLSolver::computeMUSCLFlux(StateVector leftLeftStateVector, StateVector leftStateVector, StateVector rightStateVector, StateVector rightRightStateVector,
                                             double cellSpacing, double timeStep, double bias, int limiter)
{
    StateVector rightEvolvedState = Solvers::computeEvolvedState(leftLeftStateVector, leftStateVector, rightStateVector, cellSpacing, timeStep, bias, limiter, 1);
    StateVector leftEvolvedState = Solvers::computeEvolvedState(leftStateVector, rightStateVector, rightRightStateVector, cellSpacing, timeStep, bias, limiter, 0);

    return HLLCSolver::computeHLLCFlux(rightEvolvedState, leftEvolvedState);
}

// Evolves the computational domain by one timestep, using the MUSCL-Hancock scheme.
void MUSCLSolver::computeMUSCLTimeStep(vector<StateVector> & newCells, vector<StateVector> & currentCells, double cellSpacing, double timeStep, double bias, int limiter)
{
    double cellCount = newCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        vector<double> leftFluxVector = computeMUSCLFlux(currentCells[i], currentCells[i + 1], currentCells[i + 2], currentCells[i + 3], cellSpacing, timeStep, bias, limiter);
        vector<double> rightFluxVector = computeMUSCLFlux(currentCells[i + 1], currentCells[i + 2], currentCells[i + 3], currentCells[i + 4], cellSpacing, timeStep, bias, limiter);
        vector<double> conservedVariableVector = newCells[i].computeConservedVariableVector();

        newCells[i].setConservedVariableVector(LinearAlgebra::addVectors(conservedVariableVector, LinearAlgebra::multiplyVector(
                                                                             (timeStep / cellSpacing), LinearAlgebra::subtractVectors(leftFluxVector, rightFluxVector))));
    }
}

// Evolves the computational domain until finalTime, using the MUSCL-Hancock scheme.
vector<StateVector> MUSCLSolver::solve(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int limiter)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<StateVector> currentCells = cells;

    while (currentTime < finalTime)
    {
        vector<StateVector> currentCellsWithBoundary = Solvers::insertBoundaryCells(currentCells, 2);
        double timeStep = Solvers::computeStableTimeStep(currentCellsWithBoundary, cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

        computeMUSCLTimeStep(currentCells, currentCellsWithBoundary, cellSpacing, timeStep, bias, limiter);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return currentCells;
}

// Evolves the computational domain for a given multimaterial system until finalTime, using the MUSCL-Hancock scheme and the original ghost fluid method.
MultimaterialSystem MUSCLSolver::solveWithBasicGhostFluid(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int limiter)
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

        vector<StateVector> currentLeftMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentLeftMaterialCells, 2);
        vector<StateVector> currentRightMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentRightMaterialCells, 2);

        double timeStep = Solvers::computeGhostFluidStableTimeStep(MultimaterialSystem(currentLeftMaterialCellsWithBoundary, currentRightMaterialCellsWithBoundary, levelSetFunction),
                                                                   cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

        levelSetFunction = GhostFluidMethod::updateLevelSetFunction(levelSetFunction, currentLeftMaterialCells, currentRightMaterialCells, cellSpacing, timeStep);

        computeMUSCLTimeStep(currentLeftMaterialCells, currentLeftMaterialCellsWithBoundary, cellSpacing, timeStep, bias, limiter);
        computeMUSCLTimeStep(currentRightMaterialCells, currentRightMaterialCellsWithBoundary, cellSpacing, timeStep, bias, limiter);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction);
}

// Evolves the computational domain for a given multimaterial system until finalTime, using the MUSCL-Hancock scheme and the original ghost fluid method, with level set reinitialisation
// (with the period of reinitialisation given by reinitialisationIteration).
MultimaterialSystem MUSCLSolver::solveWithBasicGhostFluidReinitialisation(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                          double bias, int limiter, int reinitialisationIteration)
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

        vector<StateVector> currentLeftMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentLeftMaterialCells, 2);
        vector<StateVector> currentRightMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentRightMaterialCells, 2);

        double timeStep = Solvers::computeGhostFluidStableTimeStep(MultimaterialSystem(currentLeftMaterialCellsWithBoundary, currentRightMaterialCellsWithBoundary, levelSetFunction),
                                                                   cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

        levelSetFunction = GhostFluidMethod::updateLevelSetFunction(levelSetFunction, currentLeftMaterialCells, currentRightMaterialCells, cellSpacing, timeStep);

        if ((currentIteration % reinitialisationIteration) == (reinitialisationIteration - 1))
        {
            if (GhostFluidMethod::computeInterfaceCount(levelSetFunction) == 2)
            {
                levelSetFunction = GhostFluidMethod::reinitialiseTwoInterfaceLevelSetFunction(levelSetFunction);
            }
            else
            {
                levelSetFunction = GhostFluidMethod::reinitialiseLevelSetFunction(levelSetFunction);
            }
        }

        computeMUSCLTimeStep(currentLeftMaterialCells, currentLeftMaterialCellsWithBoundary, cellSpacing, timeStep, bias, limiter);
        computeMUSCLTimeStep(currentRightMaterialCells, currentRightMaterialCellsWithBoundary, cellSpacing, timeStep, bias, limiter);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction);
}

// Evolves the computational domain for a given multimaterial system until finalTime, using the MUSCL-Hancock scheme and the real ghost fluid method.
MultimaterialSystem MUSCLSolver::solveWithRiemannGhostFluid(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int limiter)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<StateVector> currentLeftMaterialCells = leftMaterialCells;
    vector<StateVector> currentRightMaterialCells = rightMaterialCells;

    double timeStep = Solvers::computeGhostFluidStableTimeStep(MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction),
                                                               cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

    while (currentTime < finalTime)
    {
        MultimaterialSystem newMultimaterialSystem = GhostFluidMethod::applyRiemannGhostFluidBoundaryConditions(
                    MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction));

        currentLeftMaterialCells = newMultimaterialSystem.getLeftMaterialCells();
        currentRightMaterialCells = newMultimaterialSystem.getRightMaterialCells();

        vector<StateVector> currentLeftMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentLeftMaterialCells, 2);
        vector<StateVector> currentRightMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentRightMaterialCells, 2);

        timeStep = Solvers::computeGhostFluidStableTimeStep(MultimaterialSystem(currentLeftMaterialCellsWithBoundary, currentRightMaterialCellsWithBoundary, levelSetFunction),
                                                                   cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

        levelSetFunction = GhostFluidMethod::updateLevelSetFunction(levelSetFunction, currentLeftMaterialCells, currentRightMaterialCells, cellSpacing, timeStep);

        computeMUSCLTimeStep(currentLeftMaterialCells, currentLeftMaterialCellsWithBoundary, cellSpacing, timeStep, bias, limiter);
        computeMUSCLTimeStep(currentRightMaterialCells, currentRightMaterialCellsWithBoundary, cellSpacing, timeStep, bias, limiter);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction);
}

// Evolves the computational domain for a given multimaterial system until finalTime, using the MUSCL-Hancock scheme and the real ghost fluid method, with level set reinitialisation
// (with the period of reinitialisation given by reinitialisationIteration).
MultimaterialSystem MUSCLSolver::solveWithRiemannGhostFluidReinitialisation(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime,
                                                                            double bias, int limiter, int reinitialisationIteration)
{
    double currentTime = 0.0;
    int currentIteration = 0;
    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<StateVector> currentLeftMaterialCells = leftMaterialCells;
    vector<StateVector> currentRightMaterialCells = rightMaterialCells;

    double timeStep = Solvers::computeGhostFluidStableTimeStep(MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction),
                                                               cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

    while (currentTime < finalTime)
    {
        MultimaterialSystem newMultimaterialSystem = GhostFluidMethod::applyRiemannGhostFluidBoundaryConditions(
                    MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction));

        currentLeftMaterialCells = newMultimaterialSystem.getLeftMaterialCells();
        currentRightMaterialCells = newMultimaterialSystem.getRightMaterialCells();

        vector<StateVector> currentLeftMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentLeftMaterialCells, 2);
        vector<StateVector> currentRightMaterialCellsWithBoundary = Solvers::insertBoundaryCells(currentRightMaterialCells, 2);

        timeStep = Solvers::computeGhostFluidStableTimeStep(MultimaterialSystem(currentLeftMaterialCellsWithBoundary, currentRightMaterialCellsWithBoundary, levelSetFunction),
                                                            cellSpacing, CFLCoefficient, currentTime, finalTime, currentIteration);

        levelSetFunction = GhostFluidMethod::updateLevelSetFunction(levelSetFunction, currentLeftMaterialCells, currentRightMaterialCells, cellSpacing, timeStep);

        if ((currentIteration % reinitialisationIteration) == (reinitialisationIteration - 1))
        {
            if (GhostFluidMethod::computeInterfaceCount(levelSetFunction) == 2)
            {
                levelSetFunction = GhostFluidMethod::reinitialiseTwoInterfaceLevelSetFunction(levelSetFunction);
            }
            else
            {
                levelSetFunction = GhostFluidMethod::reinitialiseLevelSetFunction(levelSetFunction);
            }
        }

        computeMUSCLTimeStep(currentLeftMaterialCells, currentLeftMaterialCellsWithBoundary, cellSpacing, timeStep, bias, limiter);
        computeMUSCLTimeStep(currentRightMaterialCells, currentRightMaterialCellsWithBoundary, cellSpacing, timeStep, bias, limiter);

        currentTime += timeStep;
        currentIteration += 1;
    }

    return MultimaterialSystem(currentLeftMaterialCells, currentRightMaterialCells, levelSetFunction);
}
