#include "solvers.h"

// This class encapsulates the methods that are common to all solvers for the (stiffened gas) Euler equations.
Solvers::Solvers()
{
}

// Returns the computational domain, but with a given number of boundary cells inserted on each side (transmissive boundary conditions).
vector<StateVector> Solvers::insertBoundaryCells(vector<StateVector> & cells, int boundaryCells)
{
    int cellCount = cells.size();
    vector<StateVector> cellsWithBoundary(cellCount + (2 * boundaryCells));

    if (boundaryCells == 1)
    {
        cellsWithBoundary[0] = cells[0];
        cellsWithBoundary[cellCount + 1] = cells[cellCount - 1];
    }

    if (boundaryCells == 2)
    {
        cellsWithBoundary[0] = cells[1];
        cellsWithBoundary[1] = cells[0];
        cellsWithBoundary[cellCount + 2] = cells[cellCount - 1];
        cellsWithBoundary[cellCount + 3] = cells[cellCount - 2];
    }

    for (int i = 0; i < cellCount; i++)
    {
        cellsWithBoundary[i + boundaryCells] = cells[i];
    }

    return cellsWithBoundary;
}

// Computes the estimated maximum wave speed across the entire computational domain, for the purposes of timestep calculation.
double Solvers::computeMaximumWaveSpeed(vector<StateVector> & cells)
{
    double maximumWaveSpeed = 0.0;
    int cellCount = cells.size();

    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = abs(cells[i].getXVelocity()) + cells[i].computeSoundSpeed();
        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

// Computes the maximum stable timestep, in accordance with the CFL condition.
double Solvers::computeStableTimeStep(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeMaximumWaveSpeed(cells));

    if (currentIteration <= 5)
    {
        timeStep *= 0.2;
    }
    if (currentTime + timeStep > finalTime)
    {
        timeStep = finalTime - currentTime;
    }

    return timeStep;
}

// Computes the estimated maximum wave speed across the entire computational domain for a given multimaterial system (i.e. the maximum wave speed within either material),
// for the purposes of timestep calculcation.
double Solvers::computeGhostFluidMaximumWaveSpeed(MultimaterialSystem multimaterialSystem)
{
    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();

    double maximumWaveSpeed = 0.0;
    int cellCount = leftMaterialCells.size();

    for (int i = 0; i < cellCount; i++)
    {
        double waveSpeed = max(abs(leftMaterialCells[i].getXVelocity()), abs(rightMaterialCells[i].getXVelocity())) + max(leftMaterialCells[i].computeSoundSpeed(),
                                                                                                                          rightMaterialCells[i].computeSoundSpeed());
        if (waveSpeed > maximumWaveSpeed)
        {
            maximumWaveSpeed = waveSpeed;
        }
    }

    return maximumWaveSpeed;
}

// Computes the maximum stable timestep for a given multimaterial system (i.e. using the maximum wave speed within either material), in accordance with the CFL condition.
double Solvers::computeGhostFluidStableTimeStep(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                                int currentIteration)
{
    double timeStep = CFLCoefficient * (cellSpacing / computeGhostFluidMaximumWaveSpeed(multimaterialSystem));

    if (currentIteration <= 5)
    {
        timeStep *= 0.2;
    }
    if (currentTime + timeStep > finalTime)
    {
        timeStep = finalTime - currentTime;
    }

    return timeStep;
}

// Evolves the staet given by middleStateVector for half a timestep, using the slope across the neighbouring cells, given by leftStateVector, middleStateVector and rightStateVector,
// respectively.
StateVector Solvers::computeEvolvedState(StateVector leftStateVector, StateVector middleStateVector, StateVector rightStateVector, double cellSpacing, double timeStep, double bias,
                                         int limiter, int side)
{
    vector<double> slopeVector = SlopeLimiter::computeSlopeVector(leftStateVector, middleStateVector, rightStateVector, bias, limiter);

    vector<double> leftBoundaryExtrapolatedValue = LinearAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(), LinearAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> rightBoundaryExtrapolatedValue = LinearAlgebra::addVectors(middleStateVector.computeConservedVariableVector(), LinearAlgebra::multiplyVector(0.5, slopeVector));
    vector<double> fluxUpdate = LinearAlgebra::multiplyVector(0.5 * (timeStep / cellSpacing), LinearAlgebra::subtractVectors(
                                                                  StateVector::computeFluxVector(leftBoundaryExtrapolatedValue), StateVector::computeFluxVector(rightBoundaryExtrapolatedValue)));

    StateVector evolvedStateVector;

    if (side == 0)
    {
        evolvedStateVector.setConservedVariableVector(LinearAlgebra::addVectors(leftBoundaryExtrapolatedValue, fluxUpdate));
    }
    else
    {
        evolvedStateVector.setConservedVariableVector(LinearAlgebra::addVectors(rightBoundaryExtrapolatedValue, fluxUpdate));
    }

    return evolvedStateVector;
}
