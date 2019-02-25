#include "ghostfluidmethod.h"

// This class encapsulates the methods that are common to all ghost fluid solvers for the multimaterial (stiffened gas) Euler equations.
GhostFluidMethod::GhostFluidMethod()
{
}

// Returns the computational domain, but with Fedkiw's original ghost fluid boundary conditions applied.
MultimaterialSystem GhostFluidMethod::applyBasicGhostFluidBoundaryConditions(MultimaterialSystem multimaterialSystem)
{
    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<StateVector> newLeftMaterialCells = leftMaterialCells;
    vector<StateVector> newRightMaterialCells = rightMaterialCells;

    double oldLevelSetValue = levelSetFunction[0];
    int cellCount = levelSetFunction.size();
    bool ghostRegion = false;

    for (int i = 1; i < cellCount; i++)
    {
        if ((levelSetFunction[i] * oldLevelSetValue) <= 0.0 && i > 2 && i < (cellCount - 3))
        {
            for (int j = 0; j < 3; j++)
            {
                if (!ghostRegion)
                {
                    double leftGhostPressure = rightMaterialCells[i + j].getPressure();
                    double leftInterfacePressure = leftMaterialCells[i - 1].getPressure();
                    double leftInterfaceAdiabaticIndex = leftMaterialCells[i - 1].getAdiabaticIndex();
                    double leftInterfaceDensity = leftMaterialCells[i - 1].getDensity();
                    double leftGhostVelocity = rightMaterialCells[i + j].getXVelocity();
                    double leftGhostDensity = pow((leftGhostPressure / leftInterfacePressure), 1.0 / leftInterfaceAdiabaticIndex) * leftInterfaceDensity;

                    newLeftMaterialCells[i + j].setPressure(leftGhostPressure);
                    newLeftMaterialCells[i + j].setXVelocity(leftGhostVelocity);
                    newLeftMaterialCells[i + j].setDensity(leftGhostDensity);

                    double rightGhostPressure = leftMaterialCells[i - j - 1].getPressure();
                    double rightInterfacePressure = rightMaterialCells[i].getPressure();
                    double rightInterfaceAdiabaticIndex = rightMaterialCells[i].getAdiabaticIndex();
                    double rightInterfaceDensity = rightMaterialCells[i].getDensity();
                    double rightGhostVelocity = leftMaterialCells[i - j - 1].getXVelocity();
                    double rightGhostDensity = pow((rightGhostPressure / rightInterfacePressure), 1.0 / rightInterfaceAdiabaticIndex) * rightInterfaceDensity;

                    newRightMaterialCells[i - j - 1].setPressure(rightGhostPressure);
                    newRightMaterialCells[i - j - 1].setXVelocity(rightGhostVelocity);
                    newRightMaterialCells[i - j - 1].setDensity(rightGhostDensity);
                }
                else
                {
                    double rightGhostPressure = leftMaterialCells[i + j].getPressure();
                    double rightInterfacePressure = rightMaterialCells[i - 1].getPressure();
                    double rightInterfaceAdiabaticIndex = rightMaterialCells[i - 1].getAdiabaticIndex();
                    double rightInterfaceDensity = rightMaterialCells[i - 1].getDensity();
                    double rightGhostVelocity = leftMaterialCells[i + j].getXVelocity();
                    double rightGhostDensity = pow((rightGhostPressure / rightInterfacePressure), 1.0 / rightInterfaceAdiabaticIndex) * rightInterfaceDensity;

                    newRightMaterialCells[i + j].setPressure(rightGhostPressure);
                    newRightMaterialCells[i + j].setXVelocity(rightGhostVelocity);
                    newRightMaterialCells[i + j].setDensity(rightGhostDensity);

                    double leftGhostPressure = rightMaterialCells[i - j - 1].getPressure();
                    double leftInterfacePressure = leftMaterialCells[i].getPressure();
                    double leftInterfaceAdiabaticIndex = leftMaterialCells[i].getAdiabaticIndex();
                    double leftInterfaceDensity = leftMaterialCells[i].getDensity();
                    double leftGhostVelocity = rightMaterialCells[i - j - 1].getXVelocity();
                    double leftGhostDensity = pow((leftGhostPressure / leftInterfacePressure), 1.0 / leftInterfaceAdiabaticIndex) * leftInterfaceDensity;

                    newLeftMaterialCells[i - j - 1].setPressure(leftGhostPressure);
                    newLeftMaterialCells[i - j - 1].setXVelocity(leftGhostVelocity);
                    newLeftMaterialCells[i - j - 1].setDensity(leftGhostDensity);
                }
            }

            ghostRegion = !ghostRegion;
        }

        oldLevelSetValue = levelSetFunction[i];
    }

    return MultimaterialSystem(newLeftMaterialCells, newRightMaterialCells, levelSetFunction);
}

// Returns the computational domain, but with the real (i.e. Riemann problem-based) ghost fluid boundary conditions applied.
MultimaterialSystem GhostFluidMethod::applyRiemannGhostFluidBoundaryConditions(MultimaterialSystem multimaterialSystem)
{
    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<StateVector> newLeftMaterialCells = leftMaterialCells;
    vector<StateVector> newRightMaterialCells = rightMaterialCells;

    double oldLevelSetValue = levelSetFunction[0];
    int cellCount = levelSetFunction.size();
    bool ghostRegion = false;

    for (int i = 1; i < cellCount; i++)
    {
        if ((levelSetFunction[i] * oldLevelSetValue) <= 0.0 && i > 2 && i < (cellCount - 3))
        {
            for (int j = 0; j < 4; j++)
            {
                if (!ghostRegion)
                {
                    StateVector riemannProblemSolution = ExactSolver::solve(0.0, 1.0, 0.0, leftMaterialCells[i - 1], rightMaterialCells[i]);

                    double leftGhostPressure = riemannProblemSolution.getPressure();
                    double leftGhostVelocity = riemannProblemSolution.getXVelocity();
                    double leftGhostDensity = ExactSolver::computeStarRegionDensity(leftGhostPressure, leftMaterialCells[i - 1]);

                    newLeftMaterialCells[i + j - 1].setPressure(leftGhostPressure);
                    newLeftMaterialCells[i + j - 1].setXVelocity(leftGhostVelocity);
                    newLeftMaterialCells[i + j - 1].setDensity(leftGhostDensity);

                    double rightGhostPressure = riemannProblemSolution.getPressure();
                    double rightGhostVelocity = riemannProblemSolution.getXVelocity();
                    double rightGhostDensity = ExactSolver::computeStarRegionDensity(rightGhostPressure, rightMaterialCells[i]);

                    newRightMaterialCells[i - j].setPressure(rightGhostPressure);
                    newRightMaterialCells[i - j].setXVelocity(rightGhostVelocity);
                    newRightMaterialCells[i - j].setDensity(rightGhostDensity);
                }
                else
                {
                    StateVector riemannProblemSolution = ExactSolver::solve(0.0, 1.0, 0.0, rightMaterialCells[i - 1], leftMaterialCells[i]);

                    double rightGhostPressure = riemannProblemSolution.getPressure();
                    double rightGhostVelocity = riemannProblemSolution.getXVelocity();
                    double rightGhostDensity = ExactSolver::computeStarRegionDensity(rightGhostPressure, rightMaterialCells[i - 1]);

                    newRightMaterialCells[i + j - 1].setPressure(rightGhostPressure);
                    newRightMaterialCells[i + j - 1].setXVelocity(rightGhostVelocity);
                    newRightMaterialCells[i + j - 1].setDensity(rightGhostDensity);

                    double leftGhostPressure = riemannProblemSolution.getPressure();
                    double leftGhostVelocity = riemannProblemSolution.getXVelocity();
                    double leftGhostDensity = ExactSolver::computeStarRegionDensity(leftGhostPressure, leftMaterialCells[i]);

                    newLeftMaterialCells[i - j].setPressure(leftGhostPressure);
                    newLeftMaterialCells[i - j].setXVelocity(leftGhostVelocity);
                    newLeftMaterialCells[i - j].setDensity(leftGhostDensity);
                }
            }

            ghostRegion = !ghostRegion;
        }

        oldLevelSetValue = levelSetFunction[i];
    }

    return MultimaterialSystem(newLeftMaterialCells, newRightMaterialCells, levelSetFunction);
}

/*
// Returns the computational domain, but with the real (i.e. Riemann problem-based) ghost fluid boundary conditions applied.
MultimaterialSystem GhostFluidMethod::applyRiemannGhostFluidBoundaryConditions(MultimaterialSystem multimaterialSystem, double cellSpacing, double timeStep)
{
    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    vector<StateVector> newLeftMaterialCells = leftMaterialCells;
    vector<StateVector> newRightMaterialCells = rightMaterialCells;

    double oldLevelSetValue = levelSetFunction[0];
    int cellCount = levelSetFunction.size();
    bool ghostRegion = false;

    for (int i = 1; i < cellCount; i++)
    {
        if ((levelSetFunction[i] * oldLevelSetValue) <= 0 && i > 3 && i < (cellCount - 4))
        {
            for (int j = 0; j < 3; j++)
            {
                if (!ghostRegion)
                {
                    StateVector leftRiemannProblemSolution = ExactSolver::solve((cellSpacing / 2.0), timeStep, 0.0, leftMaterialCells[i - 1], rightMaterialCells[i]);
                    //StateVector leftRiemannProblemSolution = ExactSolver::solve(cellSpacing, timeStep, 0.0, leftMaterialCells[i - j], rightMaterialCells[i + j - 1]);

                    newLeftMaterialCells[i + j].setPressure(leftRiemannProblemSolution.getPressure());
                    newLeftMaterialCells[i + j].setXVelocity(leftRiemannProblemSolution.getXVelocity());
                    newLeftMaterialCells[i + j].setDensity(leftRiemannProblemSolution.getDensity());

                    StateVector rightRiemannProblemSolution = ExactSolver::solve(-(cellSpacing / 2.0), timeStep, 0.0, leftMaterialCells[i - 1], rightMaterialCells[i]);
                    //StateVector rightRiemannProblemSolution = ExactSolver::solve(-cellSpacing, timeStep, 0.0, leftMaterialCells[i - j], rightMaterialCells[i + j - 1]);

                    newRightMaterialCells[i - j - 1].setPressure(rightRiemannProblemSolution.getPressure());
                    newRightMaterialCells[i - j - 1].setXVelocity(rightRiemannProblemSolution.getXVelocity());
                    newRightMaterialCells[i - j - 1].setDensity(rightRiemannProblemSolution.getDensity());
                }
                else
                {
                    StateVector rightRiemannProblemSolution = ExactSolver::solve((cellSpacing / 2.0), timeStep, 0.0, rightMaterialCells[i - 1], leftMaterialCells[i]);
                    //StateVector rightRiemannProblemSolution = ExactSolver::solve(cellSpacing, timeStep, 0.0, rightMaterialCells[i - j], leftMaterialCells[i + j - 1]);

                    newRightMaterialCells[i + j].setPressure(rightRiemannProblemSolution.getPressure());
                    newRightMaterialCells[i + j].setXVelocity(rightRiemannProblemSolution.getXVelocity());
                    newRightMaterialCells[i + j].setDensity(rightRiemannProblemSolution.getDensity());

                    StateVector leftRiemannProblemSolution = ExactSolver::solve(-(cellSpacing / 2.0), timeStep, 0.0, rightMaterialCells[i - 1], leftMaterialCells[i]);
                    //StateVector leftRiemannProblemSolution = ExactSolver::solve(-cellSpacing, timeStep, 0.0, rightMaterialCells[i - j], leftMaterialCells[i + j - 1]);

                    newLeftMaterialCells[i - j - 1].setPressure(leftRiemannProblemSolution.getPressure());
                    newLeftMaterialCells[i - j - 1].setXVelocity(leftRiemannProblemSolution.getXVelocity());
                    newLeftMaterialCells[i - j - 1].setDensity(leftRiemannProblemSolution.getDensity());
                }
            }

            ghostRegion = !ghostRegion;
            if (levelSetFunction[i] * levelSetFunction[i + 1] <= 0)
            {
                oldLevelSetValue = levelSetFunction[i + 1];
            }
            else
            {
                oldLevelSetValue = levelSetFunction[i];
            }
        }
    }

    return MultimaterialSystem(newLeftMaterialCells, newRightMaterialCells, levelSetFunction);
}
*/

// Evolves the level set function, using the standard first-order upwind scheme.
vector<double> GhostFluidMethod::updateLevelSetFunction(vector<double> oldLevelSetFunction, vector<StateVector> leftMaterialCells, vector<StateVector> rightMaterialCells,
                                                        double cellSpacing, double timeStep)
{
    int cellCount = oldLevelSetFunction.size();
    vector<double> newLevelSetFunction = oldLevelSetFunction;
    double upwindApproximation;

    double velocity = leftMaterialCells[0].getXVelocity();
    if (velocity < 0)
    {
        upwindApproximation = (oldLevelSetFunction[1] - oldLevelSetFunction[0]) / cellSpacing;
    }
    else
    {
        upwindApproximation = (oldLevelSetFunction[0] - oldLevelSetFunction[0]) / cellSpacing;
    }
    newLevelSetFunction[0] = oldLevelSetFunction[0] - (timeStep * (velocity * upwindApproximation));

    int oldLevelSetValue = oldLevelSetFunction[0];
    bool ghostRegion = false;

    for (int i = 1; i < cellCount - 1; i++)
    {
        if ((oldLevelSetValue * oldLevelSetFunction[i]) <= 0.0)
        {
            ghostRegion = !ghostRegion;
        }

        if (!ghostRegion)
        {
            velocity = leftMaterialCells[i].getXVelocity();
        }
        else
        {
            velocity = rightMaterialCells[i].getXVelocity();
        }

        if (velocity < 0)
        {
            upwindApproximation = (oldLevelSetFunction[i + 1] - oldLevelSetFunction[i]) / cellSpacing;
        }
        else
        {
            upwindApproximation = (oldLevelSetFunction[i] - oldLevelSetFunction[i - 1]) / cellSpacing;
        }
        newLevelSetFunction[i] = oldLevelSetFunction[i] - (timeStep * (velocity * upwindApproximation));

        oldLevelSetValue = oldLevelSetFunction[i];
    }

    if (!ghostRegion)
    {
        velocity = leftMaterialCells[cellCount - 1].getXVelocity();
    }
    else
    {
        velocity = rightMaterialCells[cellCount - 1].getXVelocity();
    }

    if (velocity < 0)
    {
        upwindApproximation = (oldLevelSetFunction[cellCount - 1] - oldLevelSetFunction[cellCount - 1]) / cellSpacing;
    }
    else
    {
        upwindApproximation = (oldLevelSetFunction[cellCount - 1] - oldLevelSetFunction[cellCount - 2]) / cellSpacing;
    }
    newLevelSetFunction[cellCount - 1] = oldLevelSetFunction[cellCount - 1] - (timeStep * (velocity * upwindApproximation));

    return newLevelSetFunction;
}

// Computes the number of distinct interfaces, using the cell-centred level set function.
int GhostFluidMethod::computeInterfaceCount(vector<double> levelSetFunction)
{
    int cellCount = levelSetFunction.size();
    int interfaceCount = 0;
    double oldLevelSetValue = levelSetFunction[0];

    for (int i = 1; i < cellCount; i++)
    {
        if ((levelSetFunction[i] * oldLevelSetValue) <= 0.0 && (levelSetFunction[i - 1] * oldLevelSetValue) > 0.0)
        {
            interfaceCount += 1;
        }

        oldLevelSetValue = levelSetFunction[i];
    }

    return interfaceCount;
}

// Reinitialises a cell-centred level set function, assuming a single interface.
vector<double> GhostFluidMethod::reinitialiseLevelSetFunction(vector<double> oldLevelSetFunction)
{
    int cellCount = oldLevelSetFunction.size();
    double cellSpacing = (1.0 / cellCount);
    vector<double> newLevelSetFunction(cellCount);
    double oldLevelSetValue = oldLevelSetFunction[0];
    double interfaceLocation;

    for (int i = 1; i < cellCount; i++)
    {
        if ((oldLevelSetFunction[i] * oldLevelSetValue) <= 0.0  && (oldLevelSetFunction[i - 1] * oldLevelSetValue) > 0.0)
        {
            interfaceLocation = cellSpacing * i;
        }

        oldLevelSetValue = oldLevelSetFunction[i];
    }

    for (int i = 0; i < cellCount; i++)
    {
        newLevelSetFunction[i] = (cellSpacing * i) - interfaceLocation;
    }

    return newLevelSetFunction;
}

// Reinitialises a cell-centred level set funtion, assuming two distinct interfaces.
vector<double> GhostFluidMethod::reinitialiseTwoInterfaceLevelSetFunction(vector<double> oldLevelSetFunction)
{
    int cellCount = oldLevelSetFunction.size();
    double cellSpacing = (1.0 / cellCount);
    vector<double> newLevelSetFunction(cellCount);
    double oldLevelSetValue = oldLevelSetFunction[0];
    double interfaceLocation1;
    double interfaceLocation2;
    int interfaceCount = 0;

    for (int i = 1; i < cellCount; i++)
    {
        if ((oldLevelSetFunction[i] * oldLevelSetValue) <= 0.0 && (oldLevelSetFunction[i - 1] * oldLevelSetValue) > 0.0)
        {
            interfaceCount += 1;

            if (interfaceCount == 1)
            {
                interfaceLocation1 = cellSpacing * i;
            }
            else
            {
                interfaceLocation2 = cellSpacing * i;
            }
        }

        oldLevelSetValue = oldLevelSetFunction[i];
    }

    for (int i = 0; i < cellCount; i++)
    {
        if ((cellSpacing * i) <= (0.5 * (interfaceLocation1 + interfaceLocation2)))
        {
            newLevelSetFunction[i] = (cellSpacing * i) - interfaceLocation1;
        }
        else
        {
            newLevelSetFunction[i] = -(cellSpacing * i) + interfaceLocation2;
        }
    }

    return newLevelSetFunction;
}
