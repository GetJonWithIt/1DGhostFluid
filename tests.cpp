#include "tests.h"

// This class encapsulates the relevant single-material and multimaterial tests for the exact solver, the original ghost fluid method and the real (Riemann problem-based) ghost fluid method.
Tests::Tests()
{
}

// Solves Toro's first test using the exact solver, at the given resolution.
void Tests::solveToroTest1Exact(int cellCount)
{
    outputExactSolution(cellCount, 0.25, 0.5, StateVector(1.0, 0.0, 0.0, 0.0, 1.0, 1.4, 0.0), StateVector(0.125, 0.0, 0.0, 0.0, 0.1, 1.4, 0.0));
}

// Solves Toro's second test using the exact solver, at the given resolution.
void Tests::solveToroTest2Exact(int cellCount)
{
    outputExactSolution(cellCount, 0.15, 0.5, StateVector(1.0, -2.0, 0.0, 0.0, 0.4, 1.4, 0.0), StateVector(1.0, 2.0, 0.0, 0.0, 0.4, 1.4, 0.0));
}

// Solves Toro's third test using the exact solver, at the given resolution.
void Tests::solveToroTest3Exact(int cellCount)
{
    outputExactSolution(cellCount, 0.012, 0.5, StateVector(1.0, 0.0, 0.0, 0.0, 1000.0, 1.4, 0.0), StateVector(1.0, 0.0, 0.0, 0.0, 0.01, 1.4, 0.0));
}

// Solves Toro's fourth test using the exact solver, at the given resolution.
void Tests::solveToroTest4Exact(int cellCount)
{
    outputExactSolution(cellCount, 0.035, 0.5, StateVector(1.0, 0.0, 0.0, 0.0, 0.01, 1.4, 0.0), StateVector(1.0, 0.0, 0.0, 0.0, 100.0, 1.4, 0.0));
}

// Solves Toro's fifth test using the exact solver, at the given resolution.
void Tests::solveToroTest5Exact(int cellCount)
{
    outputExactSolution(cellCount, 0.035, 0.5, StateVector(5.99924, 19.5975, 0.0, 0.0, 460.894, 1.4, 0.0), StateVector(5.99242, -6.19633, 0.0, 0.0, 46.0950, 1.4, 0.0));
}

// Solves the test from the paper of Chinnayya et al. using the exact solver, at the given resolution.
void Tests::solveChinnayyaTestExact(int cellCount)
{
    outputExactSolution(cellCount, 237.44 * pow(10, -6), 0.7, StateVector(1000.0, 0.0, 0.0, 0.0, pow(10, 9), 4.4, 6 * pow(10, 8)),
                        StateVector(50.0, 0.0, 0.0, 0.0, pow(10, 5), 1.4, 0.0));
}

// Solves the test from the paper of Fedkiw et al. using the exact solver, at the given resolution.
void Tests::solveFedkiwTestExact(int cellCount)
{
    outputExactSolution(cellCount, 0.0002, 0.5, StateVector(1.3333, 0.3535 * sqrt(pow(10, 5)), 0.0, 0.0, 1.5 * pow(10, 5), 1.4, 0.0),
                        StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0));
}

// Solves Toro's first test using the original ghost fluid method, at the given resolution.
void Tests::solveToroTest1BasicGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(0.125, 0.0, 0.0, 0.0, 0.1, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.25, 0.0, 0));
}

// Solves Toro's first test using the original ghost fluid method with reinitialisation, with the given resolution and reinitialisation period.
void Tests::solveToroTest1BasicGhostFluidReinitialisation(int cellCount, int reinitialisationIteration)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(0.125, 0.0, 0.0, 0.0, 0.1, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluidReinitialisation(
                                MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.25, 0.0, 0, reinitialisationIteration));
}

// Solves Toro's second test using the original ghost fluid method, at the given resolution (note: this fails).
void Tests::solveToroTest2BasicGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, -2.0, 0.0, 0.0, 0.4, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(1.0, 2.0, 0.0, 0.0, 0.4, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.15, 0.0, 0));
}

// Solves Toro's third test using the original ghost fluid method, at the given resolution (note: this fails).
void Tests::solveToroTest3BasicGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1000.0, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 0.01, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.012, 0.0, 0));
}

// Solves Toro's first test using the real (Riemann problem-based) ghost fluid method, at the given resolution.
void Tests::solveToroTest1RiemannGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(0.125, 0.0, 0.0, 0.0, 0.1, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.25, 0.0, 0));
}

// Solves Toro's first test using the real (Riemann problem-based) ghost fluid method with reinitialisation, with the given resolution and reinitialisation period.
void Tests::solveToroTest1RiemannGhostFluidReinitialisation(int cellCount, int reinitialisationIteration)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(0.125, 0.0, 0.0, 0.0, 0.1, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluidReinitialisation(
                                MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.25, 0.0, 0.15, reinitialisationIteration));
}

// Solves Toro's second test using the real (Riemann problem-based) ghost fluid method, at the given resolution.
void Tests::solveToroTest2RiemannGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, -2.0, 0.0, 0.0, 0.4, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(1.0, 2.0, 0.0, 0.0, 0.4, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.15, 0.0, 0));
}

// Solves Toro's third test using the real (Riemann problem-based) ghost fluid method, at the given resolution.
void Tests::solveToroTest3RiemannGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1000.0, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 0.01, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.012, 0.0, 0));
}

// Solves Toro's third test using the real (Riemann problem-based) ghost fluid method with reinitialisation, with the given resolution and reinitialisation period.
void Tests::solveToroTest3RiemannGhostFluidReinitialisation(int cellCount, int reinitialisationIteration)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1000.0, 1.4, 0.0);
        rightMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 0.01, 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluidReinitialisation(
                                MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.012, 0.0, 0, reinitialisationIteration));
}

// Solves the test from the paper of Fedkiw et al. using the original ghost fluid method, at the given resolution.
void Tests::solveFedkiwTestBasicGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            leftMaterialCells[i] = StateVector(1.3333, 0.3535 * sqrt(pow(10, 5)), 0.0, 0.0, 1.5 * pow(10, 5), 1.4, 0.0);
        }
        else
        {
            leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.4, 0.0);
        }
        rightMaterialCells[i] = StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.0012, 0.0, 0));
}

// Solves the test from the paper of Wang et al. using the original ghost fluid method, at the given resolution.
void Tests::solveWangTestBasicGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            leftMaterialCells[i] = StateVector(1.3333, 0.3535 * sqrt(pow(10, 5)), 0.0, 0.0, 1.5 * pow(10, 5), 1.4, 0.0);
        }
        else
        {
            leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.4, 0.0);
        }
        rightMaterialCells[i] = StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0);

        if (i <= (cellCount / 2))
        {
            levelSetFunction[i] = (cellSpacing * i) - 0.4;
        }
        else
        {
            levelSetFunction[i] = -(cellSpacing * i) + 0.6;
        }
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.0014, 0.0, 0));
}

// Solves the Mach 10 test using the original ghost fluid method, at the given resolution.
void Tests::solveMach10TestBasicGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            leftMaterialCells[i] = StateVector(5.92593, 6220.51, 0.0, 0.0, 4.665 * pow(10, 7), 1.4, 0.0);
        }
        else
        {
            leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.4, 0.0);
        }
        rightMaterialCells[i] = StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0);

        if (i <= (cellCount / 2))
        {
            levelSetFunction[i] = (cellSpacing * i) - 0.4;
        }
        else
        {
            levelSetFunction[i] = -(cellSpacing * i) + 0.6;
        }
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.5, 0.0002, 0.0, 0));
}

// Solves the Mach 10 test using the original ghost fluid method with reinitialisation, with the given resolution and reinitialisation period.
void Tests::solveMach10TestBasicGhostFluidReinitialisation(int cellCount, int reinitialisationIteration)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            leftMaterialCells[i] = StateVector(5.92593, 6220.51, 0.0, 0.0, 4.665 * pow(10, 7), 1.4, 0.0);
        }
        else
        {
            leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.4, 0.0);
        }
        rightMaterialCells[i] = StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0);

        if (i <= (cellCount / 2))
        {
            levelSetFunction[i] = (cellSpacing * i) - 0.4;
        }
        else
        {
            levelSetFunction[i] = -(cellSpacing * i) + 0.6;
        }
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluidReinitialisation(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.5, 0.0002, 0.0,
                                                                                  0, reinitialisationIteration));
}

// Solves the test from the paper of Fedkiw et al. using the real (Rieman problem-based) ghost fluid method, at the given resolution.
void Tests::solveFedkiwTestRiemannGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            leftMaterialCells[i] = StateVector(1.3333, 0.3535 * sqrt(pow(10, 5)), 0.0, 0.0, 1.5 * pow(10, 5), 1.4, 0.0);
        }
        else
        {
            leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.4, 0.0);
        }
        rightMaterialCells[i] = StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.5;
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.0012, 0.0, 0));
}

// Solves the test from the paper of Wang et al. using the real (Riemann problem-based) ghost fluid method, at the given resolution.
void Tests::solveWangTestRiemannGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            leftMaterialCells[i] = StateVector(1.3333, 0.3535 * sqrt(pow(10, 5)), 0.0, 0.0, 1.5 * pow(10, 5), 1.4, 0.0);
        }
        else
        {
            leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.4, 0.0);
        }
        rightMaterialCells[i] = StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0);

        if (i <= (cellCount / 2))
        {
            levelSetFunction[i] = (cellSpacing * i) - 0.4;
        }
        else
        {
            levelSetFunction[i] = -(cellSpacing * i) + 0.6;
        }
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.0014, 0.0, 0));
}

// Solves the Mach 10 test using the real (Riemann problem-based) ghost fluid method, at the given resolution.
void Tests::solveMach10TestRiemannGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            leftMaterialCells[i] = StateVector(5.92593, 6220.51, 0.0, 0.0, 4.665 * pow(10, 7), 1.4, 0.0);
        }
        else
        {
            leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.4, 0.0);
        }
        rightMaterialCells[i] = StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0);

        if (i <= (cellCount / 2.0))
        {
            levelSetFunction[i] = (cellSpacing * i) - 0.4;
        }
        else
        {
            levelSetFunction[i] = -(cellSpacing * i) + 0.6;
        }
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.0002, 0.0, 0));
}

// Solves the Mach 10 test using the real (Riemann problem-based) ghost fluid method with reinitialisation, with the given resolution and reinitialisation period.
void Tests::solveMach10TestRiemannGhostFluidReinitialisation(int cellCount, int reinitialisationIteration)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        if (i <= (cellCount * 0.05))
        {
            leftMaterialCells[i] = StateVector(5.92593, 6220.51, 0.0, 0.0, 4.665 * pow(10, 7), 1.4, 0.0);
        }
        else
        {
            leftMaterialCells[i] = StateVector(1.0, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.4, 0.0);
        }
        rightMaterialCells[i] = StateVector(0.1379, 0.0, 0.0, 0.0, 1.0 * pow(10, 5), 1.67, 0.0);

        if (i <= (cellCount / 2.0))
        {
            levelSetFunction[i] = (cellSpacing * i) - 0.4;
        }
        else
        {
            levelSetFunction[i] = -(cellSpacing * i) + 0.6;
        }
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluidReinitialisation(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 0.0002, 0.0,
                                                                                    0, reinitialisationIteration));
}

// Solves the test from the paper of Chinnayya et al. using the original ghost fluid method, at the given resolution (note: this fails).
void Tests::solveChinnayyaTestBasicGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1000.0, 0.0, 0.0, 0.0, pow(10, 9), 4.4, 6.0 * pow(10, 8));
        rightMaterialCells[i] = StateVector(50.0, 0.0, 0.0, 0.0, pow(10, 5), 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.7;
    }

    outputMultimaterialData(MUSCLSolver::solveWithBasicGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 237.44 * pow(10, -6), 0.0, 0));
}

// Solves the test from the paper of Chinnayya et al. using the real (Riemann problem-based) ghost fluid method, at the given resolution.
void Tests::solveChinnayyaTestRiemannGhostFluid(int cellCount)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1000.0, 0.0, 0.0, 0.0, pow(10, 9), 4.4, 6.0 * pow(10, 8));
        rightMaterialCells[i] = StateVector(50.0, 0.0, 0.0, 0.0, pow(10, 5), 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.7;
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluid(MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 237.44 * pow(10, -6), 0.0, 0));
}

// Solves the test from the paper of Chinnayya et al. using the real (Riemann problem-based) ghost fluid method with reinitialisation, with the given resolution and reinitialisation period.
void Tests::solveChinnayyaTestRiemannGhostFluidReinitialisation(int cellCount, int reinitialisationIteration)
{
    vector<StateVector> leftMaterialCells(cellCount);
    vector<StateVector> rightMaterialCells(cellCount);
    vector<double> levelSetFunction(cellCount);
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        leftMaterialCells[i] = StateVector(1000.0, 0.0, 0.0, 0.0, pow(10, 9), 4.4, 6.0 * pow(10, 8));
        rightMaterialCells[i] = StateVector(50.0, 0.0, 0.0, 0.0, pow(10, 5), 1.4, 0.0);

        levelSetFunction[i] = (cellSpacing * i) - 0.7;
    }

    outputMultimaterialData(MUSCLSolver::solveWithRiemannGhostFluidReinitialisation(
                                MultimaterialSystem(leftMaterialCells, rightMaterialCells, levelSetFunction), cellSpacing, 0.8, 237.44 * pow(10, -6), 0.0, 0, reinitialisationIteration));
}

// Produces output files of density, velocity, pressure and specific internal energy, for a given exact Riemann problem solution.
void Tests::outputExactSolution(int cellCount, double finalTime, double interfaceLocation, StateVector leftStateVector, StateVector rightStateVector)
{
    ofstream densityFile("density.dat");
    ofstream velocityFile("velocity.dat");
    ofstream pressureFile("pressure.dat");
    ofstream energyFile("energy.dat");

    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        StateVector exactSolution = ExactSolver::solve(cellSpacing * i, finalTime, interfaceLocation, leftStateVector, rightStateVector);

        densityFile << (cellSpacing * i) << " " << exactSolution.getDensity() << endl;
        velocityFile << (cellSpacing * i) << " " << exactSolution.getXVelocity() << endl;
        pressureFile << (cellSpacing * i) << " " << exactSolution.getPressure() << endl;
        energyFile << (cellSpacing * i) << " " << exactSolution.computeInternalEnergy() << endl;
    }

    densityFile.close();
    velocityFile.close();
    pressureFile.close();
    energyFile.close();
}

// Produces output files of density, velocity, pressure and specific internal energy, for a given solution to a single-material Riemann problem.
void Tests::outputSingleMaterialData(vector<StateVector> cells)
{
    ofstream densityFile("density.dat");
    ofstream velocityFile("velocity.dat");
    ofstream pressureFile("pressure.dat");
    ofstream energyFile("energy.dat");
    int cellCount = cells.size();
    double cellSpacing = 1.0 / cellCount;

    for (int i = 0; i < cellCount; i++)
    {
        densityFile << (cellSpacing * i) << " " << cells[i].getDensity() << endl;
        velocityFile << (cellSpacing * i) << " " << cells[i].getXVelocity() << endl;
        pressureFile << (cellSpacing * i) << " " << cells[i].getPressure() << endl;
        energyFile << (cellSpacing * i) << " " << cells[i].computeInternalEnergy() << endl;
    }

    densityFile.close();
    velocityFile.close();
    pressureFile.close();
    energyFile.close();
}

// Produces output files of density, velocity, pressure and specific internal energy, for a given solution to a multimaterial Riemann problem.
void Tests::outputMultimaterialData(MultimaterialSystem multimaterialSystem)
{
    ofstream leftDensityFile("leftDensity.dat");
    ofstream leftVelocityFile("leftVelocity.dat");
    ofstream leftPressureFile("leftPressure.dat");
    ofstream leftEnergyFile("leftEnergy.dat");

    ofstream rightDensityFile("rightDensity.dat");
    ofstream rightVelocityFile("rightVelocity.dat");
    ofstream rightPressureFile("rightPressure.dat");
    ofstream rightEnergyFile("rightEnergy.dat");

    ofstream densityFile("density.dat");
    ofstream velocityFile("velocity.dat");
    ofstream pressureFile("pressure.dat");
    ofstream energyFile("energy.dat");

    ofstream levelSetFile("levelSet.dat");

    vector<StateVector> leftMaterialCells = multimaterialSystem.getLeftMaterialCells();
    vector<StateVector> rightMaterialCells = multimaterialSystem.getRightMaterialCells();
    vector<double> levelSetFunction = multimaterialSystem.getLevelSetFunction();

    int cellCount = leftMaterialCells.size();
    double cellSpacing = 1.0 / cellCount;
    bool ghostRegion = false;
    double oldLevelSetValue = levelSetFunction[0];

    for (int i = 0; i < cellCount; i++)
    {
        leftDensityFile << (cellSpacing * i) << " " << leftMaterialCells[i].getDensity() << endl;
        leftVelocityFile << (cellSpacing * i) << " " << leftMaterialCells[i].getXVelocity() << endl;
        leftPressureFile << (cellSpacing * i) << " " << leftMaterialCells[i].getPressure() << endl;
        leftEnergyFile << (cellSpacing * i) << " " << leftMaterialCells[i].computeInternalEnergy() << endl;

        rightDensityFile << (cellSpacing * i) << " " << rightMaterialCells[i].getDensity() << endl;
        rightVelocityFile << (cellSpacing * i) << " " << rightMaterialCells[i].getXVelocity() << endl;
        rightPressureFile << (cellSpacing * i) << " " << rightMaterialCells[i].getPressure() << endl;
        rightEnergyFile << (cellSpacing * i) << " " << rightMaterialCells[i].computeInternalEnergy() << endl;

        if ((levelSetFunction[i] * oldLevelSetValue) <= 0.0 && (levelSetFunction[i - 1] * oldLevelSetValue) > 0.0)
        {
            ghostRegion = !ghostRegion;
        }

        if (!ghostRegion)
        {
            densityFile << (cellSpacing * i) << " " << leftMaterialCells[i].getDensity() << endl;
            velocityFile << (cellSpacing * i) << " " << leftMaterialCells[i].getXVelocity() << endl;
            pressureFile << (cellSpacing * i) << " " << leftMaterialCells[i].getPressure() << endl;
            energyFile << (cellSpacing * i) << " " << leftMaterialCells[i].computeInternalEnergy() << endl;
        }
        else
        {
            densityFile << (cellSpacing * i) << " " << rightMaterialCells[i].getDensity() << endl;
            velocityFile << (cellSpacing * i) << " " << rightMaterialCells[i].getXVelocity() << endl;
            pressureFile << (cellSpacing * i) << " " << rightMaterialCells[i].getPressure() << endl;
            energyFile << (cellSpacing * i) << " " << rightMaterialCells[i].computeInternalEnergy() << endl;
        }

        levelSetFile << (cellSpacing * i) << " " << levelSetFunction[i] << endl;
        oldLevelSetValue = levelSetFunction[i];
    }

    leftDensityFile.close();
    leftVelocityFile.close();
    leftPressureFile.close();
    leftEnergyFile.close();

    rightDensityFile.close();
    rightVelocityFile.close();
    rightPressureFile.close();
    rightEnergyFile.close();

    densityFile.close();
    velocityFile.close();
    pressureFile.close();
    energyFile.close();

    levelSetFile.close();
}
