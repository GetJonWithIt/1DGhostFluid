#ifndef TESTS_H
#define TESTS_H

#include <fstream>
#include "multimaterialsystem.h"
#include "exactsolver.h"
#include "musclsolver.h"
using namespace std;

class Tests
{
public:
    Tests();

    static void solveToroTest1Exact(int cellCount);
    static void solveToroTest2Exact(int cellCount);
    static void solveToroTest3Exact(int cellCount);
    static void solveToroTest4Exact(int cellCount);
    static void solveToroTest5Exact(int cellCount);
    static void solveChinnayyaTestExact(int cellCount);
    static void solveFedkiwTestExact(int cellCount);

    static void solveToroTest1BasicGhostFluid(int cellCount);
    static void solveToroTest1BasicGhostFluidReinitialisation(int cellCount, int reinitialisationIteration);
    static void solveToroTest2BasicGhostFluid(int cellCount);
    static void solveToroTest3BasicGhostFluid(int cellCount);

    static void solveToroTest1RiemannGhostFluid(int cellCount);
    static void solveToroTest1RiemannGhostFluidReinitialisation(int cellCount, int reinitialisationIteration);
    static void solveToroTest2RiemannGhostFluid(int cellCount);
    static void solveToroTest3RiemannGhostFluid(int cellCount);
    static void solveToroTest3RiemannGhostFluidReinitialisation(int cellCount, int reinitialisationIteration);

    static void solveFedkiwTestBasicGhostFluid(int cellCount);
    static void solveWangTestBasicGhostFluid(int cellCount);
    static void solveMach10TestBasicGhostFluid(int cellCount);
    static void solveMach10TestBasicGhostFluidReinitialisation(int cellCount, int reinitialisationIteration);

    static void solveFedkiwTestRiemannGhostFluid(int cellCount);
    static void solveWangTestRiemannGhostFluid(int cellCount);
    static void solveMach10TestRiemannGhostFluid(int cellCount);
    static void solveMach10TestRiemannGhostFluidReinitialisation(int cellCount, int reinitialisationIteration);

    static void solveChinnayyaTestBasicGhostFluid(int cellCount);
    static void solveChinnayyaTestRiemannGhostFluid(int cellCount);
    static void solveChinnayyaTestRiemannGhostFluidReinitialisation(int cellCount, int reinitialisationIteration);

    static void outputExactSolution(int cellCount, double finalTime, double interfaceLocation, StateVector leftStateVector, StateVector rightStateVector);
    static void outputSingleMaterialData(vector<StateVector> cells);
    static void outputMultimaterialData(MultimaterialSystem multimaterialSystem);
};

#endif // TESTS_H
