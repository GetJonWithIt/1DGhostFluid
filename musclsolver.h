#ifndef MUSCLSOLVER_H
#define MUSCLSOLVER_H

#include "statevector.h"
#include "hllcsolver.h"
#include "solvers.h"
#include "ghostfluidmethod.h"
using namespace std;

class MUSCLSolver
{
public:
    MUSCLSolver();

    static vector<double> computeMUSCLFlux(StateVector leftLeftStateVector, StateVector leftStateVector, StateVector rightStateVector, StateVector rightRightStateVector, double cellSpacing,
                                           double timeStep, double bias, int limiter);
    static void computeMUSCLTimeStep(vector<StateVector> & newCells, vector<StateVector> & currentCells, double cellSpacing, double timeStep, double bias, int limiter);

    static vector<StateVector> solve(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int limiter);

    static MultimaterialSystem solveWithBasicGhostFluid(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int limiter);
    static MultimaterialSystem solveWithBasicGhostFluidReinitialisation(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                                     int limiter, int reinitialisationIteration);

    static MultimaterialSystem solveWithRiemannGhostFluid(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias, int limiter);
    static MultimaterialSystem solveWithRiemannGhostFluidReinitialisation(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime, double bias,
                                                                          int limiter, int reinitialisationIteration);
};

#endif // MUSCLSOLVER_H
