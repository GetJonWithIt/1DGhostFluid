#ifndef SOLVERS_H
#define SOLVERS_H

#include "statevector.h"
#include "slopelimiter.h"
#include "multimaterialsystem.h"
using namespace std;

class Solvers
{
public:
    Solvers();

    static vector<StateVector> insertBoundaryCells(vector<StateVector> & cells, int boundaryCells);

    static double computeMaximumWaveSpeed(vector<StateVector> & cells);
    static double computeStableTimeStep(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime, int currentIteration);

    static double computeGhostFluidMaximumWaveSpeed(MultimaterialSystem multimaterialSystem);
    static double computeGhostFluidStableTimeStep(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double currentTime, double finalTime,
                                                  int currentIteration);

    static StateVector computeEvolvedState(StateVector leftStateVector, StateVector middleStateVector, StateVector rightStateVector, double cellSpacing, double timeStep, double bias,
                                           int limiter, int side);
};

#endif // SOLVERS_H
