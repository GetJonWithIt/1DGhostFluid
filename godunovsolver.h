#ifndef GODUNOVSOLVER_H
#define GODUNOVSOLVER_H

#include "hllcsolver.h"
#include "solvers.h"
#include "ghostfluidmethod.h"
using namespace std;

class GodunovSolver
{
public:
    GodunovSolver();

    static void computeGodunovTimeStep(vector<StateVector> & newCells, vector<StateVector> & currentCells, double cellSpacing, double timeStep);

    static vector<StateVector> solve(vector<StateVector> & cells, double cellSpacing, double CFLCoefficient, double finalTime);
    static MultimaterialSystem solveWithBasicGhostFluid(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime);
    static MultimaterialSystem solveWithRiemannGhostFluid(MultimaterialSystem multimaterialSystem, double cellSpacing, double CFLCoefficient, double finalTime);
};

#endif // GODUNOVSOLVER_H
