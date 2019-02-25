#ifndef GHOSTFLUIDMETHOD_H
#define GHOSTFLUIDMETHOD_H

#include "statevector.h"
#include "multimaterialsystem.h"
#include "exactsolver.h"
using namespace std;

class GhostFluidMethod
{
public:
    GhostFluidMethod();

    static MultimaterialSystem applyBasicGhostFluidBoundaryConditions(MultimaterialSystem multimaterialSystem);
    static MultimaterialSystem applyRiemannGhostFluidBoundaryConditions(MultimaterialSystem multimaterialSystem);
    //static MultimaterialSystem applyRiemannGhostFluidBoundaryConditions(MultimaterialSystem multimaterialSystem, double cellSpacing, double timeStep);

    static vector<double> updateLevelSetFunction(vector<double> oldLevelSetFunction, vector<StateVector> leftMaterialCells, vector<StateVector> rightMaterialCells,
                                                 double cellSpacing, double timeStep);
    static int computeInterfaceCount(vector<double> levelSetFunction);
    static vector<double> reinitialiseLevelSetFunction(vector<double> oldLevelSetFunction);
    static vector<double> reinitialiseTwoInterfaceLevelSetFunction(vector<double> oldLevelSetFunction);
};

#endif // GHOSTFLUIDMETHOD_H
