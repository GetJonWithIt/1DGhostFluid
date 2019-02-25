#ifndef MULTIMATERIALSYSTEM_H
#define MULTIMATERIALSYSTEM_H

#include "statevector.h"

class MultimaterialSystem
{
public:
    MultimaterialSystem();
    MultimaterialSystem(vector<StateVector> newLeftMaterialCells, vector<StateVector> newRightMaterialCells, vector<double> newLevelSetFunction);

    void setLeftMaterialCells(vector<StateVector> newLeftMaterialCells);
    void setRightMaterialCells(vector<StateVector> newRightMaterialCells);
    void setLevelSetFunction(vector<double> newLevelSetFunction);

    vector<StateVector> getLeftMaterialCells();
    vector<StateVector> getRightMaterialCells();
    vector<double> getLevelSetFunction();

private:
    vector<StateVector> leftMaterialCells;
    vector<StateVector> rightMaterialCells;
    vector<double> levelSetFunction;
};

#endif // MULTIMATERIALSYSTEM_H
