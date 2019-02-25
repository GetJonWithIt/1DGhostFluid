#include "multimaterialsystem.h"

// This class encapsulates a multimaterial system, consisting of two distinct material domains, and a cell-centred level set function to separate them.
MultimaterialSystem::MultimaterialSystem()
{
}

// Constructs a multimaterial system with specified left and right material cells, and a specified cell-centred level set function.
MultimaterialSystem::MultimaterialSystem(vector<StateVector> newLeftMaterialCells, vector<StateVector> newRightMaterialCells, vector<double> newLevelSetFunction)
{
    leftMaterialCells = newLeftMaterialCells;
    rightMaterialCells = newRightMaterialCells;
    levelSetFunction = newLevelSetFunction;
}

// Sets the left material cells.
void MultimaterialSystem::setLeftMaterialCells(vector<StateVector> newLeftMaterialCells)
{
    leftMaterialCells = newLeftMaterialCells;
}

// Sets the right material cells.
void MultimaterialSystem::setRightMaterialCells(vector<StateVector> newRightMaterialCells)
{
    rightMaterialCells = newRightMaterialCells;
}

// Sets the cell-centred level set function.
void MultimaterialSystem::setLevelSetFunction(vector<double> newLevelSetFunction)
{
    levelSetFunction = newLevelSetFunction;
}

// Retrieves the left material cells.
vector<StateVector> MultimaterialSystem::getLeftMaterialCells()
{
    return leftMaterialCells;
}

// Retrieves the right material cells.
vector<StateVector> MultimaterialSystem::getRightMaterialCells()
{
    return rightMaterialCells;
}

// Retrieves the cell-centered level set function.
vector<double> MultimaterialSystem::getLevelSetFunction()
{
    return levelSetFunction;
}
