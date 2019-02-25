#ifndef EXACTSOLVER_H
#define EXACTSOLVER_H

#include "statevector.h"

class ExactSolver
{
public:
    ExactSolver();

    static double computeACoefficient(double density, double adiabaticIndex);
    static double computeBCoefficient(double pressure, double adiabaticIndex, double stiffeningParameter);

    static double computeStarRegionSoundSpeed(double intermediatePressure, StateVector stateVector);

    static vector<double> computeRarefactionPrimitiveVariableVector(double waveSpeed, StateVector stateVector, double soundSpeed);

    static double computeStarRegionShockDensity(double intermediatePressure, StateVector stateVector);
    static double computeStarRegionRarefactionDensity(double intermediatePressure, StateVector stateVector);
    static double computeStarRegionVelocity(double intermediatePressure, StateVector leftStateVector, StateVector rightStateVector);
    static double computeStarRegionPressure(StateVector leftStateVector, StateVector rightStateVector);

    static double computeWaveJumpFunctionComponent(double newPressure, StateVector stateVector);
    static double computeWaveJumpFunction(double newPressure, StateVector leftStateVector, StateVector rightStateVector);

    static double computeWaveJumpFunctionDerivativeComponent(double newPressure, StateVector stateVector);
    static double computeWaveJumpFunctionDerivative(double newPressure, StateVector leftStateVector, StateVector rightStateVector);

    static double computePressureChange(double pressure1, double pressure2);

    static double computeShockSpeed(double intermediatePressure, StateVector stateVector);

    static double computeStarRegionDensity(double intermediatePressure, StateVector stateVector);
    static StateVector solve(double position, double time, double interfaceLocation, StateVector leftStateVector, StateVector rightStateVector);

};

#endif // EXACTSOLVER_H
