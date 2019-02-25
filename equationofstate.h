#ifndef EQUATIONOFSTATE_H
#define EQUATIONOFSTATE_H

#include <cmath>
using namespace std;

class EquationOfState
{
public:
    EquationOfState();

    static double computeInternalEnergy(double density, double pressure, double adiabaticIndex, double stiffeningParameter);
    static double computePressure(double density, double internalEnergy, double adiabaticIndex, double stiffeningParameter);
    static double computeSoundSpeed(double density, double pressure, double adiabaticIndex, double stiffeningParameter);
    static double computeEntropy(double density, double pressure, double adiabaticIndex);
};

#endif // EQUATIONOFSTATE_H
