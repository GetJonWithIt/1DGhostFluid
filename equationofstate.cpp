#include "equationofstate.h"

// This class encapsulates the equation of state for a stiffened gas.
EquationOfState::EquationOfState()
{
}

// Computes the specific internal energy for a stiffened gas.
double EquationOfState::computeInternalEnergy(double density, double pressure, double adiabaticIndex, double stiffeningParameter)
{
    return (pressure + (adiabaticIndex * stiffeningParameter)) / ((adiabaticIndex - 1.0) * density);
}

// Computes the pressure for a stiffened gas.
double EquationOfState::computePressure(double density, double internalEnergy, double adiabaticIndex, double stiffeningParameter)
{
    return (internalEnergy * (adiabaticIndex - 1.0) * density) - (adiabaticIndex * stiffeningParameter);
}

// Computes the sound speed for a stiffened gas.
double EquationOfState::computeSoundSpeed(double density, double pressure, double adiabaticIndex, double stiffeningParameter)
{
    return sqrt((adiabaticIndex * (pressure + stiffeningParameter)) / density);
}

// Computes the entropy for an ideal gas, using its density, pressure and adiabatic exponent.
double EquationOfState::computeEntropy(double density, double pressure, double adiabaticIndex)
{
    return pressure / pow(density, adiabaticIndex);
}
