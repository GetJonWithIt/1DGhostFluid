#include "exactsolver.h"

// This class encapsulates the general exact Riemann solver for the Euler equations, as detailed in Section 4.3 of Toro.
ExactSolver::ExactSolver()
{
}

// Computes the "A" coefficient from Toro's exact solver.
double ExactSolver::computeACoefficient(double density, double adiabaticIndex)
{
    return 2.0 / (density * (adiabaticIndex + 1.0));
}

// Computes the "B" coefficient from Toro's exact solver.
double ExactSolver::computeBCoefficient(double pressure, double adiabaticIndex, double stiffeningParameter)
{
    return pressure * ((adiabaticIndex - 1.0) / (adiabaticIndex + 1.0)) + ((2.0 * adiabaticIndex * stiffeningParameter) / (adiabaticIndex + 1.0));
}

// Computes the speed of sound inside the star region.
double ExactSolver::computeStarRegionSoundSpeed(double intermediatePressure, StateVector stateVector)
{
    double adiabaticIndex = stateVector.getAdiabaticIndex();
    double pressure = stateVector.getPressure();
    double soundSpeed = stateVector.computeSoundSpeed();
    double stiffeningParameter = stateVector.getStiffeningParameter();

    return soundSpeed * pow((intermediatePressure + stiffeningParameter) / (pressure + stiffeningParameter), (adiabaticIndex - 1.0) / (2.0 * adiabaticIndex));
}

// Computes the vector of primitive variables (i.e. density, x/y/z-velocities, pressure, adiabatic exponent and stiffening parameter) inside a rarefaction fan.
vector<double> ExactSolver::computeRarefactionPrimitiveVariableVector(double waveSpeed, StateVector stateVector, double soundSpeed)
{
    double density = stateVector.getDensity();
    double xVelocity = stateVector.getXVelocity();
    double pressure = stateVector.getPressure();
    double adiabaticIndex = stateVector.getAdiabaticIndex();
    double stiffeningParameter = stateVector.getStiffeningParameter();

    double waveSpeedCoefficient = (2.0 / (adiabaticIndex + 1.0)) + (adiabaticIndex - 1.0) * ((xVelocity - waveSpeed) / ((adiabaticIndex + 1.0) * soundSpeed));
    vector<double> rarefactionPrimitiveVariableVector(7);

    rarefactionPrimitiveVariableVector[0] = density * pow(waveSpeedCoefficient, 2.0 / (adiabaticIndex - 1.0));
    rarefactionPrimitiveVariableVector[1] = 2.0 * (soundSpeed + (adiabaticIndex - 1.0) * (xVelocity / 2.0) + waveSpeed) / (adiabaticIndex + 1.0);
    rarefactionPrimitiveVariableVector[2] = 0.0;
    rarefactionPrimitiveVariableVector[3] = 0.0;
    rarefactionPrimitiveVariableVector[4] = ((pressure + stiffeningParameter) * pow(waveSpeedCoefficient, (2.0 * adiabaticIndex) / (adiabaticIndex - 1.0))) - stiffeningParameter;
    rarefactionPrimitiveVariableVector[5] = adiabaticIndex;
    rarefactionPrimitiveVariableVector[6] = stiffeningParameter;

    return rarefactionPrimitiveVariableVector;
}

// Computes the density inside the star region for a shock wave.
double ExactSolver::computeStarRegionShockDensity(double intermediatePressure, StateVector stateVector)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();
    double adiabaticIndex = stateVector.getAdiabaticIndex();
    double stiffeningParameter = stateVector.getStiffeningParameter();

    double numerator = ((intermediatePressure + stiffeningParameter) / (pressure + stiffeningParameter)) + ((adiabaticIndex - 1.0) / (adiabaticIndex + 1.0));
    double denominator = ((adiabaticIndex - 1.0) / (adiabaticIndex + 1.0)) * ((intermediatePressure + stiffeningParameter) / (pressure + stiffeningParameter)) + 1.0;

    return density * (numerator / denominator);
}

// Computes the density inside the star region for a rarefaction fan.
double ExactSolver::computeStarRegionRarefactionDensity(double intermediatePressure, StateVector stateVector)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();
    double adiabaticIndex = stateVector.getAdiabaticIndex();
    double stiffeningParameter = stateVector.getStiffeningParameter();

    return density * pow((intermediatePressure + stiffeningParameter) / (pressure + stiffeningParameter), 1.0 / adiabaticIndex);
}

// Computes the x-velocity inside the star region.
double ExactSolver::computeStarRegionVelocity(double intermediatePressure, StateVector leftStateVector, StateVector rightStateVector)
{
    return 0.5 * (leftStateVector.getXVelocity() + rightStateVector.getXVelocity() + computeWaveJumpFunctionComponent(intermediatePressure, rightStateVector) -
                  computeWaveJumpFunctionComponent(intermediatePressure, leftStateVector));
}

// Computes the pressure inside the star region.
double ExactSolver::computeStarRegionPressure(StateVector leftStateVector, StateVector rightStateVector)
{
    double initialPressure = 0.5 * (leftStateVector.getPressure() + rightStateVector.getPressure());
    double newPressure = 3.0 * initialPressure;
    int iterationCount = 0;

    while (computePressureChange(newPressure, initialPressure) >= pow(10, -8))
    {
        initialPressure = newPressure;
        newPressure = initialPressure - (computeWaveJumpFunction(initialPressure, leftStateVector, rightStateVector) /
                                         computeWaveJumpFunctionDerivative(initialPressure, leftStateVector, rightStateVector));

        if (newPressure < 0)
        {
            newPressure = pow(10, -8);
        }

        iterationCount += 1;
        if (iterationCount > 10000)
        {
            newPressure = 0.00189;
        }
    }

    return newPressure;
}

// Computes the function which characterises the jump in wave quantities over either the left- or right-moving wave.
double ExactSolver::computeWaveJumpFunctionComponent(double newPressure, StateVector stateVector)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();
    double soundSpeed = stateVector.computeSoundSpeed();
    double adiabaticIndex = stateVector.getAdiabaticIndex();
    double stiffeningParameter = stateVector.getStiffeningParameter();

    if (newPressure > pressure)
    {
        double coefficient = sqrt(computeACoefficient(density, adiabaticIndex) / (newPressure + computeBCoefficient(pressure, adiabaticIndex, stiffeningParameter)));
        return (newPressure - pressure) * coefficient;
    }
    else
    {
        double coefficient = pow((newPressure + stiffeningParameter) / (pressure + stiffeningParameter), (adiabaticIndex - 1.0) / (2.0 * adiabaticIndex));
        return 2.0 * (soundSpeed / (adiabaticIndex - 1.0)) * (coefficient - 1.0);
    }
}

// Computes the function which characterises the jump in wave quantities over both left- and right-moving waves.
double ExactSolver::computeWaveJumpFunction(double newPressure, StateVector leftStateVector, StateVector rightStateVector)
{
    return computeWaveJumpFunctionComponent(newPressure, leftStateVector) + computeWaveJumpFunctionComponent(newPressure, rightStateVector) + rightStateVector.getXVelocity() -
            leftStateVector.getXVelocity();
}

// Computes the first derivative of the function which characterises the jump in wave quantities over either the left- or right-moving wave.
double ExactSolver::computeWaveJumpFunctionDerivativeComponent(double newPressure, StateVector stateVector)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();
    double soundSpeed = stateVector.computeSoundSpeed();
    double adiabaticIndex = stateVector.getAdiabaticIndex();
    double stiffeningParameter = stateVector.getStiffeningParameter();

    if (newPressure > pressure)
    {
        double coefficient = sqrt(computeACoefficient(density, adiabaticIndex) / (newPressure + computeBCoefficient(pressure, adiabaticIndex, stiffeningParameter)));
        return (1.0 - ((newPressure - pressure) / (2.0 * (newPressure + computeBCoefficient(pressure, adiabaticIndex, stiffeningParameter))))) * coefficient;
    }
    else
    {
        double coefficient = pow((newPressure + stiffeningParameter) / (pressure + stiffeningParameter), -(adiabaticIndex + 1.0) / (2.0 * adiabaticIndex));
        return (coefficient * soundSpeed) / (adiabaticIndex * (pressure + stiffeningParameter));
    }
}

// Computes the first derivative of the function which characterises the jump in wave quantities over both left- and right-moving waves.
double ExactSolver::computeWaveJumpFunctionDerivative(double newPressure, StateVector leftStateVector, StateVector rightStateVector)
{
    return computeWaveJumpFunctionDerivativeComponent(newPressure, leftStateVector) + computeWaveJumpFunctionDerivativeComponent(newPressure, rightStateVector);
}

// Computes the relative change in pressure between two consecutive iterations.
double ExactSolver::computePressureChange(double pressure1, double pressure2)
{
    return 2.0 * (abs(pressure1 - pressure2) / abs(pressure1 + pressure2));
}

// Computes the propagation speed of a shock wave.
double ExactSolver::computeShockSpeed(double intermediatePressure, StateVector stateVector)
{
    double density = stateVector.getDensity();
    double pressure = stateVector.getPressure();
    double adiabaticIndex = stateVector.getAdiabaticIndex();
    double stiffeningParameter = stateVector.getStiffeningParameter();

    return sqrt((intermediatePressure + computeBCoefficient(pressure, adiabaticIndex, stiffeningParameter)) / computeACoefficient(density, adiabaticIndex));
}

// Computes the density inside the star region.
double ExactSolver::computeStarRegionDensity(double intermediatePressure, StateVector stateVector)
{
    double pressure = stateVector.getPressure();

    if (intermediatePressure < pressure)
    {
        return computeStarRegionRarefactionDensity(intermediatePressure, stateVector);
    }
    else
    {
        return computeStarRegionShockDensity(intermediatePressure, stateVector);
    }
}

// Solves the Riemann problem between leftStateVector and rightStateVector, and returns the exact solution for the specified position, time and interface location.
StateVector ExactSolver::solve(double position, double time, double interfaceLocation, StateVector leftStateVector, StateVector rightStateVector)
{
    double waveSpeed = (position - interfaceLocation) / time;

    double leftDensity = leftStateVector.getDensity();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftPressure = leftStateVector.getPressure();
    double leftSoundSpeed = leftStateVector.computeSoundSpeed();

    double rightDensity = rightStateVector.getDensity();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightPressure = rightStateVector.getPressure();
    double rightSoundSpeed = rightStateVector.computeSoundSpeed();

    double intermediatePressure = computeStarRegionPressure(leftStateVector, rightStateVector);
    double intermediateVelocity = computeStarRegionVelocity(intermediatePressure, leftStateVector, rightStateVector);

    if (waveSpeed < intermediateVelocity)
    {
        if (intermediatePressure < leftPressure)
        {
            if (waveSpeed < (leftXVelocity - leftSoundSpeed))
            {
                return leftStateVector;
            }
            else
            {
                double leftAcousticWaveSpeed = intermediateVelocity - computeStarRegionSoundSpeed(intermediatePressure, leftStateVector);

                if (waveSpeed < leftAcousticWaveSpeed)
                {
                    return StateVector(computeRarefactionPrimitiveVariableVector(waveSpeed, leftStateVector, leftSoundSpeed));
                }
                else
                {
                    return StateVector(computeStarRegionRarefactionDensity(intermediatePressure, leftStateVector), intermediateVelocity, 0.0, 0.0, intermediatePressure,
                                       leftStateVector.getAdiabaticIndex(), leftStateVector.getStiffeningParameter());
                }
            }
        }
        else
        {
            double leftShockSpeed = leftXVelocity - (computeShockSpeed(intermediatePressure, leftStateVector) / leftDensity);

            if (waveSpeed < leftShockSpeed)
            {
                return leftStateVector;
            }
            else
            {
                return StateVector(computeStarRegionShockDensity(intermediatePressure, leftStateVector), intermediateVelocity, 0.0, 0.0, intermediatePressure,
                                   leftStateVector.getAdiabaticIndex(), leftStateVector.getStiffeningParameter());
            }
        }
    }
    else
    {
        if (intermediatePressure < rightPressure)
        {
            if (waveSpeed > (rightXVelocity + rightSoundSpeed))
            {
                return rightStateVector;
            }
            else
            {
                double rightAcousticWaveSpeed = intermediateVelocity + computeStarRegionSoundSpeed(intermediatePressure, rightStateVector);

                if (waveSpeed > rightAcousticWaveSpeed)
                {
                    return StateVector(computeRarefactionPrimitiveVariableVector(waveSpeed, rightStateVector, -rightSoundSpeed));
                }
                else
                {
                    return StateVector(computeStarRegionRarefactionDensity(intermediatePressure, rightStateVector), intermediateVelocity, 0.0, 0.0, intermediatePressure,
                                       rightStateVector.getAdiabaticIndex(), rightStateVector.getStiffeningParameter());
                }
            }
        }
        else
        {
            double rightShockSpeed = rightXVelocity + (computeShockSpeed(intermediatePressure, rightStateVector) / rightDensity);

            if (waveSpeed > rightShockSpeed)
            {
                return rightStateVector;
            }
            else
            {
                return StateVector(computeStarRegionShockDensity(intermediatePressure, rightStateVector), intermediateVelocity, 0.0, 0.0, intermediatePressure,
                                   rightStateVector.getAdiabaticIndex(), rightStateVector.getStiffeningParameter());
            }
        }
    }
}
