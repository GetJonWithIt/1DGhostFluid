#include "hllcsolver.h"

// This class encapsulates the first-order Harten-Lax-van Leer with Contact (HLLC) Godunov-type solver for the (stiffened gas) Euler equations, as detailed in Toro.
HLLCSolver::HLLCSolver()
{
}

// Computes the density inside the star region, using the known density, the wave speed, the particle velocity, and the intermediate wave speed.
double HLLCSolver::computeStarRegionDensity(double density, double waveSpeed, double velocity, double intermediateWaveSpeed)
{
    return density * ((waveSpeed - velocity) / (waveSpeed - intermediateWaveSpeed));
}

// Computes the flux vector inside the star region, using the known state, the wave speed, and the intermediate wave speed.
vector<double> HLLCSolver::computeStarRegionFlux(StateVector stateVector, double waveSpeed, double intermediateWaveSpeed)
{
    double density = stateVector.getDensity();
    double xVelocity = stateVector.getXVelocity();
    double pressure = stateVector.getPressure();
    double totalEnergy = stateVector.computeTotalEnergy();

    vector<double> starRegionConservedVariableVector(7);
    double starRegionDensity = computeStarRegionDensity(density, waveSpeed, xVelocity, intermediateWaveSpeed);
    double energyConstant = (intermediateWaveSpeed - xVelocity) * (intermediateWaveSpeed + (pressure / (density * (waveSpeed - xVelocity))));

    starRegionConservedVariableVector[0] = starRegionDensity;
    starRegionConservedVariableVector[1] = starRegionDensity * intermediateWaveSpeed;
    starRegionConservedVariableVector[2] = starRegionDensity * stateVector.getYVelocity();
    starRegionConservedVariableVector[3] = starRegionDensity * stateVector.getZVelocity();
    starRegionConservedVariableVector[4] = starRegionDensity * ((totalEnergy / density) + energyConstant);
    starRegionConservedVariableVector[5] = starRegionDensity * stateVector.getAdiabaticIndex();
    starRegionConservedVariableVector[6] = starRegionDensity * stateVector.getStiffeningParameter();

    vector<double> conservedVariableVector = stateVector.computeConservedVariableVector();
    vector<double> fluxVector = stateVector.computeFluxVector();

    return LinearAlgebra::addVectors(fluxVector, LinearAlgebra::multiplyVector(waveSpeed, LinearAlgebra::subtractVectors(starRegionConservedVariableVector, conservedVariableVector)));
}

// Computes the HLLC flux at the interface between two neighbouring cells, with states given by leftStateVector and rightStateVector respectively.
vector<double> HLLCSolver::computeHLLCFlux(StateVector leftStateVector, StateVector rightStateVector)
{
    WaveSpeeds waveSpeeds = computeWaveSpeeds(leftStateVector, rightStateVector);

    if (0 <= waveSpeeds.getLeftWaveSpeed())
    {
        return leftStateVector.computeFluxVector();
    }
    else if (waveSpeeds.getRightWaveSpeed() <= 0)
    {
        return rightStateVector.computeFluxVector();
    }
    else if (0 <= waveSpeeds.getStarRegionWaveSpeed())
    {
        return computeStarRegionFlux(leftStateVector, waveSpeeds.getLeftWaveSpeed(), waveSpeeds.getStarRegionWaveSpeed());
    }
    else
    {
        return computeStarRegionFlux(rightStateVector, waveSpeeds.getRightWaveSpeed(), waveSpeeds.getStarRegionWaveSpeed());
    }
}

// Computes the weightings used within the wave speed estimates, given the pressure inside the star region, and the known state.
double HLLCSolver::computeWaveSpeedWeighting(double starRegionPressure, StateVector stateVector)
{
    double pressure = stateVector.getPressure();
    double adiabaticIndex = stateVector.getAdiabaticIndex();
    double stiffeningParameter = stateVector.getStiffeningParameter();

    if (starRegionPressure <= pressure)
    {
        return 1.0;
    }
    else
    {
        return sqrt(1.0 + ((adiabaticIndex + 1.0) / (2.0 * adiabaticIndex)) * ((starRegionPressure + stiffeningParameter) /  (pressure + stiffeningParameter) - 1.0));
    }
}

// Computes the HLLC wave speed estimates, given the two known states (leftStateVector and rightStateVector).
WaveSpeeds HLLCSolver::computeWaveSpeeds(StateVector leftStateVector, StateVector rightStateVector)
{
    WaveSpeeds waveSpeeds;

    double leftDensity = leftStateVector.getDensity();
    double leftSoundSpeed = leftStateVector.computeSoundSpeed();
    double leftPressure = leftStateVector.getPressure();
    double leftXVelocity = leftStateVector.getXVelocity();
    double leftAdiabaticIndex = leftStateVector.getAdiabaticIndex();
    double leftStiffeningParameter = leftStateVector.getStiffeningParameter();

    double rightDensity = rightStateVector.getDensity();
    double rightSoundSpeed = rightStateVector.computeSoundSpeed();
    double rightPressure = rightStateVector.getPressure();
    double rightXVelocity = rightStateVector.getXVelocity();
    double rightAdiabaticIndex = rightStateVector.getAdiabaticIndex();
    double rightStiffeningParameter = rightStateVector.getStiffeningParameter();

    double pressureCoefficient = 0.25 * (leftDensity + rightDensity) * (leftSoundSpeed + rightSoundSpeed);
    double PVRSPressure = (0.5 * (leftPressure + rightPressure)) + (0.5 * (leftXVelocity - rightXVelocity) * pressureCoefficient);
    PVRSPressure = max(0.0, PVRSPressure);

    double minimumPressure = min(leftPressure, rightPressure);
    double maximumPressure = max(leftPressure, rightPressure);
    double maximumPressureRatio = maximumPressure / minimumPressure;

    double intermediatePressure;
    double intermediateVelocity;

    if (((maximumPressureRatio <= 2) && (minimumPressure <= PVRSPressure) && (PVRSPressure <= maximumPressure)) ||
            abs(leftAdiabaticIndex - rightAdiabaticIndex) > pow(10, -8) || leftStiffeningParameter > 0 || rightStiffeningParameter > 0)
    {
        intermediatePressure = PVRSPressure;
        intermediateVelocity = (0.5 * (leftXVelocity + rightXVelocity)) + (0.5 * (leftPressure - rightPressure) / pressureCoefficient);
    }
    else
    {
        double adiabaticIndex = 0.5 * (leftAdiabaticIndex + rightAdiabaticIndex);

        if (PVRSPressure < minimumPressure)
        {
            double velocityExponent = (adiabaticIndex - 1.0) / (2.0 * adiabaticIndex);
            double pressureTerm = pow((leftPressure / rightPressure), velocityExponent);
            double velocityCoefficient = 2.0 / (adiabaticIndex - 1.0);
            intermediateVelocity = ((pressureTerm * (leftXVelocity / leftSoundSpeed)) + (rightXVelocity / rightSoundSpeed) + (velocityCoefficient * (pressureTerm - 1.0))) /
                    ((pressureTerm / leftSoundSpeed) + (1.0 / rightSoundSpeed));

            double pressureCoefficient = (adiabaticIndex - 1.0) / 2.0;
            double leftVelocityTerm = 1.0 + (pressureCoefficient * ((leftXVelocity - intermediateVelocity) / leftSoundSpeed));
            double rightVelocityTerm = 1.0 + (pressureCoefficient * ((intermediateVelocity - rightXVelocity) / rightSoundSpeed));
            double pressureExponent = 2.0 * (adiabaticIndex / (adiabaticIndex - 1.0));
            intermediatePressure = 0.5 * ((leftPressure * pow(leftVelocityTerm, pressureExponent)) + (rightPressure * pow(rightVelocityTerm, pressureExponent)));
        }
        else
        {
            double densityCoefficient = 2.0 / (adiabaticIndex + 1.0);
            double pressureCoefficient = (adiabaticIndex - 1.0) / (adiabaticIndex + 1.0);
            double leftPressureCoefficient = sqrt((densityCoefficient / leftDensity) / ((pressureCoefficient * leftPressure) + PVRSPressure));
            double rightPressureCoefficient = sqrt((densityCoefficient / rightDensity) / ((pressureCoefficient * rightPressure) + PVRSPressure));

            intermediatePressure = ((leftPressureCoefficient * leftPressure) + (rightPressureCoefficient * rightPressure) - (rightXVelocity - leftXVelocity)) /
                    (leftPressureCoefficient + rightPressureCoefficient);
            intermediateVelocity = (0.5 * (leftXVelocity + rightXVelocity)) + (0.5 * ((rightPressureCoefficient * (intermediatePressure - rightPressure)) -
                                                                                      (leftPressureCoefficient * (intermediatePressure - leftPressure))));
        }
    }

    waveSpeeds.setStarRegionWaveSpeed(intermediateVelocity);
    waveSpeeds.setLeftWaveSpeed(leftXVelocity - (computeWaveSpeedWeighting(intermediatePressure, leftStateVector) * leftSoundSpeed));
    waveSpeeds.setRightWaveSpeed(rightXVelocity + (computeWaveSpeedWeighting(intermediatePressure, rightStateVector) * rightSoundSpeed));

    return waveSpeeds;
}
