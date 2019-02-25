#include "slopelimiter.h"

// This class encapsulates the methods for slope limiting and boundary value extrapolation, as required by the MUSCL-Hancock scheme for the (stiffened gas) Euler equations)
SlopeLimiter::SlopeLimiter()
{
}

// Computes the value of the chosen slope limiter (0 = Super-Bee, 1 = Van-Leer, 2 = Min-Bee).
double SlopeLimiter::computeSlopeLimiter(double steepness, double bias, int limiter)
{
    if (limiter == 0)
    {
        return computeSuperBeeLimiter(steepness, bias);
    }
    else if (limiter == 1)
    {
        return computeVanLeerLimiter(steepness, bias);
    }
    else
    {
        return computeMinBeeLimiter(steepness, bias);
    }
}

// Computes the value of the xi(R) parameter, as detailed in Toro.
double SlopeLimiter::computeRCoefficient(double steepness, double bias)
{
    return 2.0 / ((1.0 - bias) + ((1.0 + bias) * steepness));
}

// Computes the value of the Super-Bee slope limiter.
double SlopeLimiter::computeSuperBeeLimiter(double steepness, double bias)
{
    if (steepness < 0)
    {
        return 0.0;
    }
    else if (steepness < 0.5)
    {
        return 2.0 * steepness;
    }
    else if (steepness < 1)
    {
        return 1.0;
    }
    else
    {
        return min(min(steepness, computeRCoefficient(steepness, bias)), 2.0);
    }
}

// Computes the value of the Van-Leer slope limiter.
double SlopeLimiter::computeVanLeerLimiter(double steepness, double bias)
{
    if (steepness < 0)
    {
        return 0.0;
    }
    else
    {
        return min((2.0 * steepness) / (1.0 + steepness), computeRCoefficient(steepness, bias));
    }
}

// Computes the value of the Min-Bee slope limiter.
double SlopeLimiter::computeMinBeeLimiter(double steepness, double bias)
{
    if (steepness < 0)
    {
        return 0.0;
    }
    else if (steepness < 1)
    {
        return steepness;
    }
    else
    {
        return min(1.0, computeRCoefficient(steepness, bias));
    }
}

// Computes the limited slope vector between three neighbouring cells, with states given by leftStateVector, middleStateVector and rightStateVector, respectively.
vector<double> SlopeLimiter::computeSlopeVector(StateVector leftStateVector, StateVector middleStateVector, StateVector rightStateVector, double bias, int limiter)
{
    vector<double> leftConservedVariableDifference = LinearAlgebra::subtractVectors(middleStateVector.computeConservedVariableVector(), leftStateVector.computeConservedVariableVector());
    vector<double> rightConservedVariableDifference = LinearAlgebra::subtractVectors(rightStateVector.computeConservedVariableVector(), middleStateVector.computeConservedVariableVector());

    vector<double> slopeVector = LinearAlgebra::addVectors(LinearAlgebra::multiplyVector(0.5 * (1.0 + bias), leftConservedVariableDifference),
                                                           LinearAlgebra::multiplyVector(0.5 * (1.0 - bias), rightConservedVariableDifference));
    int componentCount = slopeVector.size();

    for (int i = 0; i < componentCount; i++)
    {
        double numerator = leftConservedVariableDifference[i];
        double denominator = rightConservedVariableDifference[i];

        if (abs(numerator) < pow(10, -8))
        {
            numerator = pow(10, -8);
        }
        if (abs(denominator) < pow(10, -8))
        {
            denominator = pow(10, -8);
        }

        double steepness = numerator / denominator;
        slopeVector[i] *= computeSlopeLimiter(steepness, bias, limiter);
    }

    return slopeVector;
}
