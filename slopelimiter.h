#ifndef SLOPELIMITER_H
#define SLOPELIMITER_H

#include "statevector.h"
#include "linearalgebra.h"
using namespace std;

class SlopeLimiter
{
public:
    SlopeLimiter();

    static double computeSlopeLimiter(double steepness, double bias, int limiter);
    static double computeRCoefficient(double sttepeness, double bias);

    static double computeSuperBeeLimiter(double steepness, double bias);
    static double computeVanLeerLimiter(double steepness, double bias);
    static double computeMinBeeLimiter(double steepness, double bias);

    static vector<double> computeSlopeVector(StateVector leftStateVector, StateVector middleStateVector, StateVector rightStateVector, double bias, int limiter);
};

#endif // SLOPELIMITER_H
