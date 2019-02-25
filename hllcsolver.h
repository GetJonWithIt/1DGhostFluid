#ifndef HLLCSOLVER_H
#define HLLCSOLVER_H

#include <vector>
#include "statevector.h"
#include "wavespeeds.h"
#include "linearalgebra.h"
using namespace std;

class HLLCSolver
{
public:
    HLLCSolver();

    static double computeStarRegionDensity(double density, double waveSpeed, double velocity, double intermediateWaveSpeed);

    static vector<double> computeStarRegionFlux(StateVector stateVector, double waveSpeed, double intermediateWaveSpeed);
    static vector<double> computeHLLCFlux(StateVector leftStateVector, StateVector rightStateVector);

    static double computeWaveSpeedWeighting(double starRegionPressure, StateVector stateVector);
    static WaveSpeeds computeWaveSpeeds(StateVector leftStateVector, StateVector rightStateVector);
};

#endif // HLLCSOLVER_H
