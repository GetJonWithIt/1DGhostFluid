#include "wavespeeds.h"

// This class encapsulates the triple of wave speed estimates required by the HLLC solver, as detailed in Toro.

// Constructs a triple of wave speed estimates, with all equal to zero.
WaveSpeeds::WaveSpeeds()
{
    leftWaveSpeed = 0.0;
    rightWaveSpeed = 0.0;
    starRegionWaveSpeed = 0.0;
}

// Constructs a triple of wave speed estimates, with specified values for the left, right and star region wave speeds.
WaveSpeeds::WaveSpeeds(double newLeftWaveSpeed, double newRightWaveSpeed, double newStarRegionWaveSpeed)
{
    leftWaveSpeed = newLeftWaveSpeed;
    rightWaveSpeed = newRightWaveSpeed;
    starRegionWaveSpeed = newStarRegionWaveSpeed;
}

// Sets the left wave speed estimate.
void WaveSpeeds::setLeftWaveSpeed(double newLeftWaveSpeed)
{
    leftWaveSpeed = newLeftWaveSpeed;
}

// Sets the right wave speed estimate.
void WaveSpeeds::setRightWaveSpeed(double newRightWaveSpeed)
{
    rightWaveSpeed = newRightWaveSpeed;
}

// Sets the star region wave speed estimate.
void WaveSpeeds::setStarRegionWaveSpeed(double newStarRegionWaveSpeed)
{
    starRegionWaveSpeed = newStarRegionWaveSpeed;
}

// Retrieves the left wave speed estimate.
double WaveSpeeds::getLeftWaveSpeed()
{
    return leftWaveSpeed;
}

// Retrieves the right wave speed estimate.
double WaveSpeeds::getRightWaveSpeed()
{
    return rightWaveSpeed;
}

// Retrieves the star region wave speed estimate.
double WaveSpeeds::getStarRegionWaveSpeed()
{
    return starRegionWaveSpeed;
}
