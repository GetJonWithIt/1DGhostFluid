#ifndef WAVESPEEDS_H
#define WAVESPEEDS_H


class WaveSpeeds
{   
public:
    WaveSpeeds();
    WaveSpeeds(double newLeftWaveSpeed, double newRightWaveSpeed, double newStarRegionWaveSpeed);

    void setLeftWaveSpeed(double newLeftWaveSpeed);
    void setRightWaveSpeed(double newRightWaveSpeed);
    void setStarRegionWaveSpeed(double newStarRegionWaveSpeed);

    double getLeftWaveSpeed();
    double getRightWaveSpeed();
    double getStarRegionWaveSpeed();

private:
    double leftWaveSpeed;
    double rightWaveSpeed;
    double starRegionWaveSpeed;
};

#endif // WAVESPEEDS_H
