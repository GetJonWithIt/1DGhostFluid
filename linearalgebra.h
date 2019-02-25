#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <vector>
using namespace std;

class LinearAlgebra
{
public:
    LinearAlgebra();

    static vector<double> addVectors(vector<double> vector1, vector<double> vector2);
    static vector<double> subtractVectors(vector<double> vector1, vector<double> vector2);

    static vector<double> multiplyVector(double scalar, vector<double> vector1);
};

#endif // LINEARALGEBRA_H
