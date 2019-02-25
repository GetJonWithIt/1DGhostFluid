#include "linearalgebra.h"

// This class encapsulates the routine vector operations required by (stiffened gas) Euler equation solvers.
LinearAlgebra::LinearAlgebra()
{
}

// Adds two vectors (vector1 and vector2) together.
vector<double> LinearAlgebra::addVectors(vector<double> vector1, vector<double> vector2)
{
    int componentCount = vector1.size();
    vector<double> sumVector(componentCount);

    for (int i = 0; i < componentCount; i++)
    {
        sumVector[i] = vector1[i] + vector2[i];
    }

    return sumVector;
}

// Subtracts one vector (vector2) from another (vector1).
vector<double> LinearAlgebra::subtractVectors(vector<double> vector1, vector<double> vector2)
{
    return addVectors(vector1, multiplyVector(-1.0, vector2));
}

// Multiplies a vector (vector1) by a specified scalar quantity (scalar).
vector<double> LinearAlgebra::multiplyVector(double scalar, vector<double> vector1)
{
    int componentCount = vector1.size();
    vector<double> productVector(componentCount);

    for (int i = 0; i < componentCount; i++)
    {
        productVector[i] = scalar * vector1[i];
    }

    return productVector;
}
