//
// Created by anders on 12/25/21.
//

#ifndef FSIPROJECT_FVM_UTILITIES_HPP
#define FSIPROJECT_FVM_UTILITIES_HPP

#include <vector>
#include <iostream>
#include <math.h>

#define IX(i,j) (i*(nj+4) + j) //for whole grid including ghost points

#define IXH(i,j) (i*nj + j) //horizontal access of face variables
#define IXV(i,j) (i*(nj+1) + j)

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct vec4 {
    double u1, u2, u3, u4;
    vec4 operator+(const vec4 &rhs) { return {u1 + rhs.u1, u2 + rhs.u2, u3 + rhs.u3, u4 + rhs.u4}; }
    vec4 operator*(double rhs) { return {u1 * rhs, u2 * rhs, u3 * rhs, u4 * rhs}; }
};
/*
class FieldVec4 {
    vec4 *U;
}
public:
    FieldVec4();
    FieldVec4(const FieldVec4& rhs);
    void operator*(double rhs);
    void operator*(FieldVec4& rhs);

    ~FieldVec4(){delete[] U;}
};

class FieldScalar{
    double* u;
};
*/
enum class OdeScheme{ExplicitEuler, TVD_RK3};
enum class FluxScheme{Rusanov, HLLC};
enum class BC{Wall};
struct BoundaryConditions{
    BC west, east, south, north;
};


#endif //FSIPROJECT_FVM_UTILITIES_HPP
