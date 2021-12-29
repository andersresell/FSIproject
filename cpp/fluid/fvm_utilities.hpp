//
// Created by anders on 12/25/21.
//

#ifndef FSIPROJECT_FVM_UTILITIES_HPP
#define FSIPROJECT_FVM_UTILITIES_HPP

#include <vector>
#include <iostream>
#include <cmath>

#define IX(i,j) (i*(nj+4) + j) //for cell centers including ghost points

#define IXH(i,j) (i*nj + j) //horizontal access of face variables
#define IXV(i,j) (i*(nj+1) + j) //vertical access of face variables

#define IXR(i,j) (i*nj + j) //for cell centers, omitting ghost points

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct vec4 {
        //data structure to hold field variables using an AoS structure
    double u1, u2, u3, u4;

    vec4 operator+(const vec4 &rhs) const;
    vec4 operator-(const vec4 &rhs) const;
    //vec4 operator*(const double &rhs) const;
    friend inline vec4 operator*(const double &lhs, const vec4 &rhs);
    //friend inline vec4 operator-(const vec4& lhs, const vec4 &rhs);
    //friend inline void operator-(vec4& lhs, const vec4& rhs);
    //void operator-(const vec4 &rhs);
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
