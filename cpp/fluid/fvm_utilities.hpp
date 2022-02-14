//
// Created by anders on 12/25/21.
//

#ifndef FSIPROJECT_FVM_UTILITIES_HPP
#define FSIPROJECT_FVM_UTILITIES_HPP

#include "../../includes.hpp"

//Macros for mapping values from 2D grid to 1D arrays
#define IX(i,j) ((i)*(nj+4) + (j)) //for cell centers including ghost points
#define IXH(i,j) ((i)*nj + (j) - nj - 2) //horizontal access of extrapolated variables.
#define IXV(i,j) ((i)*(nj+2) + (j) - 2*nj - 5) //Vertical access of extrapolated variables.
#define IXF(i,j) ((i)*nj + (j) - nj - 2) //Access of horizontal fluxes. The fluxes are defined at the right cell face
#define IXG(i,j) ((i)*(nj+1) + (j) -2*nj -3) //Access of vertical fluxes. The fluxes are defined at the north cell face
#define IXR(i,j) ((i)*nj + (j) - 2*(nj+1)) //for fields only defined in the solution domain (omitting ghost points)

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

namespace fluid {

    constexpr double Gamma = 1.4;

    struct vec4 {
        //data structure to hold field variables using an AoS structure
        double u1, u2, u3, u4;

        vec4 operator+(const vec4 &rhs) const { return {u1 + rhs.u1, u2 + rhs.u2, u3 + rhs.u3, u4 + rhs.u4}; }

        vec4 operator-(const vec4 &rhs) const { return {u1 - rhs.u1, u2 - rhs.u2, u3 - rhs.u3, u4 - rhs.u4}; }

        void operator*=(double rhs) {
            u1 *= rhs;
            u2 *= rhs;
            u3 *= rhs;
            u4 *= rhs;
        }

        friend vec4 operator*(const double &lhs, const vec4 &rhs) {
            return {lhs * rhs.u1, lhs * rhs.u2, lhs * rhs.u3, lhs * rhs.u4};
        }

        friend std::ostream &operator<<(std::ostream &os, const vec4 &rhs) {
            return os << rhs.u1 << "," << rhs.u2 << "," << rhs.u3 << "," << rhs.u4;
        }

    };

    enum class OdeScheme {
        ExplicitEuler, TVD_RK3
    };
    enum class FluxScheme {
        Rusanov, HLLC
    };



    enum class BC_Type{InvicidWall, SupersonicInflow, NonreflectingOutflow};

    class ExternalBCs {
        //ExternalBC &west, &east, &south, &north;
        BC_Type west, east, south, north;
        const int ni, nj;
        vec4 U_inf; //Used in case of supersonic inflow
    public:
        //ExternalBCs(int ni, int nj, ExternalBC &west, ExternalBC &east, ExternalBC &south, ExternalBC &north);
        ExternalBCs(int ni, int nj, BC_Type west, BC_Type east, BC_Type south, BC_Type north,
                    double M_inf = 0, double p_inf = 1e5, double rho_inf = 1.2);
        void set_BCs(vec4 *U_in);
        static vec4 set_vertical_invicid_wall(const vec4& U_in);
        static vec4 set_horizontal_invicid_wall(const vec4& U_in);
    };




}

#endif //FSIPROJECT_FVM_UTILITIES_HPP
