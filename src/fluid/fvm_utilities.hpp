//
// Created by anders on 12/25/21.
//
#pragma once

#include "../includes.hpp"

// Macros for mapping values from 2D grid to 1D arrays
#define IX(i, j) ((i) * (nj + 4) + (j))             // for cell centers including ghost points
#define IXH(i, j) ((i)*nj + (j)-nj - 2)             // horizontal access of extrapolated variables.
#define IXV(i, j) ((i) * (nj + 2) + (j)-2 * nj - 5) // Vertical access of extrapolated variables.
#define IXF(i, j) ((i)*nj + (j)-nj - 2) // Access of horizontal fluxes. The fluxes are defined at the right cell face
#define IXG(i, j)                                                                                                      \
    ((i) * (nj + 1) + (j)-2 * nj - 3) // Access of vertical fluxes. The fluxes are defined at the north cell face
#define IXR(i, j) ((i)*nj + (j)-2 * (nj + 1)) // for fields only defined in the solution domain (omitting ghost points)

template <typename T> inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> inline T sqr(T val) {
    return val * val;
}

namespace fluid {

constexpr double Gamma = 1.4;

struct vec4 {
    // data structure to hold field variables using an AoS structure
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

    bool isnan() { // for debugging purposes
        return (std::isnan(u1) || std::isnan(u2) || std::isnan(u3) || std::isnan(u4));
    }
};

inline double calc_P(const vec4 &U_in) {
    return (Gamma - 1) * (U_in.u4 - 0.5 * (U_in.u2 * U_in.u2 + U_in.u3 * U_in.u3) / U_in.u1);
}

inline vec4 conserved2primitive(const vec4 &U_in) {
    return {U_in.u1, U_in.u2 / U_in.u1, U_in.u3 / U_in.u1, calc_P(U_in)};
}

inline vec4 primitive2conserved(const vec4 &V_in) {
    return {V_in.u1, V_in.u2 * V_in.u1, V_in.u3 * V_in.u1,
            V_in.u4 / (Gamma - 1) + 0.5 * V_in.u1 * (V_in.u2 * V_in.u2 + V_in.u3 * V_in.u3)};
}

inline vec4 calc_F(const vec4 &U_in) {
    double P = calc_P(U_in);
    return {U_in.u2, U_in.u2 * U_in.u2 / U_in.u1 + P, U_in.u2 * U_in.u3 / U_in.u1, (U_in.u4 + P) * U_in.u2 / U_in.u1};
}

inline vec4 calc_G(const vec4 &U_in) {
    double P = calc_P(U_in);
    return {U_in.u3, U_in.u2 * U_in.u3 / U_in.u1, U_in.u3 * U_in.u3 / U_in.u1 + P, (U_in.u4 + P) * U_in.u3 / U_in.u1};
}

inline double calc_sound_speed(const vec4 &U_in) {
    return sqrt(Gamma / U_in.u1 * calc_P(U_in));
}

inline vec4 minmod(const vec4 &a, const vec4 &b) {
    return {sgn(a.u1) * std::max(0.0, std::min(std::abs(a.u1), sgn(a.u1) * b.u1)),
            sgn(a.u2) * std::max(0.0, std::min(std::abs(a.u2), sgn(a.u2) * b.u2)),
            sgn(a.u3) * std::max(0.0, std::min(std::abs(a.u3), sgn(a.u3) * b.u3)),
            sgn(a.u4) * std::max(0.0, std::min(std::abs(a.u4), sgn(a.u4) * b.u4))};
}

inline double calc_sprad_x(const vec4 &U_in) {
    return std::abs(U_in.u2 / U_in.u1) + calc_sound_speed(U_in);
}

inline double calc_sprad_y(const vec4 &U_in) {
    return std::abs(U_in.u3 / U_in.u1) + calc_sound_speed(U_in);
}

enum class OdeScheme {
    ExplicitEuler,
    TVD_RK3
};
enum class FluxScheme {
    Rusanov,
    HLLC
};

enum class Limiter {
    Minmod,
    MC
};

enum class BC_Type {
    InvicidWall,
    SupersonicInflow,
    NonreflectingOutflow,
    TimeHistory,
    BlastLoad
};

struct TimeHistory {
    double t;
    vec4 U_inner_GP;
    vec4 U_outer_GP;
};

} // namespace fluid
