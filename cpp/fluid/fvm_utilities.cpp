//
// Created by anders on 12/25/21.
//
#include "fvm_utilities.hpp"


vec4 vec4::operator+(const vec4 &rhs) const { return {u1 + rhs.u1, u2 + rhs.u2, u3 + rhs.u3, u4 + rhs.u4}; }

vec4 vec4::operator-(const vec4 &rhs) const { return {u1 - rhs.u1, u2 - rhs.u2, u3 - rhs.u3, u4 - rhs.u4}; }

vec4 vec4::operator*(const double &rhs) const { return {u1 * rhs, u2 * rhs, u3 * rhs, u4 * rhs}; }

inline vec4 operator*(const double &lhs, const vec4 &rhs) {
    return {lhs * rhs.u1, lhs * rhs.u2, lhs * rhs.u3, lhs * rhs.u4};
}
inline void operator-(vec4& lhs, const vec4& rhs){
    lhs.u1 = - rhs.u1;
    lhs.u2 = - rhs.u2;
    lhs.u3 = - rhs.u3;
    lhs.u4 = - rhs.u4;
}

void vec4::operator-(const vec4 &rhs) {
    u1 = - rhs.u1;
    u2 = - rhs.u2;
    u3 = - rhs.u3;
    u4 = - rhs.u4;
}