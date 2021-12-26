//
// Created by anders on 12/25/21.
//
#include "structuredFVM.hpp"
FieldVec4::FieldVec4()
{
    U = new vec4[size];
}

void FieldVec4::operator*(double rhs){
    for (int i = 0; i < size; i++) {
        U[i].u0 *= rhs;
        U[i].u1 *= rhs;
        U[i].u2 *= rhs;
        U[i].u3 *= rhs;
    }
}