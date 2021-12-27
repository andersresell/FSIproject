//
// Created by anders on 12/25/21.
//
#include "fvm_utilities.hpp"
FieldVec4::FieldVec4()
{
    U = new vec4[size];
}

FieldVec4::FieldVec4(const FieldVec4& rhs) : FieldVec4(){ //Copy constructor
    for (int i{0}; i<size;i++){
        U[i] = rhs.U[i];
    }
}

void FieldVec4::operator*(double rhs){ //Scalar multiplication
    for (int i{0}; i < size; i++) {
        U[i].u0 *= rhs;
        U[i].u1 *= rhs;
        U[i].u2 *= rhs;
        U[i].u3 *= rhs;
    }
}