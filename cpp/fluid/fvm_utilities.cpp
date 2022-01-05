//
// Created by anders on 12/25/21.
//
#include "fvm_utilities.hpp"

ExternalBCs::ExternalBCs(int ni, int nj, ExternalBC& west, ExternalBC& east, ExternalBC& south, ExternalBC& north)
    : ni{ni}, nj{nj}, west{west}, east{east}, south{south}, north{north}
{
}

void ExternalBCs::set_BCs(vec4* U_in) {
    for (int j{2}; j < nj + 2; j++) {
        U_in[IX(0, j)] = west.select_bc_type(U_in[IX(3, j)]);
        U_in[IX(1, j)] = west.select_bc_type(U_in[IX(2, j)]);

        U_in[IX(ni + 3, j)] = east.select_bc_type(U_in[IX(ni, j)]);
        U_in[IX(ni + 2, j)] = east.select_bc_type(U_in[IX(ni + 1, j)]);
    }
    for (int i{2}; i < ni + 2; i++) {
        U_in[IX(i, 0)] = south.select_bc_type(U_in[IX(i, 3)]);
        U_in[IX(i, 1)] = south.select_bc_type(U_in[IX(i, 2)]);

        U_in[IX(i, nj + 3)] = north.select_bc_type(U_in[IX(i, nj)]);
        U_in[IX(i, nj + 2)] = north.select_bc_type(U_in[IX(i, nj + 1)]);
    }
}
