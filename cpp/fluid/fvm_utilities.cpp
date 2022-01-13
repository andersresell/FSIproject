//
// Created by anders on 12/25/21.
//
#include "fvm_utilities.hpp"

namespace fluid {

    ExternalBCs::ExternalBCs(int ni, int nj, ExternalBC &west, ExternalBC &east, ExternalBC &south, ExternalBC &north)
            : ni{ni}, nj{nj}, west{west}, east{east}, south{south}, north{north} {
    }

    void ExternalBCs::set_BCs(vec4 *U_in) {
        for (int j{0}; j < nj + 4; j++) {
            U_in[IX(0, j)] = west.select_bc_type(U_in[IX(3, j)]);
            U_in[IX(1, j)] = west.select_bc_type(U_in[IX(2, j)]);

            U_in[IX(ni + 3, j)] = east.select_bc_type(U_in[IX(ni, j)]);
            U_in[IX(ni + 2, j)] = east.select_bc_type(U_in[IX(ni + 1, j)]);
        }
        for (int i{0}; i < ni + 4; i++) {
            U_in[IX(i, 0)] = south.select_bc_type(U_in[IX(i, 3)]);
            U_in[IX(i, 1)] = south.select_bc_type(U_in[IX(i, 2)]);

            U_in[IX(i, nj + 3)] = north.select_bc_type(U_in[IX(i, nj)]);
            U_in[IX(i, nj + 2)] = north.select_bc_type(U_in[IX(i, nj + 1)]);
        }
        /*
        //corner ghost points are set to zero so that they don't affect the cfl condition. This is neccessary since they are
        //never updated
        U_in[IX(0,0)].set_zero();
        U_in[IX(1,0)].set_zero();
        U_in[IX(0,1)].set_zero();
        U_in[IX(1,1)].set_zero();

       U_in[IX(ni+2,0)].set_zero();
        U_in[IX(ni+3,0)].set_zero();
        U_in[IX(ni+2,1)].set_zero();
        U_in[IX(ni+3,1)].set_zero();

        U_in[IX(ni+2,nj+2)].set_zero();
        U_in[IX(ni+3,nj+2)].set_zero();
        U_in[IX(ni+2,nj+3)].set_zero();
        U_in[IX(ni+3,nj+3)].set_zero();

        U_in[IX(0,nj+2)].set_zero();
        U_in[IX(1,nj+2)].set_zero();
        U_in[IX(0,nj+3)].set_zero();
        U_in[IX(1,nj+3)].set_zero();*/
    }

    void field2console(vec4 *U_in, int ni, int nj, int fieldvar) {
        //for debugging purposes. prints the field of field variable i in a carthesian manner
        assert(fieldvar >= 0 && fieldvar <= 3);

        for (int i{0}; i < ni; i++) {
            for (int j{0}; j < nj; j++) {
                auto Uptr = (double *) &U_in[IX(i, j)];
                std::cout << Uptr[fieldvar] << " ";
            }
            std::cout << std::endl;
        }

    }


}