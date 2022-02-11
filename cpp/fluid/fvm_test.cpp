//
// Created by anders on 1/4/22.
//

#include "fvm_test.hpp"

namespace fluid {

    void set_initial_cond1(vec4 *U, int ni, int nj) {
        //Sets the pressure in the bottom left part high and the pressure in the rest low
        double rho = 1.2;
        double p_low = 1e5;
        double p_high = 10 * p_low;
        double u{0};
        double v{0};

        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                U[IX(i, j)].u1 = rho;
                U[IX(i, j)].u2 = rho * u;
                U[IX(i, j)].u3 = rho * v;
                if (i <= (ni + 3) / 2 && j <= (nj + 3) / 2) {
                    U[IX(i, j)].u4 = p_high / (fluid::Gamma - 1) + 0.5 * rho * (u * u + v * v);
                } else {
                    U[IX(i, j)].u4 = p_low / (fluid::Gamma - 1) + 0.5 * rho * (u * u + v * v);
                }
            }
        }
        std::cout << "Initial condition 1 set\n";
    }

    void set_initial_cond2(vec4* U, int ni, int nj){
        //Sets the pressure in the left part high and the pressure in the rest low
        double rho = 1.2;
        double p_low = 1e5;
        double p_high = 10 * p_low;
        double u{0};
        double v{0};

        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                U[IX(i, j)].u1 = rho;
                U[IX(i, j)].u2 = rho * u;
                U[IX(i, j)].u3 = rho * v;
                if (i <= (ni + 3) / 4 ) {
                    U[IX(i, j)].u4 = p_high / (fluid::Gamma - 1) + 0.5 * rho * (u * u + v * v);
                } else {
                    U[IX(i, j)].u4 = p_low / (fluid::Gamma - 1) + 0.5 * rho * (u * u + v * v);
                }
            }
        }
        std::cout << "Initial condition 2 set\n";
    }


    void set_initial_cond_riemann(vec4* U, int ni, int nj, const vec4& V_l, const vec4& V_r){
        if (ni%2 != 0) {
            std::cerr << "ni must be even for the riemann test\n";
            exit(1);
        }
        vec4 U_l = FVM_Solver::primitive2conserved(V_l);
        vec4 U_r = FVM_Solver::primitive2conserved(V_r);
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                if (i < ni/2+2){
                    U[IX(i,j)] = U_l;
                }
                else{
                    U[IX(i,j)] = U_r;
                }
            }
        }
        std::cout << "Initial condition Riemann set\n";
    }


}
