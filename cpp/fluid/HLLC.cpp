//
// Created by anders on 1/28/22.
//
#include "fvm_solver.hpp"
#include "fvm_utilities.hpp"

namespace fluid {
    void FVM_Solver::HLLC() {
        for (int i{1}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                if (cell_status[IX(i,j)] != CellStatus::Solid)
                    F_f[IXF(i, j)] = F_f_HLLC(U_right[IXH(i, j)], U_left[IXH(i + 1, j)]);
            }
        }
        for (int i{2}; i < ni + 2; i++) {
            for (int j{1}; j < nj + 2; j++) {
                if (cell_status[IX(i,j)] != CellStatus::Solid)
                    G_f[IXG(i, j)] = G_f_HLLC(U_up[IXV(i, j)], U_up[IXV(i, j + 1)]);
            }
        }
    }

    vec4 FVM_Solver::F_f_HLLC(const vec4 &U_L, const vec4 &U_R) {
        //The theory behind this implementation is found in Toro's book about Riemann solvers
        double rho_L = U_L.u1;
        double rho_R = U_R.u1;
        double p_L = FVM_Solver::calc_P(U_L);
        double p_R = calc_P(U_R);
        double u_L = U_L.u2 / rho_L;
        double u_R = U_R.u2 / rho_R;
        double v_L = U_L.u3 / U_L.u1;
        double v_R = U_R.u3 / U_R.u1;
        double c_L = sqrt(Gamma * p_L / rho_L);
        double c_R = sqrt(Gamma * p_R / rho_R);

        //Step 1: Estimate wave speeds. A pressure based method estimate is chosen

        //Estimate for the star region (middle wave region) pressure
        double p_pvrs = 0.5 * (p_L + p_R) + (u_L - u_R) * (rho_L + rho_R) * (c_L + c_R) / 8; //(10.67) Toro
        double p_star = std::max(0.0, p_pvrs); //Excluding unphysical pressure
        double q_L = (p_star <= p_L) ? 1 : sqrt(1 + (Gamma + 1) / (2 * Gamma) * (p_star / p_L - 1)); //(10.69) Toro
        double q_R = (p_star <= p_R) ? 1 : sqrt(1 + (Gamma + 1) / (2 * Gamma) * (p_star / p_R - 1)); //(10.69) Toro

        double S_L = u_L - c_L * q_L; //Leftmost wave speed
        double S_R = u_R + c_R * q_R; //Rightmost wave speed
        //Star region wave speed may be found analytically when the two other wave speeds are known
        double S_star = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) /
                        (rho_L * (S_L - u_L) - rho_R * (S_R - u_R));


        //(10.71) Toro
        if (S_star >= 0) {
            vec4 F_L{rho_L * u_L,
                     rho_L * u_L * u_L + p_L,
                     rho_L * u_L * v_L,
                     u_L * (U_L.u4 + p_L)};
            if (S_L >= 0) {
                return F_L;
            } else { //S_L < 0
                //(10.73)
                vec4 U_star_L = {1,
                                 S_star,
                                 v_L,
                                 U_L.u4 / rho_L + (S_star - u_L) * (S_star + p_L / (rho_L * (S_L - u_L)))};
                U_star_L *= rho_L * (S_L - u_L) / (S_L - S_star);
                return F_L + S_L * (U_star_L - U_L); //(10.72)
            }
        } else { //S_star < 0
            vec4 F_R{rho_R * u_R,
                     rho_R * u_R * u_R + p_R,
                     rho_R * u_R * v_R,
                     u_R * (U_R.u4 + p_R)};
            if (S_R > 0) {
                //(10.73)
                vec4 U_star_R = {1,
                                 S_star,
                                 v_R,
                                 U_R.u4 / rho_R + (S_star - u_R) * (S_star + p_R / (rho_R * (S_R - u_R)))};
                U_star_R *= rho_R * (S_R - u_R) / (S_R - S_star);
                return F_R + S_R * (U_star_R - U_R); //(10.72)
            } else { //S_R <= 0
                return F_R;
            }
        }
    }
    vec4 FVM_Solver::G_f_HLLC(const vec4 &U_D, const vec4 &U_U) {
        //The symmetry of the problem can be utilized to use the horizontal Riemann solver to solve the vertical flux.
        //Just switch the u and v momentum, and switch the flux components back afterwards
        vec4 tmp = F_f_HLLC({U_D.u1, U_D.u3, -U_D.u2, U_D.u4}, {U_U.u1, U_U.u3, -U_U.u2, U_U.u4});
        return {tmp.u1, -tmp.u3, tmp.u2, tmp.u4};
    }

}