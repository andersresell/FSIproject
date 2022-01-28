//
// Created by anders on 12/26/21.
//
//A carthesian FVM solver employing a semi discrete approach, using an ODE solver. It uses MUSCL extrapolation
// with minmod limiter to get second order accuracy and an approximate Riemann solver to compute the fluxes.
// The boundaries are handled by ghost points

#ifndef FSIPROJECT_FVM_SOLVER_HPP
#define FSIPROJECT_FVM_SOLVER_HPP

#include "fvm_utilities.hpp"
#include "../../includes.hpp"

namespace fluid {
//typedef vec4 (*bc_func_ptr)(const vec4& U_in); //In order to set appropriate BC's functions are passed as arguents and return types

    constexpr double Gamma = 1.4;

    class FVM_Solver {
        //Computes one timestep of the fluid domain with the explicit euler method
    public:
        const int ni, nj;
        const double L_x, L_y;
        const double dx, dy;
        const double CFL;
        OdeScheme ode_scheme;
        FluxScheme flux_scheme;
        ExternalBCs external_bcs;
        vec4 *U;
    private:
        vec4 *U_tmp, *V, *Res, *U_left, *U_right, *U_down, *U_up, *F_f, *G_f;

    public:
        FVM_Solver(int ni, int nj, double L_x, double L_y, double CFL, OdeScheme ode_scheme, FluxScheme flux_scheme,
                   ExternalBCs external_bcs);

        void write_simple_fvm_csv_file(const std::string& filename);

        double ode_step();

    private:
        double calc_timestep(double CFL) const;

        void eval_RHS(vec4 *U_in);

        void MUSCL_extrapolate(vec4 *U_in);

        void conserved2primitive(vec4 *U_in);

        //vec4 conserved2primitive(const vec4 &U_in);

        static vec4 primitive2conserved(const vec4 &V_in);

        void rusanov();

        void HLLC();

        static vec4 F_f_HLLC(const vec4& U_L,const vec4& U_R); //Computes the HLLC flux F_{i+1/2,j}

        static vec4 G_f_HLLC(const vec4& U_D,const vec4& U_U); //Computes the HLLC flux G_{i,j+1/2}

        static double calc_P(const vec4 &U_in);

        static vec4 calc_F(const vec4 &U_in);

        static vec4 calc_G(const vec4 &U_in);

        static double calc_sound_speed(const vec4 &U_in);

        static vec4 minmod(const vec4 &a, const vec4 &b);

        static double calc_sprad_x(const vec4 &U_in);

        static double calc_sprad_y(const vec4 &U_in);

    public:
        ~FVM_Solver();
    };


    /*
    inline vec4 FVM_Solver::conserved2primitive(const vec4 &U_in) {
        return {U_in.u1,
                U_in.u2 / U_in.u1,
                U_in.u3 / U_in.u1,
                calc_P(U_in)};
    }*/

    inline double FVM_Solver::calc_P(const vec4 &U_in) {
        return (Gamma - 1) * (U_in.u4 - 0.5 * (U_in.u2 * U_in.u2 + U_in.u3 * U_in.u3) / U_in.u1);
    }

    inline vec4 FVM_Solver::primitive2conserved(const vec4 &V_in) {
        return {V_in.u1,
                V_in.u2 * V_in.u1,
                V_in.u3 * V_in.u1,
                V_in.u4 / (Gamma - 1) + 0.5 * V_in.u1 * (V_in.u2 * V_in.u2 + V_in.u3 * V_in.u3)};
    }

    inline vec4 FVM_Solver::calc_F(const vec4 &U_in){
        double P = calc_P(U_in);
        return {U_in.u2,
                U_in.u2 * U_in.u2 / U_in.u1 + P,
                U_in.u2 * U_in.u3 / U_in.u1,
                (U_in.u4 + P) * U_in.u2 / U_in.u1};
    }

    inline vec4 FVM_Solver::calc_G(const vec4 &U_in) {
        double P = calc_P(U_in);
        return {U_in.u3,
                U_in.u2 * U_in.u3 / U_in.u1,
                U_in.u3 * U_in.u3 / U_in.u1 + P,
                (U_in.u4 + P) * U_in.u3 / U_in.u1};
    }

    inline double FVM_Solver::calc_sound_speed(const vec4 &U_in) {
        return sqrt(Gamma / U_in.u1 * calc_P(U_in));
    }

    inline double FVM_Solver::calc_sprad_x(const vec4 &U_in) {
        return std::abs(U_in.u2 / U_in.u1) + calc_sound_speed(U_in);
    }

    inline double FVM_Solver::calc_sprad_y(const vec4 &U_in) {
        return std::abs(U_in.u3 / U_in.u1) + calc_sound_speed(U_in);
    }



}
#endif //FSIPROJECT_FVM_SOLVER_HPP
