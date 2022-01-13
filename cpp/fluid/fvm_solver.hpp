//
// Created by anders on 12/26/21.
//

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

        void write_simple_fvm_csv_file(std::string filename);

        double ode_step();


    private:
        void eval_RHS(vec4 *U_in);

        void MUSCL_extrapolate(vec4 *U_in);

        void conserved2primitive(vec4 *U_in);

        vec4 conserved2primitive(const vec4 &U_in);

        vec4 primitive2conserved(const vec4 &V_in);

        vec4 minmod(const vec4 &a, const vec4 &b) const;

        //void calc_P(vec4* U_in);
        static double calc_P(const vec4 &U_in);

        vec4 calc_F(const vec4 &U_in) const;

        vec4 calc_G(const vec4 &U_in) const;

        void rusanov();

        double calc_sprad_x(const vec4 &U_in) const;

        double calc_sprad_y(const vec4 &U_in) const;

        double calc_sound_speed(const vec4 &U_in) const;

        double calc_timestep(double CFL) const;

        //External boundary conditions
        /*
        void set_external_BCs();
        bc_func_ptr select_bc_type(const BC& bc);
        void set_western_BCs(bc_func_ptr fptr);
        void set_eastern_BCs(bc_func_ptr fptr);
        void set_southern_BCs(bc_func_ptr fptr);
        void set_northern_BCs(bc_func_ptr fptr);
    */

        // static vec4 set_vertical_wall(const vec4& U_in);
        //static vec4 set_horizontal_wall(const vec4& U_in);

    public:
        ~FVM_Solver();
    };
}



#endif //FSIPROJECT_FVM_SOLVER_HPP
