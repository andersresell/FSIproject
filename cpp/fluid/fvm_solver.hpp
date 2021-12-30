//
// Created by anders on 12/26/21.
//

#ifndef FSIPROJECT_FVM_SOLVER_HPP
#define FSIPROJECT_FVM_SOLVER_HPP

#include "fvm_utilities.hpp"
constexpr double Gamma = 1.4;

class FVM_Solver{
    //Computes one timestep of the fluid domain with the explicit euler method
    const int ni, nj;
    double dx, dy;
    double CFL;
    OdeScheme ode_scheme;
    FluxScheme flux_scheme;
    ExternalBCs external_bcs;
    vec4 *U, *U_tmp, *V, *Res, *U_left, *U_right, *U_down, *U_up, *F_f, *G_f;

public:
    FVM_Solver(int ni, int nj, double L_x, double L_y, double CFL, OdeScheme ode_scheme, FluxScheme flux_scheme,
               ExternalBCs external_bcs);
    double ode_step();

private:
    void eval_RHS(vec4* U_in);
    void MUSCL_extrapolate(vec4* U_in);
    void conserved2primitive(vec4* U_in);
    vec4 primitive2conserved(const vec4& V_in);
    vec4 minmod(const vec4& a, const vec4& b) const;
    //void calc_P(vec4* U_in);
    double calc_P(const vec4& U_in) const;
    vec4 calc_F(const vec4& U_in) const;
    vec4 calc_G(const vec4& U_in) const;
    void rusanov();
    double calc_sprad_x(const vec4& U_in) const;
    double calc_sprad_y(const vec4& U_in) const;
    double calc_sound_speed(const vec4& U_in) const;
    void set_external_boundary_conditions();
    double calc_timestep(double CFL) const;

public:

    ~FVM_Solver();
};




#endif //FSIPROJECT_FVM_SOLVER_HPP
