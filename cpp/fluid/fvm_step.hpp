//
// Created by anders on 12/26/21.
//

#ifndef FSIPROJECT_FVM_STEP_HPP
#define FSIPROJECT_FVM_STEP_HPP

#include "fvm_utilities.hpp"
double Gamma = 1.4;

class FVM_Step{
    //Computes one timestep of the fluid domain with the explicit euler method
    const int ni, nj;
    double dx, dy;
    OdeScheme ode_scheme;
    FluxScheme flux_scheme;
    BoundaryConditions boundary_conitions;
    vec4 *U, *U_temporary, *V, *Res, *U_left, *U_right, *U_down, *U_up, *F_f, *G_f;
    double *P, *sprad;

public:
    FVM_Step(int ni, int nj, double L_x, double L_y, OdeScheme ode_scheme, FluxScheme flux_scheme, BoundaryConditions boundary_conditions);
    void ode_step();
    void eval_RHS(vec4* U_in);
    void MUSCL_extrapolate(vec4* U_in);
    void conserved2primitive(vec4* U_in);
    inline vec4 primitive2conserved(const vec4& V_in);
    inline vec4 minmod(const vec4& a, const vec4& b) const;
    void calc_P(vec4* U_in);
    inline double calc_P(const vec4& U_in) const;
    inline vec4 calc_F(const vec4& U_in) const;
    inline vec4 calc_G(const vec4& U_in) const;
    void rusanov();
    inline double calc_sprad_x(const vec4& U_in) const;
    inline double calc_sprad_y(const vec4& U_in) const;
    inline double calc_sound_speed(const vec4& U_in) const;
    void set_ghost_points();
    double calc_timestep(double CFL) const;

    ~FVM_Step();
};




#endif //FSIPROJECT_FVM_STEP_HPP
