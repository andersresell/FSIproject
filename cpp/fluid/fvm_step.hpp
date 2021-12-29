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
    OdeScheme ode_scheme;
    FluxScheme flux_scheme;
    BoundaryConditions boundary_conitions;
    vec4 *U, *U_temporary, *V, *R, *U_left, *U_right, *U_down, *U_up, *F, *G;
    double *P, *sprad;

public:
    FVM_Step(int ni, int nj, OdeScheme ode_scheme, FluxScheme flux_scheme, BoundaryConditions boundary_conditions);
    void eval_RHS(vec4* U_in);
    void MUSCL_exstrapolate(vec4* U_in);
    void conserved2primitive(vec4* U_in);
    inline vec4 primitive2conserved(const vec4& V_in);
    inline vec4 minmod(const vec4& a, const vec4& b) const;
    void calc_P(vec4* U_in);
    inline double calc_P(const vec4& U_in) const;
    void calc_F(vec4* U_hor);
    void calc_G(vec4* U_ver);
    void ode_step(double dt);
    void set_ghost_points();

    ~FVM_Step();
};




#endif //FSIPROJECT_FVM_STEP_HPP
