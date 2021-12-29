//
// Created by anders on 12/26/21.
//

#include "fvm_step.hpp"

FVM_Step::FVM_Step(int ni, int nj, OdeScheme ode_scheme, FluxScheme flux_scheme, BoundaryConditions boundary_conditions)
: ni{ni}, nj{nj}, ode_scheme{ode_scheme}, flux_scheme{flux_scheme}
{

}

void FVM_Step::eval_RHS(vec4* U_in){
    //eval U_left, U_right etc
    MUSCL_exstrapolate(U_in);
}

void FVM_Step::MUSCL_exstrapolate(vec4* U_in){
    vec4 tmpV;
    conserved2primitive(U_in);
    for (int i{1};i<ni+2;i++){
        for (int j{2};j<nj+2;j++) {
            tmpV = V[IX(i, j)] + minmod(V[IX(i, j)] - V[IX(i - 1, j)], V[IX(i + 1, j)] - V[IX(i, j)])*0.5;
            U_left[IXH(i, j)] = primitive2conserved(tmpV);

            tmpV = V[IX(i + 1, j)] - minmod(V[IX(i + 2, j)] - V[IX(i + 1, j)], V[IX(i + 1, j)] - V[IX(i, j)])*0.5;
            U_right[IXH(i, j)] = primitive2conserved(tmpV);
        }
    }
    for (int i{2};i< ni+2;i++) {
        for (int j{1}; j < nj + 2; j++) {
            tmpV = V[IX(i, j)] + minmod(V[IX(i, j)] - V[IX(i, j-1)], V[IX(i, j+1)] - V[IX(i, j)])*0.5;
            U_down[IXV(i,j)] = primitive2conserved(tmpV);

            tmpV = V[IX(i, j)] - minmod(V[IX(i, j+2)] - V[IX(i, j+1)], V[IX(i, j+1)] - V[IX(i, j)])*0.5;
            U_up[IXV(i,j)] = primitive2conserved(tmpV);
        }
    }
}

inline vec4 FVM_Step::minmod(const vec4& a, const vec4& b) const{
    return {sgn(a.u1)*max(0,min(abs(a.u1),sgn(a.u1)*b.u1)),
             sgn(a.u2)*max(0,min(abs(a.u2),sgn(a.u2)*b.u2)),
             sgn(a.u3)*max(0,min(abs(a.u3),sgn(a.u3)*b.u3)),
             sgn(a.u4)*max(0,min(abs(a.u4),sgn(a.u4)*b.u4))};
}

void FVM_Step::conserved2primitive(vec4* U_in){
    /*
    for (int i{0};i<ni+4;i++){
        for (int j{0};j<nj+4;j++){
            V[IX(i,j)].u1 = U_in[IX(i,j)].u1;
            V[IX(i,j)].u2 = U_in[IX(i,j)].u2/U_in[IX(i,j)].u1;
            V[IX(i,j)].u3 = U_in[IX(i,j)].u3/U_in[IX(i,j)].u1;
            V[IX(i,j)].u4 = calc_P(U_in[IX(i,j)]);
        }
    }
     */
    for (int i{0};i<(ni+4)*(nj+4);i++){
        V[i].u1 = U_in[i].u1;
        V[i].u2 = U_in[i].u2/U_in[i].u1;
        V[i].u3 = U_in[i].u3/U_in[i].u1;
        V[i].u4 = calc_P(U_in[i]);
    }
}

void FVM_Step::calc_P(vec4* U_in){
    for (int i{0};i<ni+4;i++){
        for (int j{0};j<nj+4;j++) {
            // p = (gamma-1)* (U4 - 0.5*(U2² + U3²)/U1
            P[XI(i, j)] = (Gamma - 1) *
                          (U_in[XI(i, j)].u4 - 0.5 * (U_in[XI(i, j)].u2 * U_in[XI(i, j)].u2 + U_in[XI(i, j)].u3 * U_in[XI(i, j)].u3)) /
                          U_in[XI(i, j)].u1;
        }
    }
}

inline double FVM_Step::calc_P(const vec4& U_in) const{
    return (Gamma - 1) * (U_in.u4 - 0.5 * (U_in.u2 * U_in.u2 + U_in.u3 * U_in.u3)) /U_in.u1;
}

vec4 FVM_Step::primitive2conserved(const vec4& V_in) {
    return {V_in.u1,
            V_in.u2 * V_in.u1,
            V_in.u3 * V_in.u1,
            V_in.u4 / (Gamma - 1) + 0.5 * V_in.u1 * (V_in.u2 * V_in.u2 + V_in.u3 * V_in.u3)};
}

void FVM_Step::calc_F(vec4* U_hor) {
    double p;
    vec4* U_ptr;
    for (int i{0}; i < ni + 1; i++) {
        for (int j{0}; j < nj; j++) {
            U_ptr = &U_in[IXH(i, j)]; //to avoid multiple lookups
            p = calc_P(*U_ptr);
            F[IXH(i,j)] = {U_ptr->u2,
                      U_ptr->u2 * U_ptr->u2 / U_ptr->u1 + p,
                      U_ptr->u2 * U_ptr->u3 / U_ptr->u1,
                      (U_ptr->u4 + p) * U_ptr->u2 / U_ptr->u1};
        }
    }
    U_ptr = nullptr;
}

void FVM_Step::calc_G(vec4* U_ver) {
    double p;
    vec4* U_ptr;
    for (int i{0}; i < ni; i++) {
        for (int j{0}; j < nj + 1; j++) {
            U_ptr = &U_in[IXH(i, j)]; //to avoid multiple lookups
            p = calc_P(*U_ptr);
            F[IXH(i,j)] = {U_ptr->u3,
                      U_ptr->u2 * U_ptr->u3 / U_ptr->u1,
                      U_ptr->u3 * U_ptr->u3 / U_ptr->u1 + p,
                      (U_ptr->u4 + p) * U_ptr->u3 / U_ptr->u1};
        }
    }
    U_ptr = nullptr;
}


FVM_Step::~FVM_Step(){
    delete[] U, U_temporary, V, R, U_left, U_right, U_down, U_up, F, G;
    delete[] P, sprad;
}