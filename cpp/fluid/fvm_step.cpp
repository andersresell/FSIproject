//
// Created by anders on 12/26/21.
//

#include "fvm_step.hpp"

FVM_Step::FVM_Step(int ni, int nj, OdeScheme ode_scheme, FluxScheme flux_scheme, BoundaryConditions boundary_conditions)
: ni{ni}, nj{nj}, ode_scheme{ode_scheme}, flux_scheme{flux_scheme}
{

}

void FVM_Step::eval_RHS(vec4* U){
    //eval U_left, U_right etc
    MUSCL_exstrapolate(U);
}

void FVM_Step::MUSCL_exstrapolate(vec4* U){
    vec4 tmpV;
    conserved2primitive(U);
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

void FVM_Step::conserved2primitive(vec4* U){
    for (int i{0};i<ni+4;i++){
        for (int j{0};j<nj+4;j++){
            V[IX(i,j)].u1 = U[IX(i,j)].u1;
            V[IX(i,j)].u2 = U[IX(i,j)].u2/U[IX(i,j)].u1;
            V[IX(i,j)].u3 = U[IX(i,j)].u3/U[IX(i,j)].u1;
            V[IX(i,j)].u4 = calc_P(U[IX(i,j)]);
        }
    }
}

void FVM_Step::calc_P(vec4* U){
    for (int i{0};i<ni+4;i++){
        for (int j{0};j<nj+4;j++) {
            // p = (gamma-1)* (U4 - 0.5*(U2² + U3²)/U1
            P[XI(i, j)] = (Gamma - 1) *
                          (U[XI(i, j)].u4 - 0.5 * (U[XI(i, j)].u2 * U[XI(i, j)].u2 + U[XI(i, j)].u3 * U[XI(i, j)].u3)) /
                          U[XI(i, j)].u1;
        }
    }
}

inline double FVM_Step::calc_P(vec4& U) const{
    return (Gamma - 1) * (U.u4 - 0.5 * (U.u2 * U.u2 + U.u3 * U.u3)) /U.u1;
}

vec4 FVM_Step::primitive2conserved(const vec4& V) {
    return {V.u1,
            V.u2 * V.u1,
            V.u3 * V.u1,
            V.u4 / (Gamma - 1) + 0.5 * V.u1 * (V.u2 * V.u2 + V.u3 * V.u3)};
}

void FVM_Step::calc_F(vec4* U_hor) {
    double p;
    uint8_t ind;
    vec4 U_tmp;
    for (int i{0}; i < ni + 1; i++) {
        for (int j{0}; j < nj + 1; j++) {
            ind = IXH(i, j);
            U_tmp = U_hor[ind];
            p = calc_P(U_tmp);
            F[ind] = {U_tmp.u2,
                      U_tmp.u2 * U_tmp.u2 / U_tmp.u1 + p,
                      U_tmp.u2 * U_tmp.u3 / U_tmp.u1,
                      (U_tmp.u4 + p) * U_tmp.u2 / U_tmp.u1};
        }
    }
}
void FVM_Step::calc_G(vec4* U_hor){

}

FVM_Step::~FVM_Step(){
    delete[] U, U_temporary, R, U_left, U_right, U_down, U_up, F, G;
    delete[] P, sprad;
}