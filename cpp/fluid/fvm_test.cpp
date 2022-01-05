//
// Created by anders on 1/4/22.
//

#include "fvm_test.hpp"

void FVM_Test::set_inital_cond1(vec4* U, int ni, int nj){
    //Sets the pressure in the bottom left part high and the pressure in the rest low
    double rho = 1.2;
    double p_low = 1e5;
    double p_high = 10*p_low;
    double u{0};
    double v{0};

    for (int i{0};i< ni+4;i++){
        for (int j{0};j< nj+4;j++){
            U[i].u1 = rho;
            U[i].u2 = rho* u;
            U[i].u3 = rho*v;
            if (i < ni/2 && j < nj/2){
                U[i].u4 = p_high/(Gamma - 1) + 0.5*rho*(u*u + v*v);
            }
            else{
                U[i].u4 = p_low/(Gamma - 1) + 0.5*rho*(u*u + v*v);
            }
        }
    }
    std::cout << "Initial condition 1 set\n";

}
