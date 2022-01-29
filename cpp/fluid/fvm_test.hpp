//
// Created by anders on 1/4/22.
//

#ifndef FSIPROJECT_FVM_TEST_HPP
#define FSIPROJECT_FVM_TEST_HPP

#include "fvm_solver.hpp"
#include "../../includes.hpp"

namespace fluid{

    void set_inital_cond1(vec4* U, int ni, int nj);

    void set_initial_cond_riemann(vec4* U, int ni, int nj,const vec4& V_l, const vec4& V_r);
}


#endif //FSIPROJECT_FVM_TEST_HPP
