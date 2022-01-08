//
// Created by anders on 12/29/21.
//

#ifndef FSIPROJECT_FSI_SOLVER_HPP
#define FSIPROJECT_FSI_SOLVER_HPP

#include "../fluid/fvm_solver.hpp"
#include "../fluid/fvm_utilities.hpp"
#include "../../includes.hpp"

constexpr unsigned int MAX_TIMESTEPS{10};

class FSI_Solver{
    FVM_Solver& fvm;

public:
    FSI_Solver(FVM_Solver& fvm);

    void solve(double t_end);

};

#endif //FSIPROJECT_FSI_SOLVER_HPP
