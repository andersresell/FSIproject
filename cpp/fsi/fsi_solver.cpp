//
// Created by anders on 12/29/21.
//

#include "fsi_solver.hpp"


FSI_Solver::FSI_Solver(FVM_Solver& fvm)
    : fvm{fvm}
{
}

void FSI_Solver::solve(double t_end){
    double t{0};
    double dt;
    while (t < t_end){
        dt = fvm.ode_step();
        t += dt;
        std::cout << "FSI_Solver::solve() t = " << t << "\n";
    }
}