//
// Created by anders on 12/29/21.
//

#include "fsi_solver.hpp"


FSI_Solver::FSI_Solver(FVM_Solver& fvm)
    : fvm{fvm}
{
}

void FSI_Solver::solve(double t_end){
    int n{0};
    double t{0};
    double dt;
    while (t < t_end) {
        std::cout << "FSI solve: n = " + std::to_string(n) + "\n";
        write_simple_fvm_csv_file("fvm_out_t" + std::to_string(n) + ".csv", fvm.U, fvm.ni, fvm.nj);
        n++;
        dt = fvm.ode_step();
        t += dt;
        if (n > MAX_TIMESTEPS) break;
    }
}