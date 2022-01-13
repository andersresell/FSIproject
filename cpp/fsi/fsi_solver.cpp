//
// Created by anders on 12/29/21.
//

#include "fsi_solver.hpp"


FSI_Solver::FSI_Solver(fluid::FVM_Solver& fvm, int write_stride)
    : fvm{fvm}, write_stride{write_stride}
{
}

int FSI_Solver::solve(){
    int n{0};
    double t{0};
    double dt;
    while (true) {

        std::cout << "FSI solve: n = " + std::to_string(n) + "\n";
        if (n % write_stride == 0) {
            std::cout << "Writing FVM output\n";
            fvm.write_simple_fvm_csv_file("fvm_out_t" + std::to_string(n) + ".csv");
        }
        dt = fvm.ode_step();


        if (stopping_crit.first == StoppingCrit::Time){
            if (stopping_crit.second >= t) break;
        }
        else if (stopping_crit.first == StoppingCrit::Timesteps){
            if (n >= (int)stopping_crit.second) break;
        }
        else {
            std::cerr << "Invalid stopping criterion\n";
            exit(1);
        }
        t += dt;
        n++;
    }
    return n;
}

void FSI_Solver::fluid_solve(int ni,
                 int nj,
                 double L_x,
                 double L_y,
                 double CFL,
                 int n_timesteps,
                 int write_stride) {
    //Should also add functionality for choosing the time scheme, flux scheme, boundaries and initital cond
    fluid::FVM_Solver fvm{ni, nj, L_x, L_y, CFL, fluid::OdeScheme::ExplicitEuler, fluid::FluxScheme::Rusanov,
                          fluid::AllWalls{ni, nj}};

    fluid::set_inital_cond1(fvm.U, fvm.ni, fvm.nj);
    FSI_Solver fsi{fvm, write_stride};
    fsi.set_timesteps(n_timesteps);
    fsi.solve();
}
