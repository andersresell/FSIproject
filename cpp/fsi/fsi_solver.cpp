//
// Created by anders on 12/29/21.
//

#include "fsi_solver.hpp"


FSI_Solver::FSI_Solver(fluid::FVM_Solver& fvm, int fvm_write_stride, std::string fvm_output_folder)
    : fvm{fvm}, fvm_write_stride{fvm_write_stride}, fvm_output_folder{std::move(fvm_output_folder)}
{
}

int FSI_Solver::solve() {
    auto start_time{std::chrono::high_resolution_clock::now()}; //Start timing
    int n{0};
    double t{0};
    double dt;
    bool breaker{false};
    while (true) {
        std::cout << "FSI solve: n = " + std::to_string(n) + "\n";
        if (n % fvm_write_stride == 0) {
            std::cout << "Writing FVM output\n";
            fvm.write_fvm_csv_out_file(fvm_output_folder, n);
        }

        dt = fvm.ode_step(); //Time stepping the fvm simulation
        std::cout << "dt = " << dt << '\n';

        if (stopping_crit.first == StoppingCrit::Time) {
            if (t >= stopping_crit.second) breaker = true;
        } else if (stopping_crit.first == StoppingCrit::Timesteps) {
            if (n >= (int) stopping_crit.second) breaker = true;
        } else {
            std::cerr << "Invalid stopping criterion\n";
            exit(1);
        }
        //Writing fvm header at last timestep, since the number of steps can't be known a priori if endtime is
        //used as stopping criterion
        if (breaker) {
            fvm.write_fvm_csv_header_file(fvm_output_folder, fvm_write_stride, n, t);
            break;
        }
        t += dt;
        n++;
    }
    auto stop_time{std::chrono::high_resolution_clock::now()};
    auto simulation_time{std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time)};
    std::cout << "End time = " << t << '\n';
    std::cout << "Computation time = " << 1.0e-3 * static_cast<double>(simulation_time.count()) << " seconds\n";
    return n;
}

void FSI_Solver::fluid_solve_test(int ni,
                 int nj,
                 double L_x,
                 double L_y,
                 double CFL,
                 int n_timesteps,
                 int fvm_write_stride,
                 std::string fvm_output_folder,
                 fluid::OdeScheme ode_scheme,
                 fluid::FluxScheme flux_scheme) {
    //Should also add functionality for choosing the time scheme, flux scheme, boundaries and initital cond
    fluid::FVM_Solver fvm{ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, fluid::AllWalls{ni, nj}};
    fluid::set_inital_cond1(fvm.U, fvm.ni, fvm.nj);
    FSI_Solver fsi{fvm, fvm_write_stride, fvm_output_folder};
    fsi.set_timesteps(n_timesteps);
    fsi.solve();
}

void FSI_Solver::fluid_solve_riemann(int ni,
                                int nj,
                                double L_x,
                                double L_y,
                                double CFL,
                                const fluid::vec4& V_l,
                                const fluid::vec4& V_r,
                                double endtime,
                                int fvm_write_stride,
                                std::string fvm_output_folder,
                                fluid::OdeScheme ode_scheme,
                                fluid::FluxScheme flux_scheme){

    fluid::FVM_Solver fvm{ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, fluid::AllWalls{ni, nj}};
    fluid::set_initial_cond_riemann(fvm.U,fvm.ni,fvm.nj, V_l, V_r);
    FSI_Solver fsi{fvm, fvm_write_stride, fvm_output_folder};
    //fsi.set_timesteps(n_timesteps);
    fsi.set_endtime(endtime);
    fsi.solve();
}