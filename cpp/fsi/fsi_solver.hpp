//
// Created by anders on 12/29/21.
//

#ifndef FSIPROJECT_FSI_SOLVER_HPP
#define FSIPROJECT_FSI_SOLVER_HPP

#include "../fluid/fvm_solver.hpp"
#include "../fluid/fvm_utilities.hpp"
#include "../fluid/fvm_test.hpp"
#include "../../includes.hpp"

enum class StoppingCrit{Time=0, Timesteps};

class FSI_Solver{
    fluid::FVM_Solver& fvm;
    std::pair<StoppingCrit,double> stopping_crit;
    const int fvm_write_stride;
    const std::string fvm_output_folder;
public:
    FSI_Solver(fluid::FVM_Solver& fvm, int write_stride, std::string fvm_output_folder);

    int solve();

    //use any either of the two functions to define the stopping criterion of the simulation
    void set_timesteps(int n_timesteps) {stopping_crit = {StoppingCrit::Timesteps, n_timesteps};}
    void set_endtime(double t_end) {stopping_crit = {StoppingCrit::Time, t_end};}


    static void fluid_solve_test(int ni,
                     int nj,
                     double L_x,
                     double L_y,
                     double CFL,
                     int n_timesteps,
                     int fvm_write_stride,
                     std::string fvm_output_folder,
                     fluid::OdeScheme ode_scheme,
                     fluid::FluxScheme flux_scheme);

    static void fluid_solve_riemann(int ni,
                                 int nj,
                                 double L_x,
                                 double L_y,
                                 double CFL,
                                 const fluid::vec4& V_l,
                                 const fluid::vec4& V_r,
                                 double entime,
                                 int fvm_write_stride,
                                 std::string fvm_output_folder,
                                 fluid::OdeScheme ode_scheme,
                                 fluid::FluxScheme flux_scheme);
};



#endif //FSIPROJECT_FSI_SOLVER_HPP
