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
    int write_stride;
public:
    FSI_Solver(fluid::FVM_Solver& fvm, int write_stride);

    int solve();

    //use any either of the two functions to define the stopping criterion of the simulation
    void set_timesteps(int n_timesteps) {stopping_crit = {StoppingCrit::Timesteps, n_timesteps};}
    void set_time(double t_end) {stopping_crit = {StoppingCrit::Time, t_end};}


    static void fluid_solve(int ni,
                     int nj,
                     double L_x,
                     double L_y,
                     double CFL,
                     int n_timesteps,
                     int write_stride,
                     fluid::FluxScheme flux_scheme);
};



#endif //FSIPROJECT_FSI_SOLVER_HPP
