//
// Created by anders on 12/29/21.
//

#ifndef FSIPROJECT_FSI_SOLVER_HPP
#define FSIPROJECT_FSI_SOLVER_HPP

#include "../../includes.hpp"
#include "../fluid/fvm_solver.hpp"
#include "../fluid/fvm_utilities.hpp"
#include "../solid/solid_body.hpp"
#include "../solid/dynamic_rigid.hpp"
#include "setup_cases.hpp"

enum class StoppingCrit{Time=0, Timesteps, Convergence};

struct Totals{
    double density_L2_norm, mass;
};

class FSI_Solver{
private: //remove later
    fluid::FVM_Solver& fvm;
    std::vector<std::shared_ptr<solid::SolidBody>> solid_bodies;
    std::pair<StoppingCrit,double> stopping_crit;
    const int fvm_write_stride;
    const std::string output_folder;

    //For convergence check and output
    double* rho_old;
    std::vector<Totals> totals_history;
public:

    void add_solid(std::shared_ptr<solid::SolidBody>&& solid_body);

    FSI_Solver(fluid::FVM_Solver& fvm, int fvm_write_stride, std::string fvm_output_folder);

    int solve();

    //use any either of the two functions to define the stopping criterion of the simulation
    void set_timesteps(int n_timesteps) {stopping_crit = {StoppingCrit::Timesteps, n_timesteps};}
    void set_endtime(double t_end) {stopping_crit = {StoppingCrit::Time, t_end};}
    void set_convergence(double tol=1e-6) {stopping_crit = {StoppingCrit::Convergence, tol};}

    //Convergence check by L2 norm of density
    double calc_density_L2_norm();
    double calc_mass();
    void set_rho_old();
    void write_totals_history();

    void write_fsi_header();
    void write_static_solid_boundaries();
    void write_movable_solid_boundaries(int n);
    void write_solid_debug_files(int n);


    ~FSI_Solver();
};





#endif //FSIPROJECT_FSI_SOLVER_HPP
