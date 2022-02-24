//
// Created by anders on 12/29/21.
//

#ifndef FSIPROJECT_FSI_SOLVER_HPP
#define FSIPROJECT_FSI_SOLVER_HPP

#include "../../includes.hpp"
#include "../fluid/fvm_solver.hpp"
#include "../fluid/fvm_utilities.hpp"
#include "../fluid/fvm_test.hpp"
#include "../solid/solid_body.hpp"

enum class StoppingCrit{Time=0, Timesteps, Convergence};

class FSI_Solver{
private: //remove later
    fluid::FVM_Solver& fvm;
    std::vector<std::shared_ptr<solid::SolidBody>> solid_bodies;
    std::pair<StoppingCrit,double> stopping_crit;
    const int fvm_write_stride;
    const std::string output_folder;

    //For convergence check and output
    double* rho_old;
    std::vector<double> convergence_history;
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
    void set_rho_old();
    void write_fvm_convergence_history();

    void write_fsi_header();
    void write_static_solid_boundaries();
    void write_movable_solid_boudnaries(int n);

    static void riemann_test();

    static void solid_test();

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

    static void wedge_verification();

    ~FSI_Solver();
};



#endif //FSIPROJECT_FSI_SOLVER_HPP
