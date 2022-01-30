
#include "includes.hpp"
//#include "pybind11/include/pybind11/pybind11.h"

#include "cpp/fsi/fsi_solver.hpp"
#include "cpp/fluid/fvm_test.hpp"
#include "cpp/fluid/fvm_utilities.hpp"



int main() {
    std::cout << "Hello, World!" << std::endl;
    /*  int ni{100};
      int nj{100};
      double L_x{10};
      double L_y{2};
      double CFL{0.8};
      int n_timesteps{500};
      std::string fvm_output_folder{"output0"};
      int write_stride{10};
      FSI_Solver::fluid_solve_test(ni,nj,L_x,L_y,CFL,n_timesteps, write_stride, fvm_output_folder,fluid::OdeScheme::TVD_RK3, fluid::FluxScheme::Rusanov);
  */
    //Testing sod's problem
    double rho_l = 3;
    double u_l = 0;
    double v_l = 0;
    double p_l = 3;
    double rho_r = 1;
    double u_r = 0;
    double v_r = 0;
    double p_r = 1;


    int ni{100};
    int nj{5};
    double length_x = 4;
    double length_y = 4;
    double CFL = 0.5;
    //int n_timesteps{100};
    int fvm_write_stride{1};
    std::string output_folder{"output_riemann_sod"};
    fluid::OdeScheme ode_scheme{fluid::OdeScheme::TVD_RK3};
    fluid::FluxScheme flux_scheme{fluid::FluxScheme::HLLC};
    //fluid::OdeScheme ode_scheme{fluid::OdeScheme::ExplicitEuler};
    //fluid::FluxScheme flux_scheme{fluid::FluxScheme::Rusanov};
    fluid::vec4 V_l{rho_l, u_l, v_l, p_l};
    fluid::vec4 V_r{rho_r, u_r, v_r, p_r};
    double endtime{1};
    FSI_Solver::fluid_solve_riemann(ni, nj, length_x, length_y, CFL, V_l, V_r, endtime, fvm_write_stride,
                                    output_folder,ode_scheme, flux_scheme);

    return 0;
}
