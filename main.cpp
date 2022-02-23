
#include "includes.hpp"
//#include "pybind11/include/pybind11/pybind11.h"

#include "cpp/fsi/fsi_solver.hpp"
#include "cpp/fluid/fvm_test.hpp"
#include "cpp/fluid/fvm_utilities.hpp"
#include "cpp/solid/solid_utilities.hpp"
#include "cpp/solid/solid_body.hpp"



int main() {
    using namespace std;
    using namespace Eigen;


    FSI_Solver::wedge_verification();
/*
    int ni = 20;
    int nj = 20;
    double L_x = 10;
    double L_y = 10;
    double CFL = 0.8;
    int n_timesteps{10};
    int fvm_write_stride{1};
    std::string output_folder{"output0"};
    fluid::FluxScheme flux_scheme{fluid::FluxScheme::HLLC};
    fluid::OdeScheme ode_scheme{fluid::OdeScheme::ExplicitEuler};
    fluid::FVM_Solver fvm{ni,nj,L_x,L_y,CFL,ode_scheme,flux_scheme,fluid::AllWalls{ni, nj}};


    vector<solid::Point> circle;
    double R = 2;
    double x0 = 5;
    double y0 = 5;
    int n = 10;
    for (int i{0}; i<n;i++){
        double theta = 2*M_PI*i/n;
        circle.push_back(solid::Point{x0+R*cos(theta),y0+R* sin(theta)});
    }

    solid::SolidBody sb{fvm,circle};


    sb.find_solid_cells();
    sb.find_ghost_cells();
    sb.find_intercepts();
    sb.interpolate_invicid_wall(fvm.U);
    sb.debug_csv();
    sb.write_boundary_csv("solid_debug_boundary");
    sb.debug_intercepts_csv();

    /*
      int ni{100};
      int nj{100};
      double L_x{10};
      double L_y{2};
      double CFL{0.8};
      int n_timesteps{500};
      std::string fvm_output_folder{"output0"};
      int write_stride{10};
      FSI_Solver::fluid_solve_test(ni,nj,L_x,L_y,CFL,n_timesteps, write_stride, fvm_output_folder,fluid::OdeScheme::TVD_RK3, fluid::FluxScheme::Rusanov);


    //Testing sod's problem
    double rho_l = 3;
    double u_l = 0;
    double v_l = 0;
    double p_l = 3;
    double rho_r = 1;
    double u_r = 0;
    double v_r = 0;
    double p_r = 1;

    //Stationary contact test
    double rho_l = 3;
    double u_l = 0;
    double v_l = 0;
    double p_l = 1;
    double rho_r = 1;
    double u_r = 0;
    double v_r = 0;
    double p_r = 1;

    double rho_l = 3;
    double u_l = 1;
    double v_l = 0;
    double p_l = 1;
    double rho_r = 1;
    double u_r = 1;
    double v_r = 0;
    double p_r = 1;


    int ni{400};
    int nj{400};
    double length_x = 4;
    double length_y = 4;
    double CFL = 0.8;
    //int n_timesteps{100};
    int fvm_write_stride{1};
    std::string output_folder{"output_riemann_sod"};
    //fluid::OdeScheme ode_scheme{fluid::OdeScheme::TVD_RK3};
    fluid::FluxScheme flux_scheme{fluid::FluxScheme::HLLC};
    fluid::OdeScheme ode_scheme{fluid::OdeScheme::ExplicitEuler};
    //fluid::FluxScheme flux_scheme{fluid::FluxScheme::Rusanov};
    fluid::vec4 V_l{rho_l, u_l, v_l, p_l};
    fluid::vec4 V_r{rho_r, u_r, v_r, p_r};
    double endtime{0.5};
    FSI_Solver::fluid_solve_riemann(ni, nj, length_x, length_y, CFL, V_l, V_r, endtime, fvm_write_stride,
                                    output_folder,ode_scheme, flux_scheme);
*/
    return 0;
}
