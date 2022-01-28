
#include "includes.hpp"
//#include "pybind11/include/pybind11/pybind11.h"

#include "cpp/fsi/fsi_solver.hpp"
#include "cpp/fluid/fvm_test.hpp"
#include "cpp/fluid/fvm_utilities.hpp"



int main() {
    std::cout << "Hello, World!" << std::endl;
    int ni{100};
    int nj{100};
    double L_x{10};
    double L_y{2};
    double CFL{0.8};
    int n_timesteps{500};
    int write_stride{10};
    FSI_Solver::fluid_solve(ni,nj,L_x,L_y,CFL,n_timesteps, write_stride,fluid::OdeScheme::TVD_RK3, fluid::FluxScheme::Rusanov);


}
