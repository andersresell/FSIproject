#include <iostream>

#include "cpp/fsi/fsi_solver.hpp"
#include "cpp/fluid/fvm_test.hpp"


int main() {
    std::cout << "Hello, World!" << std::endl;

    int ni{10};
    int nj{10};
    double Lx{10};
    double Ly{10};
    double CFL{0.5};
    FVM_Solver fvm{ni,nj,Lx,Ly,CFL, OdeScheme::ExplicitEuler, FluxScheme::Rusanov, AllWalls{ni,nj}};
    FVM_Test::set_inital_cond1(fvm.U,fvm.ni,fvm.nj);
    fvm.ode_step();
    std::cout << "hello\n";


    FSI_Solver fsi{fvm};
    //fsi.solve(1);

}
