#include <iostream>

#include "cpp/fsi/fsi_solver.hpp"
#include "cpp/fluid/fvm_test.hpp"
#include "cpp/fluid/fvm_utilities.hpp"



int main() {
    std::cout << "Hello, World!" << std::endl;

    int ni{2};
    int nj{2};
    double Lx{1};
    double Ly{1};
    double CFL{0.5};
    FVM_Solver fvm{ni,nj,Lx,Ly,CFL, OdeScheme::ExplicitEuler, FluxScheme::Rusanov, AllWalls{ni,nj}};
    FVM_Test::set_inital_cond1(fvm.U,fvm.ni,fvm.nj);


    //write_simple_fvm_csv_file("testfile.csv",fvm.U,fvm.ni,fvm.nj);
    //fvm.ode_step();
    //std::cout << "hello\n";


    FSI_Solver fsi{fvm};
    fsi.solve(1);
    //fsi.solve(1);

}
