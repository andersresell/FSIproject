#include <iostream>

#include "cpp/fsi/fsi_solver.hpp"


int main() {
    std::cout << "Hello, World!" << std::endl;

    FVM_Solver fvm{10,10,10,10,0.5,OdeScheme::ExplicitEuler,FluxScheme::Rusanov,ExternalBCs{}};

    FSI_Solver fsi{fvm};
    fsi.solve(1);

}
