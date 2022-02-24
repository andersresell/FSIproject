
#include "includes.hpp"
//#include "pybind11/include/pybind11/pybind11.h"

#include "cpp/fsi/fsi_solver.hpp"
#include "cpp/fluid/fvm_test.hpp"
#include "cpp/fluid/fvm_utilities.hpp"
#include "cpp/solid/solid_utilities.hpp"
#include "cpp/solid/solid_body.hpp"
#include "cpp/fsi/simulate.hpp"

int main() {
    //using namespace std;
    //using namespace Eigen;


    Simulate("riemann_sod.json");
    //FSI_Solver::wedge_verification();

    return 0;
}
