
#include "includes.hpp"
//#include "pybind11/include/pybind11/pybind11.h"

#include "cpp/fsi/fsi_solver.hpp"
#include "cpp/fluid/fvm_utilities.hpp"
#include "cpp/solid/solid_utilities.hpp"
#include "cpp/solid/solid_body.hpp"
#include "cpp/fsi/simulate.hpp"


//#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
int main() {
    using namespace std;

    //Simulate("supersonic_wedge.json");
    //Simulate("supersonic_wedge_experimental.json");
    //Simulate("riemann_sod.json");
    //Simulate("static_wall.json");
    //Simulate("dynamic_unconstrained.json");
    //Simulate("pressure_bubbl.json");
    //Simulate("piston_prescribed.json");
    Simulate("piston_fsi.json");
    //Simulate("dynamic_circle.json");
    //Simulate("tube.json");
    return 0;
}
