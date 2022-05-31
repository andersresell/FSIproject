
#include "includes.hpp"


#include "cpp/fsi/fsi_solver.hpp"
#include "cpp/fluid/fvm_utilities.hpp"
#include "cpp/solid/solid_utilities.hpp"
#include "cpp/solid/solid_body.hpp"
#include "cpp/fsi/simulate.hpp"


int main() {
    using namespace std;

    //Simulate("supersonic_wedge.json");
    //Simulate("supersonic_wedge_experimental.json");
    //Simulate("riemann_sod.json");
    //Simulate("static_wall.json");
    //Simulate("dynamic_unconstrained.json");
    //Simulate("pressure_bubbl.json");
    //Simulate("piston_prescribed.json");
    //Simulate("piston_fsi.json");
    //Simulate("piston_fsi_static.json");
    //Simulate("dynamic_circle_supersonic.json");
    //Simulate("circle_convergence.json");
    //Simulate("dynamic_circle_blast.json");
    //Simulate("tube.json");
    //Simulate("experiment_diamond_1_tube.json");
    //Simulate("experiment_diamond_1_box.json");
    //Simulate("experiment_square_3_tube.json");
    //Simulate("experiment_square_3_box.json");
    //Simulate("experiment_free_plate_al.json");
    //Simulate("experiment_free_plate_pc_tube.json");
    //Simulate("experiment_free_plate_al_tube.json");
    //Simulate("experiment_free_plate_st_tube.json");
    Simulate("experiment_free_plate_pc_box.json");
    Simulate("experiment_free_plate_pc_box.json");
    Simulate("experiment_free_plate_st_box.json");

    return 0;
}
