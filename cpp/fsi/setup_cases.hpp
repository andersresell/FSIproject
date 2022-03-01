//
// Created by anders on 2/23/22.
//

#ifndef FSIPROJECT_SETUP_CASES_HPP
#define FSIPROJECT_SETUP_CASES_HPP

#include "../fluid/fvm_solver.hpp"
#include "../solid/solid_body.hpp"
#include "../../includes.hpp"

namespace fluid{

    void set_initial_cond1(vec4* U, int ni, int nj);

    void set_initial_cond2(vec4* U, int ni, int nj);

    void set_constant_horizontal_flow_cond(vec4* U, int ni, int nj, double M = 0, double rho = 1.2, double p = 1e5);

    void set_initial_cond_riemann(vec4* U, int ni, int nj,const vec4& V_l, const vec4& V_r);

    void set_initial_cond_pressure_bubble(vec4* U, int ni, int nj, double L_x, double L_y);
}

namespace solid {
    std::vector<solid::Point>
    generate_wedge(double l, double half_angle_deg, double x_center, double y_center);
    std::vector<solid::Point> generate_circle(double R, int n_nodes, double x_center, double y_center);
}


#endif //FSIPROJECT_SETUP_CASES_HPP
