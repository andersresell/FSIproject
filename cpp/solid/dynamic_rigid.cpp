//
// Created by anders on 3/11/22.
//

#include "dynamic_rigid.hpp"

using namespace std;

namespace solid {

    DynamicRigid::DynamicRigid(fluid::FVM_Solver &fvm, std::vector<Point> &&boundary_in, Point CM, double M, double I)
            : SolidBody(fvm, std::move(boundary_in), SolidBodyType::Dynamic), M{M}, I{I} {
        y[0] = CM.x;
        y[1] = CM.y;
        y[2] = 0;//-200;//REMOVE LATER
        y[3] = 0;
        y[4] = 0;
        y[5] = 0;
        F_solid = {0, 0}; //Change later
        r0 = new Point[n_bound];
        for (int i{0}; i < n_bound; i++) {
            r0[i] = boundary[i] - CM;
        }
    }

    double DynamicRigid::max_boundary_speed() const {
        double tmp;
        double v_norm_max{0};
        double omega = y[5];
        for (int i{0}; i < n_bound; i++) {
            Point r = {boundary[i].x - y[0], boundary[i].y - y[1]};
            //v = v_CM + omega x r
            tmp = sqrt(squared(y[2] - omega * r.y) + squared(y[3] + omega * r.y));
            v_norm_max = std::max(v_norm_max, tmp);
        }
        return v_norm_max;
    }

    void DynamicRigid::step_solid_body(double dt) {
        update_total_fluid_force_and_moment();
        RK4_step(dt);
        update_boundary();
    }

    void DynamicRigid::update_total_fluid_force_and_moment() {
        F_fluid = {0, 0};
        tau_fluid = 0;
        Point r{};
        for (int i{0}; i < n_bound; i++) {
            F_fluid += F_boundary[i];
            r = {boundary[i].x - y[0], boundary[i].y - y[1]};
            tau_fluid += r.cross(F_boundary[i]);
        }
        //cout << "F_fluid: "<<F_fluid<<", tau_fluid "<<tau_fluid<<endl;
        //F_fluid={-1000000,0};
        //tau_fluid=0;
    }


    Vector6d DynamicRigid::evaluate_f(Vector6d y_in) {
        //rhs of the state vector derivative dy/dt = f = [u_CM, v_CM, Fx/M, Fy/M, omega, tau/I]^T
        //std::pair<Point,double> F_S_tau_S{eval_solid_forces_and_moment()};
        //cout << "y= "<<y<<endl<<"M="<<M<<", I= "<<I<<endl;
        f[0] = y_in[2];
        f[1] = y_in[3];
        f[2] = (F_fluid.x + F_solid.x) / M;
        f[3] = (F_fluid.y + F_solid.y) / M;
        f[4] = y_in[5];
        f[5] = (tau_fluid + tau_solid) / I;
        //cout << "f="<<f<<endl;
        //f << -200,0,0,0,0,0;
        return f;
    }


    void DynamicRigid::RK4_step(double dt) {
        k1 = dt * evaluate_f(y);
        k2 = dt * evaluate_f(y + k1 / 2);
        k3 = dt * evaluate_f(y + k2 / 2);
        k4 = dt * evaluate_f(y + k3);
        y += (k1 + 2 * (k2 + k3) + k4) / 6;
    }

    void DynamicRigid::update_boundary() {
        double c{cos(y[4])};
        double s{sin(y[4])};
        for (int i{0}; i < n_bound; i++) {
            //boundary_i(t) = CM(t) + R(t)*r0
            boundary[i].x = y[0] + c * r0[i].x - s * r0[i].y;
            boundary[i].y = y[1] + s * r0[i].x + c * r0[i].y;
        }
    }

    DynamicRigid::~DynamicRigid() {
        delete[] r0;
    }


}
