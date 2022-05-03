//
// Created by anders on 3/11/22.
//

#include "dynamic_rigid.hpp"

using namespace std;

namespace solid {

    RigidConstraints::RigidConstraints()
    : constrained{false}, prescribed_velocity{false,0}, material_model{MaterialModel::None}, K{0}, C{0}
    {
    }

    void RigidConstraints::setup_prescribed_velocity(double velocity_x){
        constrained = true;
        prescribed_velocity = {true, velocity_x};
        material_model = MaterialModel::None;
    }

    void RigidConstraints::setup_viscoelastic(double spring_stiffness, double damping_coefficient){
        constrained = true;
        prescribed_velocity.first = false;
        material_model = MaterialModel::ViscoElastic;
        K = spring_stiffness;
        C = damping_coefficient;
    }


    Point RigidConstraints::calc_F_solid(double deflection, double velocity) const{
        if (material_model == MaterialModel::None){
            return {0,0};
        }else if(material_model == MaterialModel::ViscoElastic){
            return {-K*deflection -C*velocity, 0};
        }else{
            std::cerr << "Illegal material model in solid force calculation\n"; exit(1);
        }
    }

    DynamicRigid::DynamicRigid(fluid::FVM_Solver &fvm, std::vector<Point> &&boundary_in, Point CM, double M, double I)
            : SolidBody(fvm, std::move(boundary_in), SolidBodyType::Dynamic), M{M}, I{I}, CM0{CM}, rigid_constraints{},
            F_fluid{0,0}, F_solid{0,0}, tau_solid{0}, omega_prev{0}, vel_CM_prev{0,0}, dt_prev{INF} {
        y[0] = CM.x;
        y[1] = CM.y;
        y[2] = 0;
        y[3] = 0;
        y[4] = 0;
        y[5] = 0;
        r0 = new Point[n_bound];
        for (int i{0}; i < n_bound; i++) {
            r0[i] = boundary[i] - CM;
            //r0[i] = Point{boundary.outer()[i].x(),boundary.outer()[i].y()} - CM;
        }
    }

    double DynamicRigid::max_boundary_speed() const {
        double tmp;
        double v_norm_max{0};
        double omega = y[5];
        for (int i{0}; i < n_bound; i++) {
            Point r = {boundary[i].x - y[0], boundary[i].y - y[1]};
            //v = v_CM + omega x r
            tmp = sqrt(sqr(y[2] - omega * r.y) + sqr(y[3] + omega * r.y));
            v_norm_max = std::max(v_norm_max, tmp);
        }
        /*for (int i{0}; i < n_bound; i++) {
            Point r = {boundary.outer()[i].x() - y[0], boundary.outer()[i].y() - y[1]};
            //v = v_CM + omega x r
            tmp = sqrt(sqr(y[2] - omega * r.y) + sqr(y[3] + omega * r.y));
            v_norm_max = std::max(v_norm_max, tmp);
        }*/
        return v_norm_max;
    }

    void DynamicRigid::step_solid_body(double dt) {
        if (rigid_constraints.prescribed_velocity.first){
            y[2] = rigid_constraints.prescribed_velocity.second;
        } else{
            update_total_fluid_force_and_moment();
        }
        if (rigid_constraints.constrained){
            F_fluid.y = 0;
            tau_fluid = 0;
        }
        cout << "F_fluid "<<F_fluid<<endl;
        omega_prev = y[5];
        vel_CM_prev = {y[2],y[3]};
        dt_prev = dt;
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
        /*for (int i{0}; i < n_bound; i++) {
            F_fluid += F_boundary[i];
            r = {boundary.outer()[i].x() - y[0], boundary.outer()[i].y() - y[1]};
            tau_fluid += r.cross(F_boundary[i]);
        }*/
    }

    Vector6d DynamicRigid::evaluate_f(Vector6d y_in) {
        //rhs of the state vector derivative dy/dt = f = [u_CM, v_CM, Fx/M, Fy/M, omega, tau/I]^T
        F_solid = rigid_constraints.calc_F_solid(y_in[0]-CM0.x,y_in[2]);
        f[0] = y_in[2];
        f[1] = y_in[3];
        f[2] = (F_fluid.x + F_solid.x) / M;
        f[3] = (F_fluid.y + F_solid.y) / M;
        f[4] = y_in[5];
        f[5] = (tau_fluid + tau_solid) / I;

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
            assert(!isnan(y[0]) && !isnan(y[1]));
            boundary[i].x = y[0] + c * r0[i].x - s * r0[i].y;
            boundary[i].y = y[1] + s * r0[i].x + c * r0[i].y;
            //boundary.outer()[i] = {y[0] + c * r0[i].x - s * r0[i].y, y[1] + s * r0[i].x + c * r0[i].y};
        }
    }

    DynamicRigid::~DynamicRigid() {
        delete[] r0;
    }


}
