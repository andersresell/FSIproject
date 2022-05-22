//
// Created by anders on 3/11/22.
//

#ifndef FSIPROJECT_DYNAMIC_RIGID_HPP
#define FSIPROJECT_DYNAMIC_RIGID_HPP

#include "solid_body.hpp"

namespace solid {

    typedef Eigen::Matrix<double, 1, 6> Vector6d;

    enum class MaterialModel {None, ViscoElastic};

    struct RigidConstraints{
        bool constrained; //removes the y and theta dofs
        std::pair<bool,double> prescribed_velocity;
        MaterialModel material_model;
        double K; //spring stiffnes
        double C; //Viscous damping
        RigidConstraints();
        void setup_prescribed_velocity(double velocity_x);
        void setup_viscoelastic(double K, double C);
        Point calc_F_solid(double deflection, double velocity) const;
    };

    class DynamicRigid : public SolidBody {
        double M;
        double I;
        Vector6d y; //state vector of the rigid body: y = [x_CM, y_CM, u_CM, v_CM, theta, omega]^T
        Vector6d k1, k2, k3, k4;
        Vector6d f;
        Point *r0; // The radius from CM to each boundary node at t=0
        const Point CM0; //CM at t=0
        Point F_fluid; //Total force from the fluid. Only updated once per timestep
        Point F_solid;
        double tau_fluid; //Total moment from the fluid. Only updated once per timestep
        double tau_solid;
        double omega_prev;
        Point vel_CM_prev;
        double dt_prev;
        RigidConstraints rigid_constraints;
    public:

        void add_rigid_constraints(RigidConstraints&& r_c){
            using namespace std;
            rigid_constraints = r_c;
            y[2] = rigid_constraints.prescribed_velocity.first ? rigid_constraints.prescribed_velocity.second : 0;
        }

        DynamicRigid(fluid::FVM_Solver &fvm, std::vector<Point> &&boundary_in, Point CM, double M, double I);

        void step_solid_body(double dt) final;

        double max_boundary_speed() const final;

        Point get_CM_velocity() const final {return {y[2],y[3]};}

    private:
        void update_total_fluid_force_and_moment();

        Vector6d evaluate_f(Vector6d y_in);

        void RK4_step(double dt);

        void update_boundary();

        void boundary_vel_and_fluid_boundary_acc(Point BI, Point n, Point u_fluid, Point &u_wall, double &a_normal) const final;

    public:
        virtual ~DynamicRigid();
    };


    inline void DynamicRigid::boundary_vel_and_fluid_boundary_acc(Point BI, Point n, Point u_fluid, Point &u_wall, double &a_normal) const {
        //This function calculates the velocity of the solid at the boundary, and estimates the fluid acceleration at
        //the boundary
        using namespace std;
        Point r = {BI.x - y[0], BI.y - y[1]};
        double omega = y[5];
        //v = v_CM + omega x r
        u_wall = {y[2] - omega * r.y, y[3] + omega * r.y};
        //Point a_CM = (F_fluid + F_solid) * (1 / M); //Won't work for prescribed velocities
        Point a_CM = {(y[2] - vel_CM_prev.x)/dt_prev, (y[3] - vel_CM_prev.y)/dt_prev}; //A backward difference to approximate acceleration
        //double alpha = (tau_fluid + tau_solid) * (1 / I);
        double alpha = (y[5] - omega_prev)/dt_prev; //A backward difference to approximate the angular acceleration.
        //a = a_CM + omega x (omega x r) + alpha x r
        //a_wall = {a_CM.x - omega * omega * r.x - alpha * r.y, a_CM.y - omega * omega * r.y + alpha * r.x};

        //Improving the estimate of the fluid vel at the wall by using the correct value for the normal velocity component
        // taken from the solid wall velocity. The tangential component is still estimated.
        double u_t = -n.y*u_fluid.x + n.x*u_fluid.y;
        double u_n = u_wall.x*n.x + u_wall.y*n.y;
        u_fluid.x = n.x*u_n - n.y*u_t;
        u_fluid.y = n.y*u_n + n.x*u_t;

        //a_wall dot n = (a_CM + alpha x r + omega x (omega x r) + 2 omega x (u_fluid - u_CM - omega x r)) dot n
        a_normal = (a_CM.x - alpha*r.y - omega*omega*r.x + 2*omega*(-u_fluid.y + y[3] + omega*r.y))*n.x
                + (a_CM.y + alpha*r.x - omega*omega*r.y + 2*omega*(u_fluid.x - y[2] - omega*r.x))*n.y;

        //cout << "u "<<y[2]<<", u_old" << vel_CM_prev.x<<", a_CM "<<a_CM<<", alpha "<<alpha<<",omega "<<omega<< ", a_wall" <<a_wall << endl;
        //a_wall = {0, 0};
    }

}

#endif //FSIPROJECT_DYNAMIC_RIGID_HPP
