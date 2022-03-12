//
// Created by anders on 3/11/22.
//

#ifndef FSIPROJECT_DYNAMIC_RIGID_HPP
#define FSIPROJECT_DYNAMIC_RIGID_HPP

#include "solid_body.hpp"

namespace solid {

    typedef Eigen::Matrix<double, 1, 6> Vector6d;

    class DynamicRigid : public SolidBody {
        double M;
        double I;
        Vector6d y; //state vector of the rigid body: y = [x_CM, y_CM, u_CM, v_CM, theta, omega]^T
        Vector6d k1, k2, k3, k4;
        Vector6d f;
        Point *r0; // The radius from CM to each boundary node at t=0
        //const static inline int n_state{6};
        Point F_fluid; //Total force from the fluid. Only updated once per timestep
        Point F_solid;
        double tau_fluid; //Total moment from the fluid. Only updated once per timestep
        double tau_solid;
    public:
        DynamicRigid(fluid::FVM_Solver &fvm, std::vector<Point> &&boundary_in, Point CM, double M, double I);

        void step_solid_body(double dt) final;

        double max_boundary_speed() const final;

    private:
        void update_total_fluid_force_and_moment();

        Vector6d evaluate_f(Vector6d y_in);

        void RK4_step(double dt);

        void update_boundary();

        void bundary_vel_and_acc(Point BI, Point &v_wall, Point &a_wall) const final;

    public:
        virtual ~DynamicRigid();
    };


    inline void DynamicRigid::bundary_vel_and_acc(Point BI, Point &v_wall, Point &a_wall) const {
        using namespace std;
        Point r = {BI.x - y[0], BI.y - y[1]};
        double omega = y[5];
        //v = v_CM + omega x r
        v_wall = {y[2] - omega * r.y, y[3] + omega * r.y};
        Point a_CM = (F_fluid + F_solid) * (1 / M);
        double alpha = (tau_fluid + tau_solid) * (1 / I);
        //a = a_CM + omega x (omega x r) + alpha x r
        a_wall = {a_CM.x - omega * omega * r.x - alpha * r.y, a_CM.y - omega * omega * r.y + alpha * r.x};
        a_wall = {0, 0};
        //cout << "y = ";
        for (int i{0}; i < 6; i++) {
            //cout << y[i] <<",";
        }//cout << endl;
    }

}

#endif //FSIPROJECT_DYNAMIC_RIGID_HPP
