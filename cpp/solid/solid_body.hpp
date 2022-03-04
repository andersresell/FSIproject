//
// Created by anders on 2/1/22.
//

#ifndef FSIPROJECT_SolidBody_HPP
#define FSIPROJECT_SolidBody_HPP

#include "../../includes.hpp"
#include "solid_utilities.hpp"
#include "../fluid/fvm_solver.hpp"



namespace solid {

    enum class SolidBodyType{Static, Dynamic};

    struct Segment{
        //Holds the ghost points, intercepts, normal vector, etc for each segment

        std::vector<Point> intercepts;
        Point n;
        bool n_is_set;
        //void sort_intercepts(int p, int q); //to be used for pressure integration
        Segment() : n_is_set{false}{}
    };


    class SolidBody {
        //set and vector have their own pros and cons here, trying to use both for different tasks for performance
        std::set<Cell> solid_cells;
        bool first_timestep;
        //A map from ghost cells to intercepts. Key is the ghost cell, 1st value is the intercept, 2nd value is the normal
        inline const static double INF = 1e6; //A large number used as inf in the calculation of points inside Boundary.
        fluid::FVM_Solver &fvm;
        int ni, nj;
        double dx, dy;

        //Variables used for interpolation
        Eigen::Matrix4d A;
        Eigen::Matrix4d A_inv_T;
    public:
        //std::map<Cell, GP_info> intercepts; //[Cell GP, {Point BI, Point n}]
        Segment *segments;
        std::map<Cell, std::pair<Point, Point>> cell2intercept; //key = Cell GP, value = {Point BI, Point n}
        Point *boundary; //A polygon defining the boundary
        Point *F_boundary; //Lumped force on each boundary node at current timestep;
        const int n_bound;
        const SolidBodyType type;


        SolidBody(fluid::FVM_Solver &fvm, std::vector<Point> &&boundary_in, SolidBodyType type);

        void find_solid_cells();

        void flag_static();

        void find_ghost_cells();

        void find_intercepts();

        double segment_length(int p) const {return (boundary[(p+1)%n_bound] - boundary[p]).norm();}

        void interpolate_fresh_points(fluid::vec4* U_in);

        void interpolate_solid(fluid::vec4 *U_in);

        void interpolate_cell(Cell GP, Point BI, Point n, fluid::vec4* U_in);

        double interpolate_dirichlet(const Eigen::Vector4d& alpha_dir, std::vector<double>&& phi,
                                     std::vector<double>&& phi_BI ,const std::vector<fluid::CellStatus>& cs);

        double interpolate_neumann(const Eigen::Vector4d& alpha_neu, std::vector<double>&& phi,
                            std::vector<double>&& phi_BI_derivative ,
                            const std::vector<fluid::CellStatus>& cs, double Delta_l);
        double interpolate_dirichlet_zero_value(const Eigen::Vector4d& alpha_dir, std::vector<double>&& phi,
                                                const std::vector<fluid::CellStatus>& cs);

        double interpolate_neumann_zero_gradient(const Eigen::Vector4d& alpha_neu, std::vector<double>&& phi,
                                   const std::vector<fluid::CellStatus>& cs);
        virtual void bundary_vel_and_acc(Point BI, Point& v_wall, Point& a_wall) const;

        bool cell_within_grid(int i, int j){ return i >=2  && i<ni+2 && j>=2 && j< nj+2;}

        void reset_containers();

        update_lumped_forces(fluid::vec4 *U_in);//Integrates the pressure over the surface and returns {F_x,F_y} and tau

    private:
        bool point_inside(Point p) const; //Check wether a poins is inside the solid boundary

        Point ind2point(int i, int j) const { return {(i - 1.5) * dx, (j - 1.5) * dy}; }

        Cell point2ind(double x, double y) const { //rounding a point down to nerest cell index
            return {(int) (x / dx + 1.5), (int) (y / dy + 1.5)};
        }

    public:
        ~SolidBody();
    };

    typedef Eigen::Matrix<double,1,6> Vector6d;

    class DynamicRigid : public SolidBody {
        double M;
        double I;
        Vector6d y; //state vector of the rigid body: y = [x_CM, y_CM, u_CM, v_CM, theta, omega]^T
        Vector6d k1,k2,k3,k4;
        //const static inline int n_state{6};
        Point F_fluid; //Total force from the fluid. Only updated once per timestep
        Point F_solid;
        double tau_fluid; //Total moment from the fluid. Only updated once per timestep
        double tau_solid;
    public:
        DynamicRigid(fluid::FVM_Solver &fvm, std::vector<Point>&& boundary_in, Point CM, double M, double I);

    private:
        std::pair<Point,double> eval_solid_forces_and_moment();
        Vector6d evaluate_f(Vector6d y_in);
        Vector6d RK4_step(double dt);

        void bundary_vel_and_acc(Point BI, Point& v_wall, Point& a_wall) const final;
    };
    void DynamicRigid::bundary_vel_and_acc(Point BI, Point& v_wall, Point& a_wall) const{
        Point r = {BI.x - y[0], BI.y - y[1]};
        double omega = y[5];
        //v = v_CM + omega x r
        v_wall = {y[2] - omega*r.y, y[3] + omega*r.y};
        Point a_CM = (F_fluid + F_solid)*(1/M);
        double alpha = (tau_fluid + tay_solid)*(1/I);
        //a = a_CM + omega x (omega x r) + alpha x r
        a_wall = {a_CM.x - omega*omega*r.x - alpha*r.y, a_CM.y - omega*omega*r.y + alpha*r.x};
    }
}

#endif //FSIPROJECT_SolidBody_HPP
