//
// Created by anders on 2/1/22.
//

#ifndef FSIPROJECT_SolidBody_HPP
#define FSIPROJECT_SolidBody_HPP

#include "../../includes.hpp"
#include "solid_utilities.hpp"
#include "../fluid/fvm_solver.hpp"



namespace solid {

    enum class SolidBodyType {
        Static, Dynamic
    };

    struct GP_data {
        Point BI;
        Point n;
        //double p;
    };

    class SolidBody {
    protected:
        //set and vector have their own pros and cons here, trying to use both for different tasks for performance
        std::set<Cell> solid_cells;
        //std::vector<Cell> solid_cells;
        bool first_timestep;
        //A map from ghost cells to intercepts. Key is the ghost cell, 1st value is the intercept, 2nd value is the normal
        inline const static double INF = 1e6; //A large number used as spatial infinity
        fluid::FVM_Solver &fvm;
        int ni, nj;
        double dx, dy;
        int timestep;

        //Variables used for interpolation
        Eigen::Matrix4d A;
        Eigen::Matrix4d A_inv_T;
    public:
        std::map<Cell, GP_data> cell2intercept; //key = Cell GP, value = {Point BI, Point n, double p}
        std::map<Cell, GP_data> cell2intercept_FP;
        Point *boundary; //A polygon defining the boundary
        //polygon_t boundary;
        Point *F_boundary; //Lumped force on each boundary node at current timestep;
        const int n_bound;
        const SolidBodyType type;


        SolidBody(fluid::FVM_Solver &fvm, std::vector<Point> &&boundary_in, SolidBodyType type);

        void step(fluid::vec4 *U_in, double dt, bool update_solid_pos);

        void find_solid_cells();

        void flag_static();

        void find_ghost_cells();

        void find_intercepts(std::map<Cell, GP_data> &intercept_map);

        double segment_length(int p) const { return (boundary[(p + 1) % n_bound] - boundary[p]).norm();}
            //const point_t& p1 = boundary.outer()[(p + 1) % n_bound];
            //const point_t& p2 = boundary.outer()[p];
            //return sqrt(sqr(p2.x() - p1.x()) + sqr(p2.y() - p1.y()));

        void interpolate_fresh_points(fluid::vec4 *U_in);

        void interpolate_ghost_points(fluid::vec4 *U_in);

        void interpolate_cell(Cell point,const std::map<Cell, GP_data> &intercept_map, fluid::vec4 *U_in,
                              bool fresh_point = false);

        double interpolate_dirichlet(const Eigen::Vector4d &alpha_dir, std::array<double, 4> phi,
                                     std::array<fluid::CellStatus, 4> cs,
                                     std::array<double, 4> phi_BI_adj = {0, 0, 0, 0});

        double interpolate_neumann(const Eigen::Vector4d &alpha_neu, std::array<double, 4> phi,
                                   std::array<fluid::CellStatus, 4> cs,
                                   std::array<double, 4> phi_derivative_BI_adj = {0, 0, 0, 0});
        /*double interpolate_dirichlet_zero_value(const Eigen::Vector4d& alpha_dir, std::vector<double>&& phi,
                                                const std::vector<fluid::CellStatus>& cs);

        double interpolate_neumann_zero_gradient(const Eigen::Vector4d& alpha_neu, std::vector<double>&& phi,
                                   const std::vector<fluid::CellStatus>& cs);*/

        //The virtual functions needed definitions to avoid vtable error, even though they are not used in the base class
        virtual void boundary_vel_and_acc(Point BI, Point &v_wall, Point &a_wall) const {};

        virtual double max_boundary_speed() const { return 0; }

        virtual void step_solid_body(double dt) {};

        bool cell_within_grid(int i, int j) { return i >= 2 && i < ni + 2 && j >= 2 && j < nj + 2; }

        bool point_within_grid(Point p) {
            return p.x >= 0.5 * dx && p.x <= (ni + 0.5) * dx && p.y >= 0.5 * dy && p.y <= (nj + 0.5) * dy;
        }

        void reset_containers();

        void update_lumped_forces(
                fluid::vec4 *U_in);//Integrates the pressure over each segment and lumps the resultant force in each node

        void write_fresh_points(int n) const;

    private:
        bool point_inside(Point p) const; //Check wether a poins is inside the solid boundary

        Point ind2point(int i, int j) const { return {(i - 1.5) * dx, (j - 1.5) * dy}; }

        Cell point2ind(double x, double y) const { //rounding a point down to nerest cell index
            return {(int) (x / dx + 1.5), (int) (y / dy + 1.5)};
        }

    public:
        virtual ~SolidBody();
    };

}
#endif //FSIPROJECT_SolidBody_HPP
