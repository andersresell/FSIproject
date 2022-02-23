//
// Created by anders on 2/1/22.
//

#ifndef FSIPROJECT_SolidBody_HPP
#define FSIPROJECT_SolidBody_HPP

#include "../../includes.hpp"
#include "solid_utilities.hpp"
#include "../fluid/fvm_solver.hpp"



namespace solid {

    enum class SolidBodyType{Static, Movable};

    class SolidBody {
        std::vector<Cell> solid_cells;
        //A map from ghost cells to intercepts. Key is the ghost cell, 1st value is the intercept, 2nd value is the normal
        std::map<Cell, std::pair<Point, Point>> intercepts;
        inline const static double INF = 1e6; //A large number used as inf in the calculation of points inside Boundary.
        fluid::FVM_Solver &fvm;
        int ni, nj;
        double dx, dy;

    public:
        Point *boundary; //A polygon defining the boundary
        const SolidBodyType type;
        const unsigned int n_bound;

        SolidBody(fluid::FVM_Solver &fvm, const std::vector<Point> &boundary_in, SolidBodyType type);

        //void set_bc(fluid::vec4* U_in);

        void debug_csv();

        void debug_intercepts_csv();

        //void write_boundary_csv(const std::string &output_folder, int solid_index, int n=-1);

        void find_solid_cells();

        void flag_static();

        void find_ghost_cells();

        //Point compute_intercept(Point GP, Point p, Point q);

        void find_intercepts();

        void interpolate_invicid_wall(fluid::vec4 *U_in);

    private:
        bool point_inside(Point p) const; //Check wether a poins is inside the solid boundary

        Point ind2point(int i, int j) const { return {(i - 1.5) * dx, (j - 1.5) * dy}; }

        Cell point2ind(double x, double y) const { //rounding a point down to nerest cell index
            return {(int) (x / dx + 1.5), (int) (y / dy + 1.5)};
        }

    public:
        ~SolidBody();
    };





    class DeformableSolidBody : SolidBody {
        
    };


}

#endif //FSIPROJECT_SolidBody_HPP
