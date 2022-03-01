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

    class SolidBody {
        std::vector<Cell> solid_cells;
        //A map from ghost cells to intercepts. Key is the ghost cell, 1st value is the intercept, 2nd value is the normal
        inline const static double INF = 1e6; //A large number used as inf in the calculation of points inside Boundary.
        fluid::FVM_Solver &fvm;
        int ni, nj;
        double dx, dy;

    public:
        std::map<Cell, std::pair<Point, Point>> intercepts;
        Point *boundary; //A polygon defining the boundary
        const unsigned int n_bound;
        const SolidBodyType type;


        SolidBody(fluid::FVM_Solver &fvm, std::vector<Point>&& boundary_in, SolidBodyType type);

        void find_solid_cells();

        void flag_static();

        void find_ghost_cells();

        void find_intercepts();

        void interpolate_invicid_wall(fluid::vec4 *U_in);

        Point integrate_pressure(fluid::vec4* U_in);//Integrates the pressure over the surface and returns {F_x,F_y}

    private:
        bool point_inside(Point p) const; //Check wether a poins is inside the solid boundary

        Point ind2point(int i, int j) const { return {(i - 1.5) * dx, (j - 1.5) * dy}; }

        Cell point2ind(double x, double y) const { //rounding a point down to nerest cell index
            return {(int) (x / dx + 1.5), (int) (y / dy + 1.5)};
        }

    public:
        ~SolidBody();
    };





    class DynamicRigid : public SolidBody {
        
    };


}

#endif //FSIPROJECT_SolidBody_HPP
