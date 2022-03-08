//
// Created by anders on 1/30/22.
//

#ifndef FSIPROJECT_SOLID_UTILITIES_HPP
#define FSIPROJECT_SOLID_UTILITIES_HPP

#include "../../includes.hpp"
#include "../fluid/fvm_utilities.hpp"

namespace solid {

    struct Point {
        //also used as a vector
        double x, y;

        Point operator+(Point rhs) const { return {x + rhs.x, y + rhs.y}; }

        void operator+=(Point rhs) {
            x += rhs.x;
            y += rhs.y;
        }

        Point operator-(Point rhs) const { return {x - rhs.x, y - rhs.y}; }

        Point operator*(double rhs) const { return {x * rhs, y * rhs}; }

        double cross(Point rhs) const { return x * rhs.y - y * rhs.x; }

        double norm() const { return sqrt(x * x + y * y); }

        void normalize() {
            double n = norm();
            x /= n;
            y /= n;
        }

        double dot(Point rhs) const { return x * rhs.x + y * rhs.y; }

        friend std::ostream &operator<<(std::ostream &ost, const Point &rhs) {
            return ost << '(' << rhs.x << ',' << rhs.y << ')';
        }
    };

    struct Cell {
        inline static int nj; //To use the IX macro.
        int i, j;

        bool operator<(Cell rhs) const { return IX(i, j) < IX(rhs.i, rhs.j); } // To use the key lookup in std::map

        friend std::ostream& operator<<(std::ostream& ost, const Cell& rhs){
            return ost << '{' << rhs.i << ',' << rhs.j << '}';
        }
    };

    bool intersection(Point p1, Point p2, Point q1, Point q2);

    int orientation(Point p1, Point p2, Point p3);

    bool on_segment(Point p, Point q, Point r);

}


#endif //FSIPROJECT_SOLID_UTILITIES_HPP
