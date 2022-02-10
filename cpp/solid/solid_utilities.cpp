//
// Created by anders on 1/30/22.
//

#include "solid_utilities.hpp"

namespace solid {

    bool intersection(Point p1, Point q1, Point p2, Point q2) {
        //Checks wether the two lines segments p1-q1 and p2-q2 intersect
        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);
        if (o1 != o2 && o3 != o4) {
            return true;
        } else if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0) { //Special case
            if (on_segment(p1, p2, q1) || on_segment(p1, q2, q1)) return true;
        }
        return false;
    }

    int orientation(Point p1, Point p2, Point p3) {
        //Checks wether three points are oriented as CCW, CW or colinear by checking the value of the
        //cross product. returns 1 for CCW 0 for colinear and -1 for CW
        double cross_prod = (p2.x - p1.x) * (p3.y - p2.y) - (p3.x - p2.x) * (p2.y - p1.y);
        if (cross_prod == 0) return 0;
        return cross_prod > 0 ? 1 : -1;
    };

    bool on_segment(Point p, Point q, Point r) {
        //Given three colinear points, the function checks wether point q is on segment p-r
        if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
            q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y)) {
            return true;
        }
        return false;
    }

}