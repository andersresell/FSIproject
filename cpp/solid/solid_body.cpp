//
// Created by anders on 2/1/22.
//

#include "solid_body.hpp"

namespace solid {

    SolidBody::SolidBody(fluid::FVM_Solver &fvm, std::vector<Point>&& boundary_in, SolidBodyType type)
            : fvm{fvm}, n_bound{static_cast<unsigned int>(boundary_in.size())}, ni{fvm.ni}, nj{fvm.nj},
              dx{fvm.L_x / ni}, dy{fvm.L_y / nj}, type{type}{
        boundary = new Point[n_bound];
        for (int i{0}; i < boundary_in.size(); i++) {
            boundary[i] = boundary_in[i];
        }
        Cell::nj = nj;
    }

    void SolidBody::find_solid_cells() {
        Point p;
        for (int i{2}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                p = ind2point(i, j);
                if (point_inside(p)) {
                    solid_cells.push_back(Cell{i,j});
                    fvm.cell_status[IX(i, j)] = fluid::CellStatus::Solid;
                }
                else if (~fvm.is_static[IX(i,j)] && fvm.solids_initialized){
                    //If the cell is not a static solid cell the cell is set back to fluid in case it has moved
                    fvm.cell_status[IX(i, j)] = fluid::CellStatus::Fluid;
                }
            }
        }
    }

    void SolidBody::flag_static(){
        for (auto& c : solid_cells){
            fvm.is_static[IX(c.i,c.j)] = true;
        }
    }

    void SolidBody::find_ghost_cells(){
        //Checks if predetermined solid cells are ghost cells, by checking if nearby cells are fluid cells
        bool ghost;
        int i, j;
        for (Cell& c : solid_cells) {
            i = c.i;
            j = c.j;
            ghost = false;
            for (int ii{i - 2}; ii <= i + 2; ii++) {
                if (fvm.cell_status[IX(ii,j)] == fluid::CellStatus::Fluid){
                    ghost = true;
                    break;
                }
            }
            for (int jj{j - 2}; jj <= j + 2; jj++) {
                if (fvm.cell_status[IX(i,jj)] == fluid::CellStatus::Fluid) {
                    ghost = true;
                    break;
                }
            }
            if (ghost){
                fvm.cell_status[IX(i,j)] = fluid::CellStatus::Ghost;
                intercepts[c]; //assigning the key so that the ghost cells are known
            }
        }
    }

    void SolidBody::find_intercepts() {
        using namespace std;
        Point r1;
        Point r2;
        Point n;
        Point GP;
        Point q;
        Point p;
        Point d;
        unsigned int lp;
        double curr_dist;
        double prev_dist;
        for (auto& e : intercepts){
            GP = ind2point(e.first.i, e.first.j);
            prev_dist = INF;
            for (unsigned int l{0}; l < n_bound; l++) {
                lp = (l + 1) % n_bound;
                p = boundary[l];
                q = boundary[lp];
                r1 = q - GP;
                r2 = p - GP;
                n = {q.y - p.y, -(q.x - p.x)};
                n.normalize();
                // cout << "c1 "<<n.cross(r1)<<", c2 "<< r2.cross(n) <<endl;
                if (n.cross(r1) >= 0 && r2.cross(n) >= 0) { //intersection is on the line segment
                    d = n * r1.dot(n);
                    curr_dist = d.norm();
                    //cout << "curr_dist = "<<curr_dist<<endl;
                    if (curr_dist < prev_dist) {
                        e.second.first = GP + d;
                        e.second.second = n;
                        prev_dist = curr_dist;
                    }

                }
            }
        }
    }

    void SolidBody::interpolate_invicid_wall(fluid::vec4* U_in) {
        using namespace Eigen;
        using namespace std;
        using CS = fluid::CellStatus;
        //Constructing Vandermonde matrix VM for the case where all interpolation points are fluid points
        Matrix4d VM;
        VM << 1, 0, 0, 0,
                1, dx, 0, 0,
                1, dx, dy, dx * dy,
                1, 0, dy, 0;
        Matrix4d VM_inv_T = VM.transpose().inverse();
        Vector4d alpha;
        Vector4d alpha_dirichlet;
        Vector4d alpha_neumann;
        Vector4d VM_IP;
        Matrix4d VM_dirichlet;
        Matrix4d VM_neumann;
        int i, j;
        double DX, DY;
        Cell bottom_left_ind;
        Point bottom_left_point;
        fluid::CellStatus S1, S2, S3, S4;
        fluid::vec4 V_IP, V_GP, V1, V2, V3, V4;
        double u_n, u_t, u_n_GP, u_t_GP; //normal and tangential velocity components
        Point BI_neighbor;
        Point n;
        std::pair<Point, Point> BI_normal_pair;
        for (auto &e: intercepts) {
            //cout << "i = "<<e.first.i << "j = "<<e.first.j<<endl;
            //finding the image point, by using that IP = GP + 2*(BI - GP) = 2*BI - GP
            Point IP = e.second.first * 2 - ind2point(e.first.i, e.first.j);
            bottom_left_ind = point2ind(IP.x, IP.y);
            i = bottom_left_ind.i;
            j = bottom_left_ind.j;
            bottom_left_point = ind2point(i, j);
            DX = IP.x - bottom_left_point.x;
            DY = IP.y - bottom_left_point.y;
            VM_IP << 1, DX, DY, DX * DY;

            //May probably be optimized by extracting the variables from fvm.V
            V1 = fluid::FVM_Solver::conserved2primitive(U_in[IX(i, j)]);
            V2 = fluid::FVM_Solver::conserved2primitive(U_in[IX(i + 1, j)]);
            V3 = fluid::FVM_Solver::conserved2primitive(U_in[IX(i + 1, j + 1)]);
            V4 = fluid::FVM_Solver::conserved2primitive(U_in[IX(i, j + 1)]);

            S1 = fvm.cell_status[IX(i, j)];
            S2 = fvm.cell_status[IX(i + 1, j)];
            S3 = fvm.cell_status[IX(i + 1, j + 1)];
            S4 = fvm.cell_status[IX(i, j + 1)];
            if (S1 == CS::Fluid && S2 == CS::Fluid && S3 == CS::Fluid && S4 == CS::Fluid) {
                //All interpolation points are fluid cells
                alpha = VM_inv_T * VM_IP;
                V_IP = alpha(0) * V1 + alpha(1) * V2 + alpha(2) * V3 + alpha(3) * V4;
                n = e.second.second;
                //Dirichlet: u_n_BI = 0 -> u_n_GP = -u_n_IP
                u_n_GP = -(V_IP.u2 * n.x + V_IP.u3 * n.y);
                //Neumann: du_t/dn = 0, d_rho/dn = 0 and dp/dn = 0: u_t_GP = u_t_IP, rho_GP = rho_IP, p_GP = p_IP
                u_t_GP = -V_IP.u2 * n.y + V_IP.u3 * n.x;
                V_GP.u1 = V_IP.u1;
                V_GP.u4 = V_IP.u4;

            } else {
                if (S1 == CS::Solid || S2 == CS::Solid || S3 == CS::Solid || S4 == CS::Solid) {
                    std::cerr << "Solid point in bilinear interpolation stencil detected\n";
                }
                //Ghost cells are part of the interpolation stencil. Dirichlet and Neumann onditions at body interface
                //is used in the interpolation
                VM_dirichlet = VM;
                VM_neumann = VM;
                if (S1 == CS::Ghost) {
                    BI_normal_pair = intercepts[Cell{i, j}];
                    BI_neighbor = BI_normal_pair.first;
                    n = BI_normal_pair.second;
                    DX = BI_neighbor.x - bottom_left_point.x;
                    DY = BI_neighbor.y - bottom_left_point.y;
                    VM_dirichlet.block<1, 4>(0, 0) << 1, DX, DY, DX * DY;
                    VM_neumann.block<1, 4>(0, 0) << 0, n.x, n.y, DY * n.x + DX * n.y;
                }
                if (S2 == CS::Ghost) {
                    BI_normal_pair = intercepts[Cell{i + 1, j}];
                    BI_neighbor = BI_normal_pair.first;
                    n = BI_normal_pair.second;
                    DX = BI_neighbor.x - bottom_left_point.x;
                    DY = BI_neighbor.y - bottom_left_point.y;
                    VM_dirichlet.block<1, 4>(1, 0) << 1, DX, DY, DX * DY;
                    VM_neumann.block<1, 4>(1, 0) << 0, n.x, n.y, DY * n.x + DX * n.y;
                }
                if (S3 == CS::Ghost) {
                    BI_normal_pair = intercepts[Cell{i + 1, j + 1}];
                    BI_neighbor = BI_normal_pair.first;
                    n = BI_normal_pair.second;
                    DX = BI_neighbor.x - bottom_left_point.x;
                    DY = BI_neighbor.y - bottom_left_point.y;
                    VM_dirichlet.block<1, 4>(2, 0) << 1, DX, DY, DX * DY;
                    VM_neumann.block<1, 4>(2, 0) << 0, n.x, n.y, DY * n.x + DX * n.y;
                }
                if (S4 == CS::Ghost) {
                    BI_normal_pair = intercepts[Cell{i, j + 1}];
                    BI_neighbor = BI_normal_pair.first;
                    n = BI_normal_pair.second;
                    DX = BI_neighbor.x - bottom_left_point.x;
                    DY = BI_neighbor.y - bottom_left_point.y;
                    VM_dirichlet.block<1, 4>(3, 0) << 1, DX, DY, DX * DY;
                    VM_neumann.block<1, 4>(3, 0) << 0, n.x, n.y, DY * n.x + DX * n.y;
                }
                VM_dirichlet.transposeInPlace();
                VM_neumann.transposeInPlace();
                alpha_dirichlet = VM_dirichlet.partialPivLu().solve(VM_IP);
                alpha_neumann = VM_neumann.partialPivLu().solve(VM_IP);

                u_n_GP = 0;
                u_t_GP = 0;
                V_GP *= 0;
                n = e.second.second;
                if (S1 == CS::Fluid) {
                    u_n_GP -= alpha_dirichlet(0) * (V1.u2 * n.x + V1.u3 * n.y);
                    u_t_GP += alpha_neumann(0) * (-V1.u2 * n.y + V1.u3 * n.x);
                    V_GP.u1 += alpha_neumann(0) * V1.u1;
                    V_GP.u4 += alpha_neumann(0) * V1.u4;
                }
                if (S2 == CS::Fluid) {
                    u_n_GP -= alpha_dirichlet(1) * (V2.u2 * n.x + V2.u3 * n.y);
                    u_t_GP += alpha_neumann(1) * (-V2.u2 * n.y + V2.u3 * n.x);
                    V_GP.u1 += alpha_neumann(1) * V2.u1;
                    V_GP.u4 += alpha_neumann(1) * V2.u4;
                }
                if (S3 == CS::Fluid) {
                    u_n_GP -= alpha_dirichlet(2) * (V3.u2 * n.x + V3.u3 * n.y);
                    u_t_GP += alpha_neumann(2) * (-V3.u2 * n.y + V3.u3 * n.x);
                    V_GP.u1 += alpha_neumann(2) * V3.u1;
                    V_GP.u4 += alpha_neumann(2) * V3.u4;
                }
                if (S4 == CS::Fluid) {
                    u_n_GP -= alpha_dirichlet(3) * (V4.u2 * n.x + V4.u3 * n.y);
                    u_t_GP += alpha_neumann(3) * (-V4.u2 * n.y + V4.u3 * n.x);
                    V_GP.u1 += alpha_neumann(3) * V4.u1;
                    V_GP.u4 += alpha_neumann(3) * V4.u4;
                }
            }
            //transform back
            V_GP.u2 = n.x * u_n_GP - n.y * u_t_GP;
            V_GP.u3 = n.y * u_n_GP + n.x * u_t_GP;

            U_in[IX(e.first.i, e.first.j)] = fluid::FVM_Solver::primitive2conserved(V_GP);

        }

    }

    bool SolidBody::point_inside(Point p) const{
        //using namespace std;
        //See https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-Polygon/ for details
        Point extreme{INF, p.y};
        unsigned int count{0};
        for (unsigned int i{0}; i < n_bound; i++) {
            unsigned int next = (i + 1) % n_bound;
            if (intersection(p, extreme, boundary[i], boundary[next])) {
                //Special case
                if (orientation(boundary[i], p, boundary[next]) == 0)
                    //I don't fully understand this one. Why should it return false if it's not on the segment?
                    return on_segment(boundary[i], p, boundary[next]);
                count++;
            }
        }
        return count % 2 == 1; //If number of intersections are odd, return true, else return false
    }

    SolidBody::~SolidBody() {
        delete[] boundary;
    }
}
