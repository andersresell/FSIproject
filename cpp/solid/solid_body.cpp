//
// Created by anders on 2/1/22.
//

#include "solid_body.hpp"

namespace solid {

    using namespace Eigen;


    SolidBody::SolidBody(fluid::FVM_Solver &fvm, std::vector<Point>&& boundary_in, SolidBodyType type)
            : fvm{fvm}, n_bound{(int)boundary_in.size()}, ni{fvm.ni}, nj{fvm.nj},
              dx{fvm.L_x / ni}, dy{fvm.L_y / nj}, type{type}, first_timestep{true}{
        boundary = new Point[n_bound];
        F_boundary = new Point[n_bound];
        for (int i{0}; i < boundary_in.size(); i++) {
            boundary[i] = boundary_in[i];
        }
        Cell::nj = nj;

        A << 1, 0, 0, 0,
                1, dx, 0, 0,
                1, dx, dy, dx * dy,
                1, 0, dy, 0;
        A_inv_T = A.transpose().inverse();
    }


    void SolidBody::update(fluid::vec4* U_in){
        if (type == SolidBodyType::Static){
            if (first_timestep){
                find_solid_cells();
                flag_static();
                find_ghost_cells();
                find_intercepts();
            }
            interpolate_solid(U_in);
        }else if (type == SolidBodyType::Dynamic){
            if(~first_timestep){
                //Step solid(dt)
            }
            find_solid_cells();
            find_ghost_cells();
            find_intercepts();
            interpolate_fresh_points(U_in);
            interpolate_solid(U_in);
            //Clean up containers
        }
        first_timestep = false;
    }

    void SolidBody::find_solid_cells() {
        Point p{};
        if (first_timestep) {
            //During the first timestep all cells has to be detected
            for (int i{2}; i < ni + 2; i++) {
                for (int j{2}; j < nj + 2; j++) {
                    p = ind2point(i, j);
                    if (point_inside(p)) {
                        //solid_cells.push_back(Cell{i, j});
                        solid_cells.emplace(Cell{i,j});
                        fvm.cell_status[IX(i, j)] = fluid::CellStatus::Solid;
                    } /*else if (~fvm.is_static[IX(i, j)] && fvm.solids_initialized) {
                        //If the cell is not a static solid cell the cell is set back to fluid in case it has moved
                        fvm.cell_status[IX(i, j)] = fluid::CellStatus::Fluid;
                    }*/
                }
            }
        } else if (type == SolidBodyType::Dynamic){
            for (const auto & e: cell2intercept){
                //The cells in the old GP's are set to fluid
                fvm.cell_status[IX(e.first.i, e.first.j)] = fluid::CellStatus::Fluid;
            }
            //Exploiting that the CFL condition only allows the boundary to travel one cell to find the new solid cells
            int i,j;
            for (const auto & e: cell2intercept){
                for (int ii{-1};ii<=1;ii++){
                    for (int jj{-1};jj<=1;jj++){
                        i = e.first.i;
                        j = e.first.j;
                        //Continuing if the indices are outside the computational domain or if the cell is allready set
                        if (~cell_within_grid(i+ii,j+jj) || fvm.cell_status[IX(i+ii, j+jj)] == fluid::CellStatus::Solid)
                            continue;
                        p = ind2point(i+ii,j+jj);
                        if (point_inside(p)){
                            fvm.cell_status[IX(i+ii, j+jj)] = fluid::CellStatus::Solid;
                        }else{
                            solid_cells.erase({i+ii,j+jj});
                        }
                    }
                }
            }
        }else{
            std::cerr << "Attempting to find solid cells for a non dynamic body after the first timestep. "
                         "<<Unnecessary operation\n"; exit(1);
        }
    }

    void SolidBody::flag_static(){
        for (auto& c : solid_cells){
            fvm.is_static[IX(c.i,c.j)] = true;
        }
    }

    void SolidBody::find_ghost_cells() {
        //Checks if predetermined solid cells are ghost cells, by checking if nearby cells are fluid cells

        //Copying all elements to a vector since it loops faster.
        std::vector<Cell> solid_cells_copy;
        solid_cells_copy.reserve(solid_cells.size());
        std::copy(solid_cells.begin(), solid_cells.end(), solid_cells_copy.begin());

        bool ghost;
        int i, j;
        for (Cell c: solid_cells_copy) {
            i = c.i;
            j = c.j;
            ghost = false;
            for (int ii{i - 2}; ii <= i + 2; ii++) {
                if (fvm.cell_status[IX(ii, j)] == fluid::CellStatus::Fluid) {
                    ghost = true;
                    break;
                }
            }
            for (int jj{j - 2}; jj <= j + 2; jj++) {
                if (fvm.cell_status[IX(i, jj)] == fluid::CellStatus::Fluid) {
                    ghost = true;
                    break;
                }
            }
            if (ghost) {
                fvm.cell_status[IX(i, j)] = fluid::CellStatus::Ghost;
                cell2intercept[c];
            }
        }
    }
    /*
    void SolidBody::find_intercepts() {
        using namespace std;
        Point r1;
        Point r2;
        Point n;
        Point GP;
        Point q;
        Point p;
        Point d;
        int ip;
        double curr_dist;
        double prev_dist;
        for (auto& e : intercepts){
            GP = ind2point(e.first.i, e.first.j);
            prev_dist = INF;
            for (int i{0}; i < n_bound; i++) {
                ip = (i + 1) % n_bound;
                p = boundary[i];
                q = boundary[ip];
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
                        e.second.BI = GP + d;
                        e.second.n = n;
                        e.second.boundary_ind = i;
                        prev_dist = curr_dist;
                    }

                }
            }
        }
    }
*/
    void SolidBody::find_intercepts() {
        using namespace std;
        Point r1;
        Point r2;
        Point n;
        Point GP;
        Point q;
        Point p;
        Point d;
        int jp;
        double curr_dist;
        double prev_dist;
        for (auto& e : cell2intercept){
            GP = ind2point(e.first.i,e.first.j);
            prev_dist = INF;
            for (int j{0}; j < n_bound; j++) {
                //A different idea is to store the length from the previous iteration and use this as an estimate.
                segments[j].segment_data_vec.clear(); //Maybe unneccessary
                int n_GPs_estimate = (int)(2* segment_length(j)/std::min(dx,dy));
                segments[j].segment_data_vec.reserve(n_GPs_estimate); //Will remove all of this if it doesn't give any performance gain

                jp = (j + 1) % n_bound;
                p = boundary[j];
                q = boundary[jp];
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
                        Point BI = GP + d;
                        segments[j].segment_data_vec.emplace_back(e.first,(BI - boundary[j]).norm());
                        e.second.BI = GP + d;
                        if (~segments[j].n_is_set){
                            segments[j].n = n;
                            segments[j].n_is_set = true;
                        }
                        prev_dist = curr_dist;
                    }

                }
            }

        }
    }

    void SolidBody::interpolate_fresh_points(fluid::vec4* U_in) {
        //Checks if the old ghost points are outside the boundary. If so they are flagged as fresh points
        //The interpolation will use the values of the old fluid points and boudnary normals,
        // but the body intercept properties (boundary velocity and acceleration from the new fluid points), since this
        //procedure creates no additional computations
        for (auto& e : cell2intercept){
            if (~point_inside(ind2point(e.first.i,e.first.j))){
                interpolate_cell(e.first, e.second.BI, e.second.n, U_in, true);
            }
        }
    }


    void SolidBody::interpolate_solid(fluid::vec4* U_in) {
        using namespace Eigen;
        using namespace std;
        //Handle fresh points

        for (const auto& e : cell2intercept){
            interpolate_cell(e.first, e.second.BI, e.second.n, U_in);
        }
    }


    void SolidBody::interpolate_cell(Cell GP, Point BI, Point n, fluid::vec4* U_in, bool fresh_point) {
        std::vector<fluid::CellStatus> cs(4);
        Matrix4d A_dir = A;
        Matrix4d A_neu = A;
        Point IP = BI * 2 - ind2point(GP.i, GP.j);
        double Delta_l = (IP - BI).norm();
        Cell bottom_left_ind = point2ind(IP.x, IP.y);
        int i_bl = bottom_left_ind.i;
        int j_bl = bottom_left_ind.j;
        Point bottom_left_point = ind2point(i_bl, j_bl);
        double DX = IP.x - bottom_left_point.x;
        double DY = IP.y - bottom_left_point.y;
        Vector4d A_IP{1, DX, DY, DX * DY};

        std::vector<double> rho; rho.reserve(4);
        std::vector<double> u_n; u_n.reserve(4);
        std::vector<double> u_t; u_t.reserve(4);
        std::vector<double> p; p.reserve(4);

        std::vector<double> u_n_BI; u_n_BI.reserve(4);
        std::vector<double> p_BI_derivative; p_BI_derivative.reserve(4);

        Point u_wall{}; Point a_wall{};
        if (type == SolidBodyType::Dynamic){
            //Only estimating one vel and acc for interpolation
            bundary_vel_and_acc(BI, u_wall, a_wall);
        }
        double rho_approx{0};
        bool rho_approx_set{false};
        bool all_fluid{true};
        for (int ii{0}; ii < 2; ii++) {
            for (int jj{0}; jj < 2; jj++) {
                int ind = ii+2*jj;
                cs[ind] = fvm.cell_status[IX(i_bl + ii, j_bl + jj)];

                if (cs[ind] == fluid::CellStatus::Fluid){
                    fluid::vec4 V = fluid::FVM_Solver::conserved2primitive(U_in[IX(i_bl + ii, j_bl + jj)]);
                    u_n[ind] = -(V.u2 * n.x + V.u3 * n.y);
                    u_t[ind] = -V.u2 * n.y + V.u3 * n.x;
                    rho[ind] = V.u1;
                    p[ind] = V.u4;
                }
                else if (cs[ii + 2 * jj] == fluid::CellStatus::Ghost) {
                    all_fluid = false;
                    Point BI_adj = cell2intercept[{i_bl + ii, j_bl + jj}].BI;
                    Point n_adj = cell2intercept[{i_bl + ii, j_bl + jj}].n;
                    DX = BI_adj.x - bottom_left_point.x;
                    DY = BI_adj.y - bottom_left_point.y;
                    A_dir.block<1, 4>(ii + 2 * jj, 0) << 1, DX, DY, DX * DY;
                    A_dir.block<1, 4>(ii + 2 * jj, 0) << 0, n_adj.x, n_adj.y, DY * n_adj.x + DX * n_adj.y;
                    if (type == SolidBodyType::Dynamic){
                        u_n_BI[ind] = u_wall.x*n_adj.x + u_wall.y * n_adj.y;
                        int counter{0};
                        //Approximating the density by an average of surrounding fluid nodes. For pressure gradient calculation
                        if (~rho_approx_set){
                            for (int iii{0};iii<2;iii++){
                                for(int jjj{0};jjj<2;jjj++){
                                    if (cs[iii+2*jjj] == fluid::CellStatus::Fluid){
                                        rho_approx += U_in[IX(i_bl + iii, j_bl + jjj)].u1;
                                        counter++;
                                    }
                                }
                            }
                            rho_approx/=counter;
                            rho_approx_set = true;
                        }
                        p_BI_derivative[ind] = -rho_approx*(a_wall.x*n_adj.x + a_wall.y*n_adj.y);
                    }
                }
                else {
                    std::cerr << "Solid point in bilinear interpolation stencil detected\n"; exit(1);
                }
            }
        }

        Vector4d alpha_dir;
        Vector4d alpha_neu;
        if (~all_fluid){
            A_dir.transposeInPlace();
            A_neu.transposeInPlace();
            alpha_dir = A_dir.partialPivLu().solve(A_IP);
            alpha_neu = A_neu.partialPivLu().solve(A_IP);
        }else{
            alpha_dir = A_inv_T*A_IP;
            alpha_neu = alpha_dir;
        }

        double p_BI_derivative0 = p_BI_derivative[0];
        fluid::vec4 V_GP{};
        double u_n_GP{0};
        double u_t_GP{0};
        u_t_GP = interpolate_neumann_zero_gradient(alpha_neu, std::move(u_t), cs);
        if (type == SolidBodyType::Static) {
            V_GP.u1 = interpolate_neumann_zero_gradient(alpha_neu, std::move(rho), cs);
            u_n_GP = interpolate_dirichlet_zero_value(alpha_dir, std::move(u_n), cs);
            V_GP.u4 = interpolate_neumann_zero_gradient(alpha_neu, std::move(p), cs);
        }
        else if (type == SolidBodyType::Dynamic){
            u_n_GP = interpolate_dirichlet(alpha_dir, std::move(u_n), std::move(u_n_BI), cs);
            V_GP.u4 = interpolate_neumann(alpha_neu, std::move(p), std::move(p_BI_derivative), cs, Delta_l);
            V_GP.u1 = interpolate_neumann_zero_gradient(alpha_neu, std::move(rho), cs);
            V_GP.u1*=(1-Delta_l*p_BI_derivative0);
        } else{
            std::cerr << "Solid type is neither Static, nor Dynamic\n"; exit(1);
        }

        V_GP.u2 = n.x * u_n_GP - n.y * u_t_GP;
        V_GP.u3 = n.y * u_n_GP + n.x * u_t_GP;
        U_in[IX(GP.i, GP.j)] = fluid::FVM_Solver::primitive2conserved(V_GP);
        if (~fresh_point) cell2intercept[GP].p = V_GP.u4 + 0.5*Delta_l*p_BI_derivative0;
    }

    double SolidBody::interpolate_dirichlet(const Vector4d& alpha_dir, std::vector<double>&& phi,
                                            std::vector<double>&& phi_BI ,const std::vector<fluid::CellStatus>& cs) {
        double phi_GP = phi_BI[0];
        for (int i{0}; i < 4;i++){
            if (cs[i] == fluid::CellStatus::Fluid){
                phi_GP -= alpha_dir[i] * phi[i];
            }else{ //Ghost
                phi_GP -= alpha_dir[i] * phi_BI[i];
            }
        }
        return phi_GP;
    }
    double SolidBody::interpolate_neumann(const Vector4d& alpha_neu, std::vector<double>&& phi,
                                          std::vector<double>&& phi_BI_derivative ,
                                          const std::vector<fluid::CellStatus>& cs, double Delta_l) {
        double phi_GP = -Delta_l * phi_BI_derivative[0];
        for (int i{0}; i < 4; i++) {
            if (cs[i] == fluid::CellStatus::Fluid) {
                phi_GP += alpha_neu[i] * phi[i];
            } else { //Ghost
                phi_GP += alpha_neu[i] * phi_BI_derivative[i];
            }
        }
        return phi_GP;
    };
    double SolidBody::interpolate_dirichlet_zero_value(const Eigen::Vector4d& alpha_dir, std::vector<double>&& phi,
                                            const std::vector<fluid::CellStatus>& cs){
        double phi_GP{0};
        for (int i{0}; i < 4;i++){
            if (cs[i] == fluid::CellStatus::Fluid) phi_GP -= alpha_dir[i] * phi[i];
        }
        return phi_GP;
    }

    double SolidBody::interpolate_neumann_zero_gradient(const Eigen::Vector4d& alpha_neu, std::vector<double>&& phi,
                                             const std::vector<fluid::CellStatus>& cs){
        double phi_GP{0};
        for (int i{0}; i < 4; i++) {
            if (cs[i] == fluid::CellStatus::Fluid) {
                phi_GP += alpha_neu[i] * phi[i];
            }
        }
        return phi_GP;
    }

    /*
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
        GP_info this_GP_info;
        for (auto &e: intercepts) {
            //cout << "i = "<<e.first.i << "j = "<<e.first.j<<endl;
            //finding the image point, by using that IP = GP + 2*(BI - GP) = 2*BI - GP
            Point IP = e.second.BI * 2 - ind2point(e.first.i, e.first.j);
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
                n = e.second.n;
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
                //Ghost cells are part of the interpolation stencil. Dirichlet and Neumann conditions at body interface
                //is used in the interpolation
                VM_dirichlet = VM;
                VM_neumann = VM;
                if (S1 == CS::Ghost) {
                    this_GP_info = intercepts[Cell{i, j}];
                    BI_neighbor = this_GP_info.BI;
                    n = this_GP_info.n;
                    DX = BI_neighbor.x - bottom_left_point.x;
                    DY = BI_neighbor.y - bottom_left_point.y;
                    VM_dirichlet.block<1, 4>(0, 0) << 1, DX, DY, DX * DY;
                    VM_neumann.block<1, 4>(0, 0) << 0, n.x, n.y, DY * n.x + DX * n.y;
                }
                if (S2 == CS::Ghost) {
                    this_GP_info = intercepts[Cell{i + 1, j}];
                    BI_neighbor = this_GP_info.BI;
                    n = this_GP_info.n;
                    DX = BI_neighbor.x - bottom_left_point.x;
                    DY = BI_neighbor.y - bottom_left_point.y;
                    VM_dirichlet.block<1, 4>(1, 0) << 1, DX, DY, DX * DY;
                    VM_neumann.block<1, 4>(1, 0) << 0, n.x, n.y, DY * n.x + DX * n.y;
                }
                if (S3 == CS::Ghost) {
                    this_GP_info = intercepts[Cell{i + 1, j + 1}];
                    BI_neighbor = this_GP_info.BI;
                    n = this_GP_info.n;
                    DX = BI_neighbor.x - bottom_left_point.x;
                    DY = BI_neighbor.y - bottom_left_point.y;
                    VM_dirichlet.block<1, 4>(2, 0) << 1, DX, DY, DX * DY;
                    VM_neumann.block<1, 4>(2, 0) << 0, n.x, n.y, DY * n.x + DX * n.y;
                }
                if (S4 == CS::Ghost) {
                    this_GP_info = intercepts[Cell{i, j + 1}];
                    BI_neighbor = this_GP_info.BI;
                    n = this_GP_info.n;
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
                n = e.second.n;
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
     */

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

    void SolidBody::update_lumped_forces(fluid::vec4* U_in) {
        int ip;
        double ds;
        double xi;
        double p_j, p_jp;
        for (int i{0}; i < n_bound; i++) {
            double F1{0};
            double F2{0};
            ip = (i + 1) % n_bound;
            double segment_length_inverse = 1 / (boundary[ip] - boundary[i]).norm();
            segments[i].sort_intercepts();
            std::vector<Segment_data> &v{segments[i].segment_data_vec};
            for (int j{0}; j < v.size() - 1; j++) {
                ds = v[j + 1].s - v[j].s;
                xi = (v[j].s + 0.5 * ds) * segment_length_inverse;
                p_j = cell2intercept[v[j].GP].p;
                p_jp = cell2intercept[v[j].GP].p;
                F1 += 0.5 * (p_jp + p_j) * ds * (1 - xi);
                F2 += 0.5 * (p_jp + p_j) * ds * xi;
            }
            F_boundary[i] += segments[i].n * F1;
            F_boundary[ip] += segments[ip].n * F2;
        }
    }

    void SolidBody::reset_containers(){
        //Call before every timestep to clear containers of variable size
        for (int i{0}; i < n_bound; i++){
            segments[i].segment_data_vec.clear();
            segments[i].n_is_set = false;
            F_boundary[i] = {0,0};
        }
    }

    SolidBody::~SolidBody() {
        delete[] boundary;
        delete[] F_boundary;
    }

    DynamicRigid::DynamicRigid(fluid::FVM_Solver &fvm, std::vector<Point>&& boundary_in,Point CM, double M, double I)
    : SolidBody(fvm,std::move(boundary_in),SolidBodyType::Dynamic){
        y[0] = CM.x;
        y[1] = CM.y;
        y[2] = 0;
        y[3] = 0;
        y[4] = 0;
        y[5] = 0;
    }

    void DynamicRigid::update_total_fluid_force_and_moment(){
        F_fluid = {0,0};
        tau_fluid = 0;
        Point r{};
        for (int i{0};i< n_bound;i++){
            F_fluid += F_boundary[i];
            r = {boundary[i].x- y[0], boundary[i].y - y[1]};
            tau_fluid += r.cross(F_boundary[i]);
        }
    }


    Vector6d DynamicRigid::evaluate_f(Vector6d y_in){
        //rhs of the state vector derivative dy/dt = f = [u_CM, v_CM, Fx/M, Fy/M, omega, tau/I]^T
        Vector6d f;
        //std::pair<Point,double> F_S_tau_S{eval_solid_forces_and_moment()};
        f[0] = y_in[2];
        f[1] = y_in[3];
        f[2] = (F_fluid.x + F_solid.x)/M;
        f[3] = (F_fluid.y + F_solid.y)/M;
        f[4] = y_in[5];
        f[5] = (tau_fluid + tau_solid)/I;
        return f;
    }


    Vector6d DynamicRigid::RK4_step(double dt){
        k1 = dt * evaluate_f(y);
        k2 = dt * evaluate_f(y + k1/2);
        k3 = dt * evaluate_f(y + k2/2);
        k4 = dt * evaluate_f(y + k3);
        return y + (k1 + 2*(k2+k3) + k4)/6;
    }





}
