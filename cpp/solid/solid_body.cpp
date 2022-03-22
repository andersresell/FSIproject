//
// Created by anders on 2/1/22.
//

#include "solid_body.hpp"

namespace solid {

    using namespace Eigen;
    using namespace std;


    SolidBody::SolidBody(fluid::FVM_Solver &fvm, std::vector<Point>&& boundary_in, SolidBodyType type)
            : fvm{fvm}, n_bound{(int)boundary_in.size()}, ni{fvm.ni}, nj{fvm.nj},
              dx{fvm.L_x / ni}, dy{fvm.L_y / nj}, type{type}, first_timestep{true}{
        boundary = new Point[n_bound];
        F_boundary = new Point[n_bound];
        //segments = new Segment[n_bound];
        for (int i{0}; i < boundary_in.size(); i++) {
            boundary[i] = boundary_in[i];
        }
        Cell::nj = nj;

        A << 1, 0, 0, 0,
                1, dx, 0, 0,
                1, dx, dy, dx * dy,
                1, 0, dy, 0;
        A_inv_T = A.transpose().inverse();
        timestep = 0;
    }


    void SolidBody::step(fluid::vec4* U_in, double dt, bool update_solid_pos){
        if (first_timestep && type == SolidBodyType::Static) {
            find_solid_cells();
            flag_static();
            find_ghost_cells();
            find_intercepts(cell2intercept);
        }
        else if (type == SolidBodyType::Dynamic){
            if (first_timestep){
                find_solid_cells();
                find_ghost_cells();
                find_intercepts(cell2intercept);
            }else{
                if (update_solid_pos) {
                    update_lumped_forces(U_in);
                    step_solid_body(dt);
                    interpolate_fresh_points(U_in);
                    //write_fresh_points(++timestep);
                    find_solid_cells();
                    reset_containers();
                    find_ghost_cells();
                    find_intercepts(cell2intercept);

                }
            }
        }
        interpolate_ghost_points(U_in);
        first_timestep = false;
    }

    void SolidBody::find_solid_cells() {
        //DEBUGGING
        /* solid_cells.clear();
         Point p{};
         for (int i{2}; i < ni + 2; i++) {
             for (int j{2}; j < nj + 2; j++) {
                 p = ind2point(i, j);
                 if (point_inside(p)) {
                     solid_cells.emplace(Cell{i, j});
                     //solid_cells.emplace(Cell{i,j});
                     fvm.cell_status[IX(i, j)] = fluid::CellStatus::Solid;
                 } else {
                     fvm.cell_status[IX(i, j)] = fluid::CellStatus::Fluid;
                 }
             }
         }
     }
 */
        Point p{};
        if (first_timestep) {
            //During the first timestep all cells has to be detected
            for (int i{2}; i < ni + 2; i++) {
                for (int j{2}; j < nj + 2; j++) {
                    p = ind2point(i, j);
                    if (point_inside(p)) {
                        //solid_cells.push_back(Cell{i, j});
                        solid_cells.emplace(Cell{i, j});
                        fvm.cell_status[IX(i, j)] = fluid::CellStatus::Solid;
                    }
                }
            }
        } else if (type == SolidBodyType::Dynamic) {
            for (const auto &e: cell2intercept) {
                //The cells in the old GP's are set to fluid
                fvm.cell_status[IX(e.first.i, e.first.j)] = fluid::CellStatus::Fluid;
            }
            //Exploiting that the CFL condition only allows the boundary to travel one cell to find the new solid cells
            int i, j;
            for (const auto &e: cell2intercept) {
                for (int ii{-1}; ii <= 1; ii++) {
                    for (int jj{-1}; jj <= 1; jj++) {
                        i = e.first.i;
                        j = e.first.j;
                        //Continuing if the indices are outside the computational domain or if the cell is allready set
                        if (fvm.cell_status[IX(i + ii, j + jj)] == fluid::CellStatus::Solid || !cell_within_grid(i + ii, j + jj))
                            continue;
                        p = ind2point(i + ii, j + jj);
                        if (point_inside(p)) {
                            fvm.cell_status[IX(i + ii, j + jj)] = fluid::CellStatus::Solid;
                            solid_cells.emplace(Cell{i + ii, j + jj});
                        } else {
                            solid_cells.erase({i + ii, j + jj});
                            fvm.cell_status[IX(i + ii, j + jj)] = fluid::CellStatus::Fluid;
                        }
                    }
                }
            }
        } else {
            std::cerr << "Attempting to find solid cells for a non dynamic body after the first timestep. "
                         "<<Unnecessary operation\n";
            exit(1);
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
        /*
        std::vector<Cell> solid_cells_copy;
        solid_cells_copy.reserve(solid_cells.size());
        std::copy(solid_cells.begin(), solid_cells.end(), solid_cells_copy.begin());*/

        bool ghost;
        int i, j;
        for (Cell c: solid_cells) {
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

    void SolidBody::find_intercepts(std::map<Cell, GP_data>& intercept_map) {
        using namespace std;
        Point r_p;
        Point r_q;
        Point n;
        Point point;
        Point q;
        Point p;
        Point d;
        Point shortest_n;
        Point shortest_d;
        int jp;
        double curr_dist;
        for (auto& e : intercept_map){
            point = ind2point(e.first.i,e.first.j);
            double prev_dist = INF;
            double shortest_dist = INF;
            bool intercept_found{false};
            for (int j{0}; j < n_bound; j++) {
                //A different idea is to store the length from the previous iteration and use this as an estimate.
                /*segments[j].segment_data_vec.clear(); //Maybe unneccessary
                int n_GPs_estimate = (int)(2* segment_length(j)/std::min(dx,dy));
                segments[j].segment_data_vec.reserve(n_GPs_estimate); //Will remove all of this if it doesn't give any performance gain*/

                jp = (j + 1) % n_bound;
                p = boundary[j];
                q = boundary[jp];
                r_p = p - point;
                r_q = q - point;
                n = {q.y - p.y, -(q.x - p.x)};
                n.normalize();
                d = n * r_p.dot(n);
                curr_dist = d.norm();

                if (curr_dist < shortest_dist){ //Used for special case where no intersection is detected
                    shortest_d = d;
                    shortest_n = n;
                }
                if (r_p.cross(n) >= 0 && n.cross(r_q) >= 0) { //intersection is on the line segment
                    if (curr_dist < prev_dist) {
                        intercept_found = true;
                        Point BI = point + d;
                        e.second.BI = point + d;
                        e.second.n = n;

                        /*if (!fresh_point) { //Ghost point
                            segments[j].segment_data_vec.emplace_back(e.first, (BI - boundary[j]).norm());
                            if (!segments[j].n_is_set) {
                                segments[j].n = n;
                                segments[j].n_is_set = true;
                            }
                        }*/
                        prev_dist = curr_dist;
                    }

                }
            }
            if (!intercept_found){
                e.second.BI = point + shortest_d;
                e.second.n = shortest_n;
            }

        }
    }

    void SolidBody::interpolate_fresh_points(fluid::vec4* U_in) {
        //DESCRIPTION OUT OF DATE!
        //Checks if the old ghost points are outside the boundary. If so they are flagged as fresh points
        //The interpolation will use the values of the old fluid points and boundary normals,
        // but the body intercept properties (boundary velocity and acceleration from the new fluid points), since this
        //procedure creates no additional computations
        for (auto& e : cell2intercept){
            if (!point_inside(ind2point(e.first.i,e.first.j))){
                cell2intercept_FP[e.first]; //Detecting the fresh points, by checking which old GP's are now outside the boundary
            }
        }
        find_intercepts(cell2intercept_FP);
        for (auto& e : cell2intercept_FP) {
            interpolate_cell(e.first, cell2intercept_FP, U_in,true);
        }
    }

    void SolidBody::interpolate_ghost_points(fluid::vec4* U_in) {
        for (auto& e : cell2intercept){
            interpolate_cell(e.first, cell2intercept, U_in);
        }

    }


    void SolidBody::interpolate_cell(Cell point, std::map<Cell, GP_data>& intercept_map,
                                     fluid::vec4* U_in,bool fresh_point) {
        using namespace std;
        Point BI = intercept_map[point].BI;
        Point n = intercept_map[point].n;
        std::array<fluid::CellStatus,4> cs{};
        Matrix4d A_dir = A;
        Matrix4d A_neu = A;;
        Point IP{};
        double Delta_l;
        Cell bottom_left_ind{};
        if (fresh_point) {
            bottom_left_ind = point;
            IP = ind2point(point.i, point.j); //IP denotes the fresh point
            //Delta_l = (IP - BI).norm() * 2;
            if (n.x < 0) bottom_left_ind.i--;
            if (n.y < 0) bottom_left_ind.j--;
           // cout << "FPpoint "<<point<<", bl"<<bottom_left_ind<<", n "<<n<<endl;
        } else{
            IP = BI * 2 - ind2point(point.i, point.j);
            Delta_l = (IP - BI).norm() * 2;
            //cout << "IP = "<<IP<<", point "<<point<<", BI "<<BI<<endl;
            //Cell bottom_left_ind = point2ind(IP.x, IP.y);
            bottom_left_ind = point2ind(IP.x, IP.y);
            if (Delta_l == 0){
                if (n.x < 0) bottom_left_ind.i--;
                if (n.y < 0) bottom_left_ind.j--;
            }
        }
        int i_bl = bottom_left_ind.i;
        int j_bl = bottom_left_ind.j;
        Point bottom_left_point = ind2point(i_bl, j_bl);

        double DX = IP.x - bottom_left_point.x;
        double DY = IP.y - bottom_left_point.y;
        Vector4d A_IP{1, DX, DY, DX * DY};

        std::array<double,4> rho{};
        std::array<double,4> u_n{};
        std::array<double,4> u_t{};
        std::array<double,4> p{};
        std::array<double,4> u_n_BI_adj{};
        std::array<double,4> p_derivative_BI_adj{};
        //std::array<double,4> rho_derivative_BI_adj{};
        double u_n_BI, p_derivative_BI;
        double rho_approx{0};
        double p_approx{0};

        Point u_wall{};
        Point a_wall{};
        if (type == SolidBodyType::Dynamic) {
            //Only estimating one vel and acc for interpolation
            bundary_vel_and_acc(BI, u_wall, a_wall);
            //int counter{0};
            //Approximating the density by an average of surrounding fluid nodes. For pressure gradient calculation
            for (int ii{0}; ii < 2; ii++) {
                for (int jj{0}; jj < 2; jj++) {
                    //if (cs[ii + 3 * jj - 2 * ii * jj] == fluid::CellStatus::Fluid) {
                        rho_approx += U_in[IX(i_bl + ii, j_bl + jj)].u1;
                        //p_approx += fluid::FVM_Solver::conserved2primitive(U_in[IX(i_bl + ii, j_bl + jj)]).u4; //Experimental
                        //counter++;
                    //}
                }
            }
            rho_approx *= 0.25;
            //p_approx /= counter;
            //cout << "GP "<<point<<", n"<<n<<", rho_approx "<<rho_approx<< ", bl "<<bottom_left_ind<<endl;
            u_n_BI = u_wall.x * n.x + u_wall.y * n.y;
            p_derivative_BI = -rho_approx * (a_wall.x * n.x + a_wall.y * n.y);
            //cout << "a_wall "<<a_wall<<"u_wall "<<u_wall<<", rho_approx = "<<rho_approx<<endl;
        }
        bool all_fluid{true};
        for (int jj{0}; jj < 2; jj++) {
            for (int ii{0}; ii < 2; ii++) {
                int ind = ii + 3*jj - 2*ii*jj;
                cs[ind] = fvm.cell_status[IX(i_bl + ii, j_bl + jj)];
                //In the case that the cell is an external GP, the value is interpolated normally for simplicity
                if (cs[ind] == fluid::CellStatus::Fluid || !cell_within_grid(i_bl+ii,j_bl+jj)){
                    fluid::vec4 V = fluid::FVM_Solver::conserved2primitive(U_in[IX(i_bl + ii, j_bl + jj)]);
                    //if (!cell_within_grid(i_bl+ii,j_bl+jj)) cout << i_bl+ii<<","<<j_bl+jj<<"point "<<point<<", BI "<<BI<<", n "<<n<<endl;
                    u_n[ind] = V.u2 * n.x + V.u3 * n.y;
                    u_t[ind] = -V.u2 * n.y + V.u3 * n.x;
                    rho[ind] = V.u1;
                    p[ind] = V.u4;
                }
                else if (cs[ind] == fluid::CellStatus::Ghost) {
                    //cout << "point "<<point<<", n "<<n<<", bl "<<bottom_left_ind<<", dist "<< (IP-BI).norm()<<endl;
                    all_fluid = false;
                    Point BI_adj = intercept_map[{i_bl + ii, j_bl + jj}].BI;
                    Point n_adj = intercept_map[{i_bl + ii, j_bl + jj}].n;
                    //if (n_adj.x == 0 && n_adj.y == 0) cout << "FP"<<fresh_point<<", point "<<point<<", n "<<n<<", bl "<<bottom_left_ind<<", dist "<< (IP-BI).norm()<<endl;

                    double DX_BI_adj = BI_adj.x - bottom_left_point.x;
                    double DY_BI_adj = BI_adj.y - bottom_left_point.y;
                    A_dir.block<1, 4>(ind, 0) << 1, DX_BI_adj, DY_BI_adj, DX_BI_adj * DY_BI_adj;
                    A_neu.block<1, 4>(ind, 0) << 0, n_adj.x, n_adj.y, DY_BI_adj * n_adj.x + DX_BI_adj * n_adj.y;
                    if (type == SolidBodyType::Dynamic){
                        u_n_BI_adj[ind] = u_wall.x*n_adj.x + u_wall.y * n_adj.y;
                        p_derivative_BI_adj[ind] = -rho_approx*(a_wall.x*n_adj.x + a_wall.y*n_adj.y);
                        //rho_derivative_BI_adj[ind] = p_derivative_BI_adj[ind] * rho_approx/p_approx;
                        //cout << "p_deriv "<<p_derivative_BI_adj[ind]<<", p_approx "<<p_approx<<", rho_deriv "<<rho_derivative_BI_adj[ind]<<endl;
                    }
                }
                else if (cs[ind] == fluid::CellStatus::Solid){
                    cout << "FP" << fresh_point<<", point "<<point<<", n "<<n<<", bl "<<bottom_left_ind<<endl;
                    std::cerr << "Solid point in bilinear interpolation stencil detected\n"; //exit(1);
                } else{

                    std::cerr << "Unidentified point in bilinear interpolation stencil detected\n"; //exit(1);
                }
            }
        }
        Vector4d alpha_dir;
        Vector4d alpha_neu;
        if (!all_fluid){
            A_dir.transposeInPlace();
            A_neu.transposeInPlace();
            alpha_dir = A_dir.partialPivLu().solve(A_IP);
            alpha_neu = A_neu.partialPivLu().solve(A_IP);
        }else{
            alpha_dir = A_inv_T*A_IP;
            alpha_neu = alpha_dir;
        }
        //cout <<"alpha_n "<<alpha_neu<<"\n" "alpha_d "<<alpha_dir<<endl;
        //for (int i{0};i<4;i++)  cout << "rho "<<rho[i] <<", p "<<p[i]<<endl;
        fluid::vec4 V_point{};
        double u_n_point;
        double u_t_point{0};

        if (fresh_point){
            u_t_point = interpolate_neumann(alpha_neu, u_t, cs);
            u_n_point = interpolate_dirichlet(alpha_dir, u_n, cs, u_n_BI_adj);
            V_point.u4 = interpolate_neumann(alpha_neu, p, cs, p_derivative_BI_adj);
            V_point.u1 = interpolate_neumann(alpha_neu, rho, cs);//, rho_derivative_BI_adj);
        }else{
            u_t_point = interpolate_neumann(alpha_neu, u_t, cs);
            if (type == SolidBodyType::Static) {
                V_point.u1 = interpolate_neumann(alpha_neu, rho, cs);
                u_n_point = - interpolate_dirichlet(alpha_dir, u_n, cs);
                V_point.u4 = interpolate_neumann(alpha_neu, p, cs);
            }
            else if (type == SolidBodyType::Dynamic){
                u_n_point = 2*u_n_BI - interpolate_dirichlet(alpha_dir, u_n, cs, u_n_BI_adj);
                V_point.u4 = - Delta_l*p_derivative_BI + interpolate_neumann(alpha_neu, p, cs, p_derivative_BI_adj);
                V_point.u1 = interpolate_neumann(alpha_neu, rho, cs);//, rho_derivative_BI_adj);
                V_point.u1*=(1-Delta_l*p_derivative_BI);
            } else{
                std::cerr << "Solid type is neither Static, nor Dynamic\n"; //exit(1);
            }
        }
        V_point.u2 = n.x * u_n_point - n.y * u_t_point;
        V_point.u3 = n.y * u_n_point + n.x * u_t_point;
        U_in[IX(point.i, point.j)] = fluid::FVM_Solver::primitive2conserved(V_point);
        /*if (V_point.isnan()){
            cout << "FP"<<fresh_point<<", V_point "<<V_point << endl;
            cout << "a_dir "<<alpha_dir<<", a_neu "<<alpha_neu<<", A_IP "<<A_IP<<endl;
            cout <<"A_dir "<<A_dir<<endl<<"A_neu "<<A_neu<<endl;
            cout << "cell" << point<<", bl "<<bottom_left_ind<<endl;
            for (int i{0};i<4;i++) cout << "cs "<<(int)cs[i]<<endl;
            for (int jj{0}; jj < 2; jj++) {
                for (int ii{0}; ii < 2; ii++) {
                    cout << "CS, " << (int) fvm.cell_status[IX(i_bl + ii, j_bl + jj)] <<endl;
                    cout << "U_in "<<U_in[IX(i_bl + ii, j_bl + jj)]<<endl;
                }
            }
            cout <<endl;
            cout << "GP point "<<ind2point(point.i,point.j)<<", BI" <<BI <<", IP "<<IP<<endl;
            cout << "dx "<< dx<<endl;
            exit(1);
        }*/
    }

    double SolidBody::interpolate_dirichlet(const Vector4d& alpha_dir, std::array<double,4> phi,
                                        std::array<fluid::CellStatus,4> cs, std::array<double,4> phi_BI_adj) {

        double phi_point{0};
        for (int i{0}; i < 4; i++) {
            if (cs[i] == fluid::CellStatus::Fluid) {
                phi_point += alpha_dir[i] * phi[i];
            } else { //Ghost
                phi_point += alpha_dir[i] * phi_BI_adj[i];
            }
        }
        return phi_point;
    }
    double SolidBody::interpolate_neumann(const Eigen::Vector4d& alpha_neu, std::array<double,4> phi,
                                          std::array<fluid::CellStatus,4> cs,
                                          std::array<double,4> phi_derivative_BI_adj){
        double phi_point{0};
        for (int i{0}; i < 4; i++) {
            if (cs[i] == fluid::CellStatus::Fluid) {
                phi_point += alpha_neu[i] * phi[i];
            } else { //Ghost
                phi_point += alpha_neu[i] * phi_derivative_BI_adj[i];
            }
        }
        return phi_point;
    };

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
/*
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
            for (size_t j{0}; j < v.size() - 1; j++) {
                //cout << "v_size "<<v.size()<<", i "<<i<<",j "<<j<<endl;
                ds = v[j + 1].s - v[j].s;
                xi = (v[j].s + 0.5 * ds) * segment_length_inverse;
                p_j = cell2intercept[v[j].GP].p;
                p_jp = cell2intercept[v[j].GP].p;
                F1 -= 0.5 * (p_jp + p_j) * ds * (1 - xi);
                F2 -= 0.5 * (p_jp + p_j) * ds * xi;
            }
            F_boundary[i] += segments[i].n * F1;
            F_boundary[ip] += segments[ip].n * F2;
        }
        for (int  i{0};i<n_bound;i++) cout<< "force b_i = "<<F_boundary[i]<<endl;
    }*/
    void SolidBody::update_lumped_forces(fluid::vec4* U_in) {
        //uses the trapezoidal method for integrationand lumps the forces into each boudnary node
        double xi;
        size_t ip;
        double ds{std::min(dx, dy)};
        double dF;
        Vector4d alpha, P;
        double DX, DY;
        for (size_t i{0}; i < n_bound; i++) {
            ip = (i + 1) % n_bound;
            Point p = boundary[i];
            Point q = boundary[ip];
            double F_p{0};
            double F_q{0};
            double seg_length = (q - p).norm();
            Point n = {q.y - p.y, -(q.x - p.x)};
            n.normalize();
            Point t = {-n.y, n.x};
            auto n_ds = static_cast<size_t>(seg_length / ds) + 1;
            //cout <<"p "<<p<<",q "<<q<<endl;
            ds = seg_length / (double) n_ds;
            //cout << "dx "<<dx<<", dy "<<dy<<endl;
            //cout << "ds " <<ds<<", seg_len "<<seg_length<< ", _nds" <<n_ds<<endl;
            std::vector<double> pressures(n_ds+1);
            //calculating pressures at all evaluation points
            for (int j{0};j< n_ds+1;j++){
                double s = ds*j;
                Point point = boundary[i] + t*s;
                if (point_within_grid(point)) {
                    Cell bl = point2ind(point.x, point.y);
                    Point bl_point = ind2point(bl.i, bl.j);
                    P[0] = fluid::FVM_Solver::calc_P(U_in[IX(bl.i, bl.j)]);
                    P[1] = fluid::FVM_Solver::calc_P(U_in[IX(bl.i + 1, bl.j)]);
                    P[2] = fluid::FVM_Solver::calc_P(U_in[IX(bl.i + 1, bl.j + 1)]);
                    P[3] = fluid::FVM_Solver::calc_P(U_in[IX(bl.i, bl.j + 1)]);
                    DX = point.x - bl_point.x;
                    DY = point.y - bl_point.y;
                    Vector4d A_point = {1, DX, DY, DX * DY};
                    alpha = A_inv_T * A_point;
                    pressures[j] = alpha.dot(P);
                }
            }
            //integrating the pressures and load lumping them
            for (int j{0};j<n_ds;j++){
                dF = -0.5*(pressures[j+1] + pressures[j])*ds;
                xi = ds*(j + 0.5)/seg_length;
                F_p += dF * (1 - xi);
                F_q += dF * xi;
                //cout << "xi "<<xi<<endl;
            }

            /*for (int j{0}; j < n_ds; j++) {
                double s = ds*(0.5 + (double) j);
                Point point = boundary[i] + t * s;
                xi = s / seg_length;
                Cell bl = point2ind(point.x, point.y);
                Point bl_point = ind2point(bl.i, bl.j);
                //cout << "U_00 "<< U_in[IX(bl.i, bl.j)]<<endl;
                //cout << "U_10 "<< U_in[IX(bl.i+1, bl.j)]<<endl;
                //cout << "U_11 "<< U_in[IX(bl.i+1, bl.j+1)]<<endl;
                //cout << "U_01 "<< U_in[IX(bl.i, bl.j+1)]<<endl;

                P[0] = fluid::FVM_Solver::calc_P(U_in[IX(bl.i, bl.j)]);
                P[1] = fluid::FVM_Solver::calc_P(U_in[IX(bl.i + 1, bl.j)]);
                P[2] = fluid::FVM_Solver::calc_P(U_in[IX(bl.i + 1, bl.j + 1)]);
                P[3] = fluid::FVM_Solver::calc_P(U_in[IX(bl.i, bl.j + 1)]);
                DX = point.x - bl_point.x;
                DY = point.y - bl_point.y;
                Vector4d A_point = {1, DX, DY, DX * DY};
                alpha = A_inv_T * A_point;
                double P_point = alpha.dot(P);
                //for (int k{0};k<4;k++) cout << "P_k "<< P[k] <<endl;
                //cout <<"P_point "<< P_point<<endl;
                dF = -ds * P_point;
                F_p += dF * (1 - xi);
                F_q += dF * xi;
            }*/
            F_boundary[i] += n * F_p;
            F_boundary[ip] += n * F_q;
        }
        //for (int i{0};i<n_bound;i++)cout << "b_i "<<boundary[i]<< ", F_b "<<F_boundary[i]<<endl;
    }

    void SolidBody::write_fresh_points(int n) const{
        std::ofstream ost{"python/output_folders/"+fvm.output_folder+"/debug_fresh_points_t" + std::to_string(n) + ".csv"};
        if (!ost) std::cerr << "error, couldn't open solid debug fresh points csv file\n";
        ost << "#x_i,y_i\n";
        for (const auto& e: cell2intercept_FP){
            solid::Point p = {(e.first.i - 1.5) * fvm.dx, (e.first.j - 1.5) * fvm.dy};
            ost << p.x << ',' << p.y << '\n';
        }
    }

    void SolidBody::reset_containers(){
        //Call before every timestep to clear containers of variable size

        cell2intercept.clear();
        cell2intercept_FP.clear();
        for (int i{0}; i < n_bound; i++){
            //segments[i].segment_data_vec.clear();
            //segments[i].n_is_set = false;
            F_boundary[i] = {0,0};
        }
    }

    SolidBody::~SolidBody() {
        delete[] boundary;
        delete[] F_boundary;
        //delete[] segments;
    }



}
