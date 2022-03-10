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
        segments = new Segment[n_bound];
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


    void SolidBody::step(fluid::vec4* U_in, double dt, bool update_solid_pos){
        if (type == SolidBodyType::Static && first_timestep) {
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
                    for (int i{0};i<n_bound;i++){
                        cout << "B_i: "<<boundary[i]<<endl;
                    }
                    cout <<endl;
                    interpolate_fresh_points(U_in);//Shady, some containers should maybe be reset before
                    reset_containers();
                    find_solid_cells();
                    find_ghost_cells();
                    find_intercepts(cell2intercept);
                }
            }
        }
        //interpolate_solid_old(U_in);
        interpolate_solid(U_in);
        first_timestep = false;
    }

    void SolidBody::find_solid_cells() {
        //DEBUGGING
        Point p;
        for (int i{2}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                p = ind2point(i, j);
                if (point_inside(p)) {
                    //solid_cells.push_back(Cell{i, j});
                    solid_cells.emplace(Cell{i,j});
                    fvm.cell_status[IX(i, j)] = fluid::CellStatus::Solid;
                }
                else {
                    fvm.cell_status[IX(i, j)] = fluid::CellStatus::Fluid;
                }
            }
        }
/*
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
                    }
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
        */
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

    void SolidBody::find_intercepts(std::map<Cell, GP_data>& intercept_map, bool fresh_point) {
        using namespace std;
        Point r_p;
        Point r_q;
        Point n;
        Point point;
        Point q;
        Point p;
        Point d;
        int jp;
        double curr_dist;
        double prev_dist;
        for (auto& e : intercept_map){
            point = ind2point(e.first.i,e.first.j);
            prev_dist = INF;
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
                // cout << "c1 "<<n.cross(r1)<<", c2 "<< r2.cross(n) <<endl;
                if (r_p.cross(n) >= 0 && n.cross(r_q) >= 0) { //intersection is on the line segment
                    d = n * r_p.dot(n);
                    curr_dist = d.norm();
                    //cout << "curr_dist = "<<curr_dist<<endl;
                    if (curr_dist < prev_dist) {
                        Point BI = point + d;
                        e.second.BI = point + d;
                        e.second.n = n;
                        if (!fresh_point) { //Ghost point
                            segments[j].segment_data_vec.emplace_back(e.first, (BI - boundary[j]).norm());
                            if (~segments[j].n_is_set) {
                                segments[j].n = n;
                                segments[j].n_is_set = true;
                            }
                        }
                        prev_dist = curr_dist;
                    }

                }
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
            if (~point_inside(ind2point(e.first.i,e.first.j))){
                cell2intercept_FP[e.first]; //Detecting the fresh pointsm, by checking which old GP's are now outside the boundary
            }
        }
        find_intercepts(cell2intercept_FP,true);
        for (auto& e : cell2intercept_FP) {

        }

        //interpolate_cell(e.first, e.second.BI, e.second.n, U_in, true);
    }


    void SolidBody::interpolate_solid(fluid::vec4* U_in) {
        using namespace Eigen;
        using namespace std;
        //Handle fresh points

        for (auto& e : cell2intercept){
            interpolate_cell(e.first, cell2intercept, U_in);
        }
    }


    void SolidBody::interpolate_cell(Cell point, std::map<Cell, GP_data>& intercept_map,
                                     fluid::vec4* U_in,bool fresh_point) {
        using namespace std;
        Point BI = intercept_map[point].BI;
        Point n = intercept_map[point].n;
        std::vector<fluid::CellStatus> cs(4);
        Matrix4d A_dir = A;
        Matrix4d A_neu = A;
        Point IP = BI * 2 - ind2point(point.i, point.j);
        double Delta_l = (IP - BI).norm()*2;
        Cell bottom_left_ind = point2ind(IP.x, IP.y);
        //cout <<"Fresh point: "<<fresh_point<< ", GP: "<<GP<<"BI: "<<BI<<"n: "<<n <<", bl " << bottom_left_ind<<endl;
        int i_bl = bottom_left_ind.i;
        int j_bl = bottom_left_ind.j;
        Point bottom_left_point = ind2point(i_bl, j_bl);
        double DX = IP.x - bottom_left_point.x;
        double DY = IP.y - bottom_left_point.y;
        Vector4d A_IP{1, DX, DY, DX * DY};

        std::vector<double> rho(4);
        std::vector<double> u_n(4);
        std::vector<double> u_t(4);
        std::vector<double> p(4);
        std::vector<double> u_n_BI_adj(4);
        std::vector<double> p_derivative_BI_adj(4);
        double rho_BI_derivative, u_n_BI, p_derivative_BI;



        Point u_wall{}; Point a_wall{};
        if (type == SolidBodyType::Dynamic){
            //Only estimating one vel and acc for interpolation
            bundary_vel_and_acc(BI, u_wall, a_wall);
        }
        double rho_approx{0};
        bool rho_approx_set{false};
        bool all_fluid{true};
        for (int jj{0}; jj < 2; jj++) {
            for (int ii{0}; ii < 2; ii++) {
                int ind = ii + 3*jj - 2*ii*jj;
                cs[ind] = fvm.cell_status[IX(i_bl + ii, j_bl + jj)];
                //cout <<"ind: "<<ind<< ", ("<<i_bl+ii<<","<<j_bl+jj<<")";
                //cout << ", S: " <<(int)cs[ind] << ", ind: "<<IX(i_bl + ii, j_bl + jj)<<endl;
                if (cs[ind] == fluid::CellStatus::Fluid){
                    fluid::vec4 V = fluid::FVM_Solver::conserved2primitive(U_in[IX(i_bl + ii, j_bl + jj)]);
                    u_n[ind] = V.u2 * n.x + V.u3 * n.y;
                    u_t[ind] = -V.u2 * n.y + V.u3 * n.x;
                    /*cout << "tull = "<<-V.u2 * n.y + V.u3 * n.x<<endl;
                    cout <<"n "<< n << endl;
                    cout << "V "<<V<<endl;
                    cout << "ind = "<< ind<< ", u_t_ind "<< u_t[ind]<<endl;*/
                    rho[ind] = V.u1;
                    p[ind] = V.u4;
                }
                else if (cs[ind] == fluid::CellStatus::Ghost) {
                    all_fluid = false;
                    Point BI_adj = cell2intercept[{i_bl + ii, j_bl + jj}].BI;
                    Point n_adj = cell2intercept[{i_bl + ii, j_bl + jj}].n;
                    double DX_BI = BI_adj.x - bottom_left_point.x;
                    double DY_BI = BI_adj.y - bottom_left_point.y;
                    A_dir.block<1, 4>(ind, 0) << 1, DX_BI, DY_BI, DX_BI * DY_BI;
                    A_neu.block<1, 4>(ind, 0) << 0, n_adj.x, n_adj.y, DY_BI * n_adj.x + DX_BI * n_adj.y;
                    if (type == SolidBodyType::Dynamic){
                        u_n_BI_adj[ind] = u_wall.x*n_adj.x + u_wall.y * n_adj.y;
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
                        p_derivative_BI_adj[ind] = -rho_approx*(a_wall.x*n_adj.x + a_wall.y*n_adj.y);
                    }
                }
                else if (cs[ind] == fluid::CellStatus::Solid){
                    std::cerr << "Solid point in bilinear interpolation stencil detected\n"; exit(1);
                } else{
                    std::cerr << "Unidentified point in bilinear interpolation stencil detected\n"; exit(1);
                }
            }
        }
      //  cout << "all fluid: "<<all_fluid;
        //cout << endl;
        Vector4d alpha_dir;
        Vector4d alpha_neu;
        if (!all_fluid){
            A_dir.transposeInPlace();
            A_neu.transposeInPlace();
            alpha_dir = A_dir.partialPivLu().solve(A_IP);
            alpha_neu = A_neu.partialPivLu().solve(A_IP);

            //cout << "not all fluid\n";
            //cout << "alpha_dir = " << alpha_dir << endl;
            //cout << "alpha_neu = " << alpha_neu << endl<<endl;
        }else{
            alpha_dir = A_inv_T*A_IP;
            alpha_neu = alpha_dir;
            //cout << "all fluid\n";
            //cout << "alpha = " << alpha_dir << endl<<endl;
        }


        double p_BI_derivative0 = p_derivative_adj[0];
        fluid::vec4 V_point{};
        double u_n_point{0};
        double u_t_point{0};
        //cout << "alpha_neu = "<<alpha_neu;
        //cout << "u_t = ";
        //for (int r{0};r<4;r++) cout << u_t.at(r)<< ",";
        //cout <<endl;
        u_t_point = interpolate_neumann_zero_gradient(alpha_neu, std::move(u_t), cs);
        //cout << "u_t_GP= "<<u_t_GP<<endl;
        if (type == SolidBodyType::Static) {
            V_point.u1 = interpolate_neumann_zero_gradient(alpha_neu, std::move(rho), cs);
            u_n_point = interpolate_dirichlet_zero_value(alpha_dir, std::move(u_n), cs);
            V_point.u4 = interpolate_neumann_zero_gradient(alpha_neu, std::move(p), cs);
        }
        else if (type == SolidBodyType::Dynamic){
            u_n_point = interpolate_dirichlet(alpha_dir, std::move(u_n), std::move(u_n_BI_adj), cs);
            V_point.u4 = interpolate_neumann(alpha_neu, std::move(p), std::move(p_derivative_BI_adj), cs, Delta_l);
            V_point.u1 = interpolate_neumann_zero_gradient(alpha_neu, std::move(rho), cs);
            V_point.u1*=(1-Delta_l*p_BI_derivative0);
        } else{
            std::cerr << "Solid type is neither Static, nor Dynamic\n"; exit(1);
        }

        V_point.u2 = n.x * u_n_point - n.y * u_t_point;
        V_point.u3 = n.y * u_n_point + n.x * u_t_point;
        U_in[IX(point.i, point.j)] = fluid::FVM_Solver::primitive2conserved(V_point);
        if (!fresh_point) cell2intercept[point].p = V_point.u4 + 0.5*Delta_l*p_BI_derivative0;
    }

    double SolidBody::interpolate_dirichlet(const Vector4d& alpha_dir, std::vector<double>&& phi,
                            std::vector<double>&& phi_BI ,const std::vector<fluid::CellStatus>& cs, bool fresh_point) {
        if (!fresh_point) {
            double phi_GP = phi_BI[0];
            for (int i{0}; i < 4; i++) {
                if (cs[i] == fluid::CellStatus::Fluid) {
                    phi_GP -= alpha_dir[i] * phi[i];
                } else { //Ghost
                    phi_GP -= alpha_dir[i] * phi_BI[i];
                }
            }
            return phi_GP;
        }else{

        }
    }
    double SolidBody::interpolate_neumann(const Vector4d& alpha_neu, std::vector<double>&& phi,
                                          std::vector<double>&& phi_BI_derivative ,
                                          const std::vector<fluid::CellStatus>& cs, double Delta_l, bool fresh_point) {
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
                                             const std::vector<fluid::CellStatus>& cs, bool fresh_point){
        if (!fresh_point) {
            double phi_GP{0};
            for (int i{0}; i < 4; i++) {
                if (cs[i] == fluid::CellStatus::Fluid) {
                    phi_GP += alpha_neu[i] * phi[i];
                }
            }
            return phi_GP;
        }else{

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
        for (int  i{0};i<n_bound;i++) cout<< "force b_i = "<<F_boundary[i]<<endl;
    }

    void SolidBody::reset_containers(){
        //Call before every timestep to clear containers of variable size
        cell2intercept.clear();
        cell2intercept_FP.clear();
        for (int i{0}; i < n_bound; i++){
            segments[i].segment_data_vec.clear();
            segments[i].n_is_set = false;
            F_boundary[i] = {0,0};
        }
    }

    SolidBody::~SolidBody() {
        delete[] boundary;
        delete[] F_boundary;
        delete[] segments;
    }

    DynamicRigid::DynamicRigid(fluid::FVM_Solver &fvm, std::vector<Point>&& boundary_in,Point CM, double M, double I)
    : SolidBody(fvm,std::move(boundary_in),SolidBodyType::Dynamic), M{M}, I{I}{
        y[0] = CM.x;
        y[1] = CM.y;
        y[2] = 0;
        y[3] = 0;
        y[4] = 0;
        y[5] = 0;
        F_solid = {0,0}; //Change later
        r0 = new Point[n_bound];
        for (int i{0};i<n_bound;i++){
            r0[i] = boundary[i] - CM;
        }
    }

    double DynamicRigid::max_boundary_speed() const{
        double tmp;
        double v_norm_max{0};
        double omega = y[5];
        for (int i{0}; i<n_bound;i++){
            Point r = {boundary[i].x - y[0], boundary[i].y - y[1]};
            //v = v_CM + omega x r
            tmp = sqrt(squared(y[2] - omega*r.y) + squared(y[3] + omega*r.y));
            v_norm_max = std::max(v_norm_max,tmp);
        }
        return v_norm_max;
    }

    void DynamicRigid::step_solid_body(double dt) {
        update_total_fluid_force_and_moment();
        RK4_step(dt);
        update_boundary();
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
        cout << "F_fluid: "<<F_fluid<<", tau_fluid"<<tau_fluid<<endl;
        //F_fluid={-1000000,0};
        //tau_fluid=0;
    }


    Vector6d DynamicRigid::evaluate_f(Vector6d y_in){
        //rhs of the state vector derivative dy/dt = f = [u_CM, v_CM, Fx/M, Fy/M, omega, tau/I]^T
        //std::pair<Point,double> F_S_tau_S{eval_solid_forces_and_moment()};
        //cout << "y= "<<y<<endl<<"M="<<M<<", I= "<<I<<endl;
        f[0] = y_in[2];
        f[1] = y_in[3];
        f[2] = (F_fluid.x + F_solid.x)/M;
        f[3] = (F_fluid.y + F_solid.y)/M;
        f[4] = y_in[5];
        f[5] = (tau_fluid + tau_solid)/I;
        //cout << "f="<<f<<endl;
        return f;
    }


    void DynamicRigid::RK4_step(double dt){
        k1 = dt * evaluate_f(y);
        k2 = dt * evaluate_f(y + k1/2);
        k3 = dt * evaluate_f(y + k2/2);
        k4 = dt * evaluate_f(y + k3);
        y += (k1 + 2*(k2+k3) + k4)/6;
    }

    void DynamicRigid::update_boundary(){
        double c{cos(y[4])};
        double s{sin(y[4])};
        for (int i{0};i<n_bound;i++){
            //boundary_i(t) = CM(t) + R(t)*r0
            boundary[i].x = y[0] + c*r0[i].x - s*r0[i].y;
            boundary[i].y = y[1] + s*r0[i].x + c*r0[i].y;
        }
    }

    DynamicRigid::~DynamicRigid(){
        delete[] r0;
    }




}
