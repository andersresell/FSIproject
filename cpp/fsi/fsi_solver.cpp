//
// Created by anders on 12/29/21.
//

#include "fsi_solver.hpp"


FSI_Solver::FSI_Solver(fluid::FVM_Solver& fvm, int fvm_write_stride, std::string output_folder)
    : fvm{fvm}, fvm_write_stride{fvm_write_stride}, output_folder{std::move(output_folder)}
{
    rho_old = new double[(fvm.ni + 4) * (fvm.nj + 4)];
}

void FSI_Solver::add_solid(std::shared_ptr<solid::SolidBody>&& solid_body) {
    solid_bodies.push_back(solid_body);
    fvm.solid_bodies.push_back(std::move(solid_body));
}

int FSI_Solver::solve() {
    auto start_time{std::chrono::high_resolution_clock::now()}; //Start timing
    int n{0};
    double t{0};
    double dt{0};
    write_static_solid_boundaries();
    bool breaker{false};
    double res_norm, res_norm0;
    double mass_0, mass;
    while (true) {
        std::cout << "FSI solve: n = " + std::to_string(n) + "\n";

        if (n % fvm_write_stride == 0) {
            std::cout << "Writing output\n";
            fvm.write_fvm_output(output_folder, n);
            write_movable_solid_boundaries(n);
            write_solid_debug_files(n);
            set_rho_old();
        }

        dt = fvm.ode_step(dt); //Time stepping the fvm simulation
        std::cout << "dt = " << dt << ", t = " << t+dt << '\n';

        if (n == 0) {
            res_norm0 = calc_density_L2_norm();
            mass_0 = calc_mass();
            totals_history.push_back({res_norm0,mass_0});
        } else if (n % fvm_write_stride == 0) {
            //Checking if ||rho_{n+1} - rho_n|| <= tol * ||rho_1 - rho_0||
            res_norm = calc_density_L2_norm();
            mass = calc_mass();
            totals_history.push_back({res_norm,mass});
            std::cout << "Convergence check: res_norm0 = " << res_norm0 << ", res_norm = " << res_norm << '\n';
            std::cout << "mass_0 = " << mass_0 << ", mass = " << mass << ", mass change = " << mass - mass_0 << '\n';
        }

        if (stopping_crit.first == StoppingCrit::Time) {
            if (t >= stopping_crit.second) breaker = true;
        }
        else if (stopping_crit.first == StoppingCrit::Timesteps) {
            if (n >= (int) stopping_crit.second) breaker = true;
        }
        else if (stopping_crit.first == StoppingCrit::Convergence) {
            if (n % fvm_write_stride == 0 && n > 0) {
                if (res_norm <= stopping_crit.second * res_norm0) {
                    breaker = true;
                }
                //std::cout <<"SC2 = "<< stopping_crit.second << std::endl;
                std::cout << "Desired res norm is " << 100 * stopping_crit.second * res_norm0 / res_norm
                          << "% of current res norm\n";
            }
        } else {
            std::cerr << "Invalid stopping criterion\n";
            exit(1);
        }
        //Writing fvm header at last timestep, since the number of steps can't be known a priori if endtime is
        //used as stopping criterion
        if (breaker) {
            if (n % fvm_write_stride != 0) fvm.write_fvm_output(output_folder, n); //making sure that last step is read
            fvm.write_fvm_header(output_folder, fvm_write_stride, n, t);
            write_fsi_header();
            write_totals_history();
            write_movable_solid_boundaries(n);
            write_solid_debug_files(n);
            break;
        }
        t += dt;
        n++;
        std::cout << '\n';
    }
    auto stop_time{std::chrono::high_resolution_clock::now()};
    auto simulation_time{std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time)};
    std::cout << "End time = " << t << '\n';
    std::cout << "Computation time = " << 1.0e-3 * static_cast<double>(simulation_time.count()) << " seconds\n";
    return n;
}


double FSI_Solver::calc_density_L2_norm(){
    int nj = fvm.nj; //needed for the IX macro
    double res_norm{0};
    for (int i{2}; i < fvm.ni + 2; i++) {
        for (int j{2}; j < nj + 4; j++) {
            if (fvm.cell_status[IX(i,j)] == fluid::CellStatus::Fluid) res_norm += squared(fvm.U[IX(i,j)].u1 - rho_old[IX(i,j)]);
        }
    }
    return sqrt(fvm.dx*fvm.dy*res_norm);
}
double FSI_Solver::calc_mass(){
    int nj = fvm.nj;
    double tmp{0};
    for (int i{2}; i < fvm.ni + 2; i++) {
        for (int j{2}; j < nj + 4; j++) {
            if (fvm.cell_status[IX(i,j)] == fluid::CellStatus::Fluid) tmp += fvm.U[IX(i,j)].u1;
        }
    }
    return tmp*fvm.dx*fvm.dy;
}


void FSI_Solver::set_rho_old() {
    int nj = fvm.nj; //needed for the IX macro
    for (int i{2}; i < fvm.ni + 2; i++) {
        for (int j{2}; j < nj + 4; j++) {
            rho_old[IX(i, j)] = fvm.U[IX(i,j)].u1;
        }
    }
}

void FSI_Solver::write_totals_history() {
    std::ofstream ost{"python/output_folders/" + output_folder + "/fvm_convergence_history.csv"};
    if (!ost) {
        std::cerr << "Error: couldn't open fvm convergence csv output file\n";
        exit(1);
    }
    ost << "#n,norm,mass\n";
    for (int i{0}; i < totals_history.size(); i++) {
        ost << i * fvm_write_stride << ',' << totals_history[i].density_L2_norm << ',' << totals_history[i].mass << '\n';
    }
}


void FSI_Solver::write_fsi_header() {
    std::ofstream ost{"python/output_folders/" + output_folder + "/fsi_header.csv"};
    if (!ost) {
        std::cerr << "Error: couldn't open fsi csv header file\n";
        exit(1);
    }
    ost << "#n_static_solids,n_movable_solids\n";
    int n_static{0};
    int n_movable{0};
    for (auto &s: solid_bodies) {
        if (s->type == solid::SolidBodyType::Static) {
            n_static++;
        } else if (s->type == solid::SolidBodyType::Dynamic) {
            n_movable++;
        } else {
            std::cerr << "Error: Solid Body type is neither Static or Movable (in write fsi header)\n";
            exit(1);
        }
    }
    ost << n_static << ',' << n_movable << '\n';
}

void FSI_Solver::write_static_solid_boundaries() {
    int solid_ind{0};
    for (auto &s: solid_bodies) {
        if (s->type == solid::SolidBodyType::Static) {
            std::ofstream ost{
                    "python/output_folders/" + output_folder + "/static_boundary" + std::to_string(solid_ind) +
                    ".csv"};
            if (!ost) {
                std::cerr << "error: couldn't open static boundary csv file\n";
                exit(1);
            }
            ost << "#x,y\n";
            for (int i{0}; i < s->n_bound; i++) {
                ost << s->boundary[i].x << ',' << s->boundary[i].y << '\n';
            }
            solid_ind++;
        }
    }
}

void FSI_Solver::write_movable_solid_boundaries(int n) {
    int solid_ind{0};
    for (auto &s: solid_bodies) {
        if (s->type == solid::SolidBodyType::Dynamic) {
            std::ofstream ost{
                    "python/output_folders/" + output_folder + "/movable_boundary" + std::to_string(solid_ind) +
                    "_t" + std::to_string(n) + ".csv"};
            if (!ost) {
                std::cerr << "error: couldn't open movable boundary csv file\n";
                exit(1);
            }
            ost << "#x,y\n";
            for (int i{0}; i < s->n_bound; i++) {
                ost << s->boundary[i].x << ',' << s->boundary[i].y << '\n';
            }
            solid_ind++;
        }
    }
}
void FSI_Solver::write_solid_debug_files(int n){
    //Writes the status of the solid bodies. One file for the node cell status and one for the intercepts
    int nj = fvm.nj;
    std::ofstream ost1{"python/output_folders/" + output_folder + "/debug_nodes_t" + std::to_string(n) + ".csv"};
    if (!ost1) std::cerr << "error, couldn't open node debug csv file\n";
    ost1 << "#type,x,y\n";
    for (int i{0}; i < fvm.ni + 4; i++) {
        for (int j{0}; j < nj + 4; j++) {
            solid::Point p = {(i - 1.5) * fvm.dx, (j - 1.5) * fvm.dy};
            ost1 << static_cast<int>(fvm.cell_status[IX(i, j)]) << ',' << p.x << ',' << p.y << '\n';
        }
    }

    std::ofstream ost2{"python/output_folders/"+output_folder+"/debug_intercepts_t" + std::to_string(n) + ".csv"};
    if (!ost2) std::cerr << "error, couldn't open solid debug intercepts csv file\n";
    ost2 << "#x_i,y_i\n";
    for (auto& s: solid_bodies) {
        for(auto&e : s->cell2intercept){
            ost2 << e.second.BI.x << ',' <<e.second.BI.y << '\n';
        }
    }



}

FSI_Solver::~FSI_Solver(){
    delete[] rho_old;
}
