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
    double dt;
    fvm.initialize_solids();
    write_static_solid_boundaries();
    bool breaker{false};
    double res_norm, res_norm0;
    while (true) {
        std::cout << "FSI solve: n = " + std::to_string(n) + "\n";

        if (n % fvm_write_stride == 0) {
            std::cout << "Writing output\n";
            fvm.write_fvm_output(output_folder, n);
            write_movable_solid_boudnaries(n);
            set_rho_old();
        }

        dt = fvm.ode_step(); //Time stepping the fvm simulation
        std::cout << "dt = " << dt << '\n';

        if (n == 0) {
            res_norm0 = calc_density_L2_norm();
            convergence_history.push_back(res_norm0);
        } else if (n % fvm_write_stride == 0) {
            //Checking if ||rho_{n+1} - rho_n|| <= tol * ||rho_1 - rho_0||
            res_norm = calc_density_L2_norm();
            convergence_history.push_back(res_norm);
            std::cout << "Convergence check: res_norm0 = " << res_norm0 << ", res_norm = " << res_norm << '\n';
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
            write_fvm_convergence_history();
            break;
        }
        t += dt;
        n++;
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
            res_norm += squared(fvm.U[IX(i,j)].u1 - rho_old[IX(i,j)]);
        }
    }
    return sqrt(fvm.dx*fvm.dy*res_norm);
}

void FSI_Solver::set_rho_old() {
    int nj = fvm.nj; //needed for the IX macro
    for (int i{2}; i < fvm.ni + 2; i++) {
        for (int j{2}; j < nj + 4; j++) {
            rho_old[IX(i, j)] = fvm.U[IX(i,j)].u1;
        }
    }
}

void FSI_Solver::write_fvm_convergence_history() {
    std::ofstream ost{"python/output_folders/" + output_folder + "/fvm_convergence_history.csv"};
    if (!ost) {
        std::cerr << "Error: couldn't open fvm convergence csv output file\n";
        exit(1);
    }
    ost << "#n,norm\n";
    for (int i{0}; i < convergence_history.size(); i++) {
        ost << i * fvm_write_stride << ',' << convergence_history[i] << '\n';
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
        } else if (s->type == solid::SolidBodyType::Movable) {
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
        }
    }
}

void FSI_Solver::write_movable_solid_boudnaries(int n) {
    int solid_ind{0};
    for (auto &s: solid_bodies) {
        if (s->type == solid::SolidBodyType::Movable) {
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
        }
    }
}

void FSI_Solver::riemann_test(){
    int ni{100};
    int nj{100};
    double L_x = 4;
    double L_y = 4;
    double CFL = 0.8;
    int n_timesteps{100};
    int fvm_write_stride{1};
    std::string output_folder{"output_riemann_sod"};
    //fluid::OdeScheme ode_scheme{fluid::OdeScheme::TVD_RK3};
    fluid::FluxScheme flux_scheme{fluid::FluxScheme::HLLC};
    fluid::OdeScheme ode_scheme{fluid::OdeScheme::ExplicitEuler};
    //fluid::FluxScheme flux_scheme{fluid::FluxScheme::Rusanov};

    //Testing sod's problem
    double rho_l = 3;
    double u_l = 0;
    double v_l = 0;
    double p_l = 3;
    double rho_r = 1;
    double u_r = 0;
    double v_r = 0;
    double p_r = 1;

    fluid::vec4 V_l{rho_l, u_l, v_l, p_l};
    fluid::vec4 V_r{rho_r, u_r, v_r, p_r};
    double endtime{0.5};
    fluid::BC_Type wall = fluid::BC_Type::InvicidWall;
    fluid::ExternalBCs bcs{ni,nj,wall,wall,wall,wall};

    fluid::FVM_Solver fvm{ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, bcs};
    fluid::set_initial_cond_riemann(fvm.U,fvm.ni,fvm.nj, V_l, V_r);
    FSI_Solver fsi{fvm, fvm_write_stride, output_folder};
    //fsi.set_timesteps(n_timesteps);
    fsi.set_endtime(endtime);
    fsi.solve();



}

void FSI_Solver::solid_test(){
    using namespace std;
    int ni = 250;
    int nj = 250;
    double L_x = 10;
    double L_y = 10;
    double CFL = 0.5;
    int n_timesteps{500};
    int fvm_write_stride{10};
    //std::string output_folder{"output0"};
    std::string output_folder{"output1"};
    fluid::FluxScheme flux_scheme{fluid::FluxScheme::Rusanov};
    //fluid::OdeScheme ode_scheme{fluid::OdeScheme::ExplicitEuler};
    //fluid::FluxScheme flux_scheme{fluid::FluxScheme::HLLC};
    fluid::OdeScheme ode_scheme{fluid::OdeScheme::TVD_RK3};
    fluid::BC_Type wall = fluid::BC_Type::InvicidWall;
    fluid::ExternalBCs bcs{ni,nj,wall,wall,wall,wall};
    double M_inf{3};
    //fluid::ExternalBCs bcs{ni,nj,fluid::BC_Type::SupersonicInflow, fluid::BC_Type::NonreflectingOutflow,wall, wall, M_inf};
    fluid::FVM_Solver fvm{ni,nj,L_x,L_y,CFL,ode_scheme,flux_scheme,bcs};
    fluid::set_initial_cond2(fvm.U, fvm.ni, fvm.nj);
    //fluid::set_constant_horizontal_flow_cond(fvm.U,ni,nj,M_inf);
    FSI_Solver fsi{fvm, fvm_write_stride, output_folder};
    fsi.set_timesteps(n_timesteps);
    {
        vector<solid::Point> circle1;
        vector<solid::Point> circle2;
        vector<solid::Point> circle3;
        double R = 1.99999;
        double x0 = 5;
        double y0 = 5;
        int n = 4;
        for (int i{0}; i < n; i++) {
            double theta = 2 * M_PI * i / n;
            circle1.push_back(solid::Point{x0 + R * cos(theta), R * sin(theta)});
            circle2.push_back(solid::Point{x0 + R * cos(theta), L_y + R * sin(theta)});
        }
        n = 10;
        R = 0.76;
        for (int i{0}; i < n; i++) {
            double theta = 2 * M_PI * i / n;
            circle3.push_back(solid::Point{x0 + R * cos(theta), y0 + R * sin(theta)});
        }
        vector<solid::Point> tri;
        tri.push_back(solid::Point{4,5});
        tri.push_back(solid::Point{6,4});
        tri.push_back(solid::Point{6,6});

        fsi.add_solid(std::make_unique<solid::SolidBody>(fvm,circle1,solid::SolidBodyType::Static));
        fsi.add_solid(std::make_unique<solid::SolidBody>(fvm,circle2,solid::SolidBodyType::Static));
        fsi.add_solid(std::make_unique<solid::SolidBody>(fvm, circle3,solid::SolidBodyType::Static));
        //fsi.add_solid(std::move(std::make_shared<solid::SolidBody>(fvm,tri,solid::SolidBodyType::Static)));
    }
    fsi.solve();

    if (fsi.solid_bodies.size() > 0) {
        fsi.solid_bodies[0]->debug_intercepts_csv();
        fsi.solid_bodies[0]->debug_csv();
        //fsi.solid_bodies[0]->write_boundary_csv("solid_debug_boundary");
    }



}

void FSI_Solver::fluid_solve_test(int ni,
                 int nj,
                 double L_x,
                 double L_y,
                 double CFL,
                 int n_timesteps,
                 int fvm_write_stride,
                 std::string output_folder,
                 fluid::OdeScheme ode_scheme,
                 fluid::FluxScheme flux_scheme) {
    //Should also add functionality for choosing the time scheme, flux scheme, boundaries and initital cond
    fluid::BC_Type wall = fluid::BC_Type::InvicidWall;
    fluid::ExternalBCs wall_BCs{ni,nj,wall,wall,wall,wall};
    fluid::FVM_Solver fvm{ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, wall_BCs};
    fluid::set_initial_cond1(fvm.U, fvm.ni, fvm.nj);
    FSI_Solver fsi{fvm, fvm_write_stride, output_folder};
    fsi.set_timesteps(n_timesteps);
    fsi.solve();
}

void FSI_Solver::fluid_solve_riemann(int ni,
                                int nj,
                                double L_x,
                                double L_y,
                                double CFL,
                                const fluid::vec4& V_l,
                                const fluid::vec4& V_r,
                                double endtime,
                                int fvm_write_stride,
                                std::string output_folder,
                                fluid::OdeScheme ode_scheme,
                                fluid::FluxScheme flux_scheme){
    fluid::BC_Type wall = fluid::BC_Type::InvicidWall;
    fluid::ExternalBCs wall_BCs{ni,nj,wall,wall,wall,wall};
    fluid::FVM_Solver fvm{ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, wall_BCs};
    fluid::set_initial_cond_riemann(fvm.U,fvm.ni,fvm.nj, V_l, V_r);
    FSI_Solver fsi{fvm, fvm_write_stride, output_folder};
    //fsi.set_timesteps(n_timesteps);
    fsi.set_endtime(endtime);
    fsi.solve();
}

void FSI_Solver::wedge_verification(){
    using namespace std;
    int ni = 250;
    int nj = 250;
    double L_x = 10;
    double L_y = 10;
    double CFL = 0.5;
    int n_timesteps{10000};
    int fvm_write_stride{100};

    //std::string output_folder{"wedge_verification_attached"};
    std::string output_folder{"wedge_verification_detached"};
    fluid::FluxScheme flux_scheme{fluid::FluxScheme::Rusanov};
    //fluid::OdeScheme ode_scheme{fluid::OdeScheme::ExplicitEuler};
    //fluid::FluxScheme flux_scheme{fluid::FluxScheme::HLLC};
    fluid::OdeScheme ode_scheme{fluid::OdeScheme::TVD_RK3};
    fluid::BC_Type wall = fluid::BC_Type::InvicidWall;
    fluid::BC_Type ss_outflow = fluid::BC_Type::NonreflectingOutflow;
    double M_inf{3};
    fluid::ExternalBCs bcs{ni,nj,fluid::BC_Type::SupersonicInflow,ss_outflow,ss_outflow,ss_outflow,M_inf};

    fluid::FVM_Solver fvm{ni,nj,L_x,L_y,CFL,ode_scheme,flux_scheme,bcs};
    fluid::set_constant_horizontal_flow_cond(fvm.U,ni,nj,M_inf);
    FSI_Solver fsi{fvm, fvm_write_stride, output_folder};
    fsi.set_timesteps(n_timesteps);
    //fsi.set_convergence();
    {
        vector<solid::Point> circle1;
        vector<solid::Point> circle2;
        vector<solid::Point> circle3;

        //Constructing a wedge around origin
        vector<solid::Point> wedge;
        double theta = 35 * M_PI / 180; //half angle of the wedge
        double l = 3;
        double h = l * tan(theta);
        wedge.push_back(solid::Point{L_x / 2 - l, L_y / 2});
        wedge.push_back(solid::Point{L_x / 2, L_y / 2 - h});
        wedge.push_back(solid::Point{L_x / 2 + l, L_y / 2});
        wedge.push_back(solid::Point{L_x / 2, L_y / 2 + h});

        fsi.add_solid(std::make_unique<solid::SolidBody>(fvm, wedge, solid::SolidBodyType::Static));
    }
    fsi.solve();

    if (~fsi.solid_bodies.empty()) {
        fsi.solid_bodies[0]->debug_intercepts_csv();
        fsi.solid_bodies[0]->debug_csv();
        //fsi.solid_bodies[0]->write_boundary_csv("solid_debug_boundary",0,0);
    }
}

FSI_Solver::~FSI_Solver(){
    delete[] rho_old;
}