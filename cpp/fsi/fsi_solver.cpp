//
// Created by anders on 12/29/21.
//

#include "fsi_solver.hpp"


FSI_Solver::FSI_Solver(fluid::FVM_Solver& fvm, int fvm_write_stride, std::string fvm_output_folder)
    : fvm{fvm}, fvm_write_stride{fvm_write_stride}, fvm_output_folder{std::move(fvm_output_folder)}
{
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
    bool breaker{false};
    while (true) {
        std::cout << "FSI solve: n = " + std::to_string(n) + "\n";

        /*
        for (auto& sb_ptr : solid_bodies){ //Setting solid body bc's only once per iteration
            sb_ptr->set_bc();
        }*/

        if (n % fvm_write_stride == 0) {
            std::cout << "Writing FVM output\n";
            fvm.write_fvm_csv_out_file(fvm_output_folder, n);
        }

        dt = fvm.ode_step(); //Time stepping the fvm simulation
        std::cout << "dt = " << dt << '\n';

        if (stopping_crit.first == StoppingCrit::Time) {
            if (t >= stopping_crit.second) breaker = true;
        } else if (stopping_crit.first == StoppingCrit::Timesteps) {
            if (n >= (int) stopping_crit.second) breaker = true;
        } else {
            std::cerr << "Invalid stopping criterion\n";
            exit(1);
        }
        //Writing fvm header at last timestep, since the number of steps can't be known a priori if endtime is
        //used as stopping criterion
        if (breaker) {
            fvm.write_fvm_csv_header_file(fvm_output_folder, fvm_write_stride, n, t);
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

void FSI_Solver::solid_test(){
    using namespace std;
    int ni = 300;
    int nj = 300;
    double L_x = 10;
    double L_y = 10;
    double CFL = 0.5;
    int n_timesteps{500};
    int fvm_write_stride{10};
    //std::string fvm_output_folder{"output0"};
    std::string fvm_output_folder{"output1"};
    //fluid::FluxScheme flux_scheme{fluid::FluxScheme::Rusanov};
    //fluid::OdeScheme ode_scheme{fluid::OdeScheme::ExplicitEuler};
    fluid::FluxScheme flux_scheme{fluid::FluxScheme::HLLC};
    fluid::OdeScheme ode_scheme{fluid::OdeScheme::TVD_RK3};
    fluid::BC_Type wall = fluid::BC_Type::InvicidWall;
    //fluid::ExternalBCs bcs{ni,nj,wall,wall,wall,wall};
    double M_inf{3};
    fluid::ExternalBCs bcs{ni,nj,fluid::BC_Type::SupersonicInflow, fluid::BC_Type::NonreflectingOutflow,wall, wall, M_inf};
    fluid::FVM_Solver fvm{ni,nj,L_x,L_y,CFL,ode_scheme,flux_scheme,bcs};
    //fluid::set_initial_cond2(fvm.U, fvm.ni, fvm.nj);
    fluid::set_constant_horizontal_flow_cond(fvm.U,ni,nj,M_inf);
    FSI_Solver fsi{fvm, fvm_write_stride, fvm_output_folder};
    fsi.set_timesteps(n_timesteps);
    {
        vector<solid::Point> circle1;
        vector<solid::Point> circle2;
        vector<solid::Point> circle3;
        double R = 2.99999999;
        double x0 = 5;
        double y0 = 5;
        int n = 4;
        for (int i{0}; i < n; i++) {
            double theta = 2 * M_PI * i / n;
            circle1.push_back(solid::Point{x0 + R * cos(theta), R * sin(theta)});
            circle2.push_back(solid::Point{x0 + R * cos(theta), L_y + R * sin(theta)});
        }
        n = 20;
        R = 0.76;
        for (int i{0}; i < n; i++) {
            double theta = 2 * M_PI * i / n;
            circle3.push_back(solid::Point{x0 + R * cos(theta), y0 + R * sin(theta)});
        }
        vector<solid::Point> tri;
        tri.push_back(solid::Point{4,5});
        tri.push_back(solid::Point{6,4});
        tri.push_back(solid::Point{6,6});

        //fsi.add_solid(std::make_unique<solid::SolidBody>(fvm,circle1));
        //fsi.add_solid(std::make_unique<solid::SolidBody>(fvm,circle2));
        //fsi.add_solid(std::make_unique<solid::SolidBody>(fvm, circle3));
        fsi.add_solid(std::move(std::make_shared<solid::SolidBody>(fvm,tri)));
    }
    fsi.solve();

    if (fsi.solid_bodies.size() > 0) {
        fsi.solid_bodies[0]->debug_intercepts_csv();
        fsi.solid_bodies[0]->debug_csv();
        fsi.solid_bodies[0]->write_boundary_csv("solid_debug_boundary");
    }



}

void FSI_Solver::fluid_solve_test(int ni,
                 int nj,
                 double L_x,
                 double L_y,
                 double CFL,
                 int n_timesteps,
                 int fvm_write_stride,
                 std::string fvm_output_folder,
                 fluid::OdeScheme ode_scheme,
                 fluid::FluxScheme flux_scheme) {
    //Should also add functionality for choosing the time scheme, flux scheme, boundaries and initital cond
    fluid::BC_Type wall = fluid::BC_Type::InvicidWall;
    fluid::ExternalBCs wall_BCs{ni,nj,wall,wall,wall,wall};
    fluid::FVM_Solver fvm{ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, wall_BCs};
    fluid::set_initial_cond1(fvm.U, fvm.ni, fvm.nj);
    FSI_Solver fsi{fvm, fvm_write_stride, fvm_output_folder};
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
                                std::string fvm_output_folder,
                                fluid::OdeScheme ode_scheme,
                                fluid::FluxScheme flux_scheme){
    fluid::BC_Type wall = fluid::BC_Type::InvicidWall;
    fluid::ExternalBCs wall_BCs{ni,nj,wall,wall,wall,wall};
    fluid::FVM_Solver fvm{ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, wall_BCs};
    fluid::set_initial_cond_riemann(fvm.U,fvm.ni,fvm.nj, V_l, V_r);
    FSI_Solver fsi{fvm, fvm_write_stride, fvm_output_folder};
    //fsi.set_timesteps(n_timesteps);
    fsi.set_endtime(endtime);
    fsi.solve();
}