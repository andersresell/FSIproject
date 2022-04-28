//
// Created by anders on 2/24/22.
//
#include "simulate.hpp"


Simulate::Simulate(const std::string& input_file) {
    using namespace std;

    namespace pt = boost::property_tree;
    pt::ptree root;
    pt::read_json("input_files/" + input_file, root);

    auto output_folder = root.get<std::string>("output_folder");
    auto ni = root.get<int>("ni");
    auto nj = root.get<int>("nj");
    auto L_x = root.get<double>("L_x");
    auto L_y = root.get<double>("L_y");
    auto CFL = root.get<double>("CFL");

    fluid::FluxScheme flux_scheme;
    if (root.get<std::string>("flux_scheme") == "HLLC"){
        flux_scheme = fluid::FluxScheme::HLLC;
    }else if (root.get<std::string>("flux_scheme") == "Rusanov"){
        flux_scheme = fluid::FluxScheme::Rusanov;
    }else{
        std::cerr << "Invalid flux scheme in " + input_file + '\n';
        exit(1);
    }

    fluid::OdeScheme ode_scheme;
    if (root.get<std::string>("ode_scheme") == "ExplicitEuler"){
        ode_scheme = fluid::OdeScheme::ExplicitEuler;
    } else if (root.get<std::string>("ode_scheme") == "TVD_RK3"){
        ode_scheme = fluid::OdeScheme::TVD_RK3;
    } else{
        std::cerr << "Invalid ode scheme in " + input_file + '\n';
        exit(1);
    }
    auto fvm_write_stride = root.get<int>("fvm_write_stride");

    fluid::BC_Type west_bc, east_bc, south_bc, north_bc;
    west_bc = bc_type_from_str(std::move(root.get<std::string>("west_bc.type")));
    double M_inf;
    if (west_bc == fluid::BC_Type::SupersonicInflow) M_inf = root.get<double>("west_bc.M_inf");
    east_bc = bc_type_from_str(std::move(root.get<std::string>("east_bc.type")));
    south_bc = bc_type_from_str(std::move(root.get<std::string>("south_bc.type")));
    north_bc = bc_type_from_str(std::move(root.get<std::string>("north_bc.type")));
    fluid::ExternalBCs bcs{ni,nj,west_bc,east_bc,south_bc,north_bc,M_inf};

    fluid::FVM_Solver fvm{ni,nj,L_x,L_y,CFL,ode_scheme,flux_scheme,bcs,output_folder};

    auto initial_cond = root.get<std::string>("initial_cond.case");
    if (initial_cond == "riemann_problem") {
        fluid::vec4 V_l{}, V_r{};
        V_l.u1 = root.get<double>("initial_cond.rho_l");
        V_l.u2 = root.get<double>("initial_cond.u_l");
        V_l.u3 = root.get<double>("initial_cond.v_l");
        V_l.u4 = root.get<double>("initial_cond.p_l");
        V_r.u1 = root.get<double>("initial_cond.rho_r");
        V_r.u2 = root.get<double>("initial_cond.u_r");
        V_r.u3 = root.get<double>("initial_cond.v_r");
        V_r.u4 = root.get<double>("initial_cond.p_r");
        fluid::set_initial_cond_riemann(fvm.U,ni,nj,V_l,V_r);
    }else if (initial_cond == "constant_horizontal_flow") {
        fluid::set_constant_horizontal_flow_cond(fvm.U,ni,nj,M_inf);
    }else if (initial_cond == "pressure_bubble") {
        fluid::set_initial_cond_pressure_bubble(fvm.U, ni, nj, L_x, L_y);
    }else if (initial_cond == "ambient"){
        fluid::set_initial_cond_ambient(fvm.U,ni,nj);
    }else if (initial_cond == "initial_cond1"){
        fluid::set_initial_cond1(fvm.U,ni,nj);
    }else if (initial_cond == "initial_cond2"){
        fluid::set_initial_cond2(fvm.U,ni,nj);
    }else if (initial_cond == "initial_cond_shock_tube_experiment"){
        fluid::vec4 V_l{0,0,0,0}, V_r{0,0,0,0};
        V_l.u1 = root.get<double>("initial_cond.rho_l");
        V_l.u4 = root.get<double>("initial_cond.p_l");
        V_r.u1 = root.get<double>("initial_cond.rho_r");
        V_r.u4 = root.get<double>("initial_cond.p_r");
        fluid::set_initial_cond_shock_tube_experiment(fvm.U,ni,nj,L_x,V_l,V_r);
    }else if (initial_cond == "constant_data"){
        fluid::vec4 V{};
        V.u1 = root.get<double>("initial_cond.rho");
        V.u2 = root.get<double>("initial_cond.u");
        V.u3 = root.get<double>("initial_cond.v");
        V.u4 = root.get<double>("initial_cond.p");
        fluid::set_initial_cond_constant_data(fvm.U,ni,nj,V);
    }
    else{
        std::cerr << "Invalid initial cond \"" + initial_cond + "\" in " + input_file + '\n';
        exit(1);
    }

    FSI_Solver fsi{fvm, fvm_write_stride, output_folder};

    auto stopping_crit = root.get<std::string>("stopping_crit.type");
    auto stopping_crit_value = root.get<double>("stopping_crit.value");
    if (stopping_crit == "Time") {
        fsi.set_endtime(stopping_crit_value);
    }else if (stopping_crit == "Timesteps") {
        fsi.set_timesteps((int)stopping_crit_value);
    }else if (stopping_crit == "Convergence") {
        fsi.set_convergence(stopping_crit_value);
    }else {
        std::cerr << "Illegal stopping criterion read from " + input_file + '\n';
    }

    for (pt::ptree::value_type &solid : root.get_child("solids")) {
        //Definately better ways to do this, this is mainly because I'm bad at json
        string solid_str = "solids." + solid.first;
        auto is_static = root.get<bool>(solid_str + ".is_static");
        solid::SolidBodyType solid_body_type;
        double rho;
        double M, I;
        solid::Point CM{};
        if (is_static) {
            solid_body_type = solid::SolidBodyType::Static;
        } else {
            solid_body_type = solid::SolidBodyType::Dynamic;
            rho = root.get<double>(solid_str + ".rho");
        }

        string geom_str = solid_str + ".geometry";
        vector<solid::Point> boundary;
        if (root.get<string>(geom_str + ".case") == "circle") {
            auto R = root.get<double>(geom_str + ".R");
            auto n_nodes = root.get<int>(geom_str + ".n_nodes");
            auto x_center = root.get<double>(geom_str + ".x_center");
            auto y_center = root.get<double>(geom_str + ".y_center");
            boundary = solid::generate_circle(R, n_nodes, x_center, y_center);
            if (solid_body_type == solid::SolidBodyType::Dynamic) {
                M = M_PI * sqr(R) * rho;
                I = 0.5 * M * sqr(R);
                CM = {x_center, y_center};
            }
        } else if (root.get<string>(geom_str + ".case") == "wedge") {
            auto l = root.get<double>(geom_str + ".l");
            auto half_angle_deg = root.get<double>(geom_str + ".half_angle_deg");
            auto x_center = root.get<double>(geom_str + ".x_center");
            auto y_center = root.get<double>(geom_str + ".y_center");
            boundary = solid::generate_wedge(l, half_angle_deg, x_center, y_center);
            if (solid_body_type == solid::SolidBodyType::Dynamic) {
                double h = l * sin(half_angle_deg * M_PI / 180);
                M = 2 * h * l * rho;
                I = rho / 3 * h * l *
                    (sqr(h) + sqr(l)); //Calculated by solving a double integral, might be wrong
                CM = {x_center, y_center};
            }
        } else if (root.get<string>(geom_str + ".case") == "rectangle") {
            auto W = root.get<double>(geom_str + ".W");
            auto H = root.get<double>(geom_str + ".H");
            auto rotation_angle_deg = root.get<double>(geom_str + ".rotation_angle_deg");
            auto x_center = root.get<double>(geom_str + ".x_center");
            auto y_center = root.get<double>(geom_str + ".y_center");
            boundary = solid::generate_rectangle(W, H, rotation_angle_deg, x_center, y_center);
            if (solid_body_type == solid::SolidBodyType::Dynamic) {
                M = W * H * rho;
                I = M / 12 * (sqr(W) + sqr(H));
                CM = {x_center, y_center};
            }
        } else {
            std::cerr << "Invalid solid geometry in solid " + solid.first + " from " + input_file + '\n';
            exit(1);
        }
        if (is_static) {
            fsi.add_solid(std::move(std::make_shared<solid::SolidBody>(fvm, std::move(boundary), solid_body_type)));
        } else {
            std::shared_ptr<solid::DynamicRigid> dynamic_rigid{
                    std::make_shared<solid::DynamicRigid>(fvm, std::move(boundary), CM, M, I)};
            auto is_constrained = root.get<bool>(solid_str + ".is_constrained");
            if (is_constrained) {
                solid::RigidConstraints rigid_constraints{};
                string constraints_str = solid_str + ".constraints";
                if (root.get<string>(constraints_str + ".type") == "ViscoElastic") {
                    auto K = root.get<double>(constraints_str + ".K");
                    auto C = root.get<double>(constraints_str + ".C");
                    rigid_constraints.setup_viscoelastic(K,C);
                }else if(root.get<string>(constraints_str + ".type") == "PrescribedVelocity"){
                    auto vel = root.get<double>(constraints_str + ".vel");
                    rigid_constraints.setup_prescribed_velocity(vel);
                }
                dynamic_rigid->add_rigid_constraints(std::move(rigid_constraints));
            }
            fsi.add_solid(std::move(dynamic_rigid));
        }

    }

    std::cout << "Setup complete\n";
    fsi.solve();
}


fluid::BC_Type Simulate::bc_type_from_str(std::string&& bc_type){
    if (bc_type == "InvicidWall"){
        return fluid::BC_Type::InvicidWall;
    } else if (bc_type == "NonreflectingOutflow"){
        return fluid::BC_Type::NonreflectingOutflow;
    } else if (bc_type == "SupersonicInflow"){
        return fluid::BC_Type::SupersonicInflow;
    } else{
        std::cerr << "Invalid wall condition \"" + bc_type + "\" in json file\n";
        exit(1);
    }
}

