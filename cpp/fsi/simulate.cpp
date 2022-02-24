//
// Created by anders on 2/24/22.
//
#include "simulate.hpp"


Simulate::Simulate(const std::string& input_file) {
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

    fluid::FVM_Solver fvm{ni,nj,L_x,L_y,CFL,ode_scheme,flux_scheme,bcs};

    auto initial_cond = root.get<std::string>("initial_cond.case");
    if (initial_cond == "riemann_problem") {
        fluid::vec4 V_l, V_r;
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
    } else{
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