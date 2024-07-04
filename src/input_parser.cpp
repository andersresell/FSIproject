#include "input_parser.hpp"
#include "fluid/fvm_solver.hpp"

InputParser::InputParser(const string &input_filename) : input_filename{input_filename} {
    using namespace std;
    try {
        root_node = YAML::LoadFile(input_filename);
    } catch (const exception &e) {
        if (input_filename.size() > 4) {
            string yml_extension = input_filename.substr(input_filename.size() - 4);
            if (yml_extension != ".yml") {
                throw runtime_error("Error loading input file '" + input_filename + "'\n" +
                                    "Don't forget .yml extension in name");
            }
        }

        throw runtime_error("Error loading input file '" + input_filename + "', " + string(e.what()) +
                            "\nNote that the input file has to be located in the configurations directory.\n");
    }
    cout << "Input file: " << input_filename << endl;
    base_dir = filesystem::current_path();
    if (base_dir[base_dir.size() - 1] != '/')
        base_dir += '/';
    cout << "Simulation directory: \n" << base_dir << endl;

    // Create output directory
    output_dir = base_dir + "output/";
    if (filesystem::exists(output_dir)) {
        if (!filesystem::remove_all(output_dir)) {
            throw runtime_error{"Failed to remove old output directory: " + output_dir};
        }
    }
    if (!filesystem::create_directory(output_dir)) {
        throw runtime_error{string("Failed to create output directory:" + output_dir + "\n")};
    }
}

void InputParser::create_solvers(unique_ptr<fluid::FVM_Solver> &fvm, unique_ptr<FSI_Solver> &fsi) {
    const string root_name_setup = "setup";
    const string root_name_ext_bcs = "external_bcs";
    const string root_name_initial_cond = "initial_cond";
    const string root_name_solids = "solids";
    try {

        /*--------------------------------------------------------------------
        Setup
        --------------------------------------------------------------------*/
        const size_t num_threads = read_optional_option<size_t>(root_name_setup, "num_threads", 1);
        if (num_threads > 8) {
            throw runtime_error("Specify < 8 threads for open mp\n");
        }
        omp_set_num_threads(num_threads);
        const int ni = read_required_option<int>(root_name_setup, "nx");
        const int nj = read_required_option<int>(root_name_setup, "ny");
        const double L_x = read_required_option<double>(root_name_setup, "L_x");
        const double L_y = read_required_option<double>(root_name_setup, "L_y");
        const double CFL = read_required_option<double>(root_name_setup, "CFL");
        const fluid::FluxScheme flux_scheme = read_optional_enum_option<fluid::FluxScheme>(
            root_name_setup, "flux_scheme", flux_scheme_from_string, fluid::FluxScheme::HLLC);
        const fluid::OdeScheme ode_scheme = read_optional_enum_option<fluid::OdeScheme>(
            root_name_setup, "ode_scheme", ode_scheme_from_string, fluid::OdeScheme::TVD_RK3);
        const fluid::Limiter limiter = read_optional_enum_option<fluid::Limiter>(
            root_name_setup, "limiter", limiter_from_string, fluid::Limiter::MC);
        const int fvm_write_stride = read_required_option<int>(root_name_setup, "fvm_write_stride");

        const StoppingCrit stopping_crit_type =
            read_required_enum_option<StoppingCrit>(root_name_setup, "stopping_crit_type", stopping_crit_from_string);
        const double stopping_crit_value = read_required_option<double>(root_name_setup, "stopping_crit_val");

        /*--------------------------------------------------------------------
        External boundary conditions
        --------------------------------------------------------------------*/
        const fluid::BC_Type west_bc =
            read_required_enum_option<fluid::BC_Type>(root_name_ext_bcs, "type_west", bc_type_from_string);
        double M_inf, p_inf, rho_inf;
        if (west_bc == fluid::BC_Type::SupersonicInflow) {
            M_inf = read_required_option<double>(root_name_ext_bcs, "M_inf",
                                                 "Freestream mach number must be specified for supersonic inflow");

            p_inf = read_required_option<double>(root_name_ext_bcs, "p_inf",
                                                 "Freestream pressure must be specified for supersonic inflow");

            rho_inf = read_required_option<double>(root_name_ext_bcs, "rho_inf",
                                                   "Freestream density must be specified for supersonic inflow");
        }
        const fluid::BC_Type east_bc =
            read_required_enum_option<fluid::BC_Type>(root_name_ext_bcs, "type_east", bc_type_from_string);

        const fluid::BC_Type south_bc =
            read_required_enum_option<fluid::BC_Type>(root_name_ext_bcs, "type_south", bc_type_from_string);

        const fluid::BC_Type north_bc =
            read_required_enum_option<fluid::BC_Type>(root_name_ext_bcs, "type_north", bc_type_from_string);

        if (west_bc == fluid::BC_Type::TimeHistory || east_bc == fluid::BC_Type::TimeHistory ||
            south_bc == fluid::BC_Type::TimeHistory || north_bc == fluid::BC_Type::TimeHistory) {
            throw runtime_error("Time history boundary condition type is currently not supported");
        }

        fluid::ExternalBCs bcs = {ni, nj, west_bc, east_bc, south_bc, north_bc, M_inf, p_inf, rho_inf, ""};

        fvm = make_unique<fluid::FVM_Solver>(ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, limiter, bcs, output_dir);

        /*--------------------------------------------------------------------
        Initial condition
        --------------------------------------------------------------------*/

        const string initial_cond = read_required_option<string>(root_name_initial_cond, "case");
        if (initial_cond == "riemann_problem") {
            fluid::vec4 V_l, V_r;
            V_l.u1 = read_required_option<double>(root_name_initial_cond, "rho_l");
            V_l.u2 = read_required_option<double>(root_name_initial_cond, "u_l");
            V_l.u3 = read_required_option<double>(root_name_initial_cond, "v_l");
            V_l.u4 = read_required_option<double>(root_name_initial_cond, "p_l");
            V_r.u1 = read_required_option<double>(root_name_initial_cond, "rho_r");
            V_r.u2 = read_required_option<double>(root_name_initial_cond, "u_r");
            V_r.u3 = read_required_option<double>(root_name_initial_cond, "v_r");
            V_r.u4 = read_required_option<double>(root_name_initial_cond, "p_r");
            fluid::set_initial_cond_riemann(fvm->U, ni, nj, V_l, V_r);
        }

        else if (initial_cond == "pressure_bubble") {
            double x_c = read_required_option<double>(root_name_initial_cond, "initial_cond.x_c");
            double y_c = read_required_option<double>(root_name_initial_cond, "initial_cond.y_c");
            double radius = read_required_option<double>(root_name_initial_cond, "initial_cond.radius");
            fluid::set_initial_cond_pressure_bubble(fvm->U, ni, nj, L_x, L_y, x_c, y_c, radius);
        }

        else if (initial_cond == "shock_tube_experiment") {
            double driver_length = read_required_option<double>(root_name_initial_cond, "driver_length");
            fluid::vec4 V_l{0, 0, 0, 0}, V_r{0, 0, 0, 0};
            V_l.u1 = read_required_option<double>(root_name_initial_cond, "rho_l");
            V_l.u4 = read_required_option<double>(root_name_initial_cond, "p_l");
            V_r.u1 = read_required_option<double>(root_name_initial_cond, "rho_r");
            V_r.u4 = read_required_option<double>(root_name_initial_cond, "p_r");
            fluid::set_initial_cond_shock_tube_experiment(fvm->U, ni, nj, L_x, driver_length, V_l, V_r);
        } else if (initial_cond == "constant_data") {
            fluid::vec4 V{};
            V.u1 = read_required_option<double>(root_name_initial_cond, "rho");
            V.u2 = read_required_option<double>(root_name_initial_cond, "u");
            V.u3 = read_required_option<double>(root_name_initial_cond, "v");
            V.u4 = read_required_option<double>(root_name_initial_cond, "p");
            fluid::set_initial_cond_constant_data(fvm->U, ni, nj, V);
        } else {
            const string available_cases = "riemann_problem, pressure_bubble,  shock_tube_experiment, constant_data";
            throw runtime_error("Invalid initial condition case, available cases are \'" + available_cases + "\'");
        }

        /*--------------------------------------------------------------------
        Construct FSI solver object
        --------------------------------------------------------------------*/
        fsi = make_unique<FSI_Solver>(*fvm, fvm_write_stride, output_dir);
        switch (stopping_crit_type) {
        case StoppingCrit::Time:
            fsi->set_endtime(stopping_crit_value);
            break;
        case StoppingCrit::Timesteps:
            fsi->set_timesteps((int)stopping_crit_value);
            break;
        case StoppingCrit::Convergence:
            fsi->set_convergence(stopping_crit_value);
            break;
        default:
            assert(false);
        }
        /*--------------------------------------------------------------------
        Construct solids
        --------------------------------------------------------------------*/

        YAML::Node solid_nodes = root_node[root_name_solids];
        if (solid_nodes) {
            int counter = 0;
            for (const YAML::Node &solid_node : solid_nodes) {
                counter++;
                try {
                    double M{0}, rho{0}, I{0};
                    solid::Point CM;
                    vector<solid::Point> boundary;
                    solid::SolidBodyType type = solid::SolidBodyType::Static;
                    if (solid_node["is_static"] && solid_node["is_static"].as<bool>() == false) {
                        type = solid::SolidBodyType::Dynamic;
                        rho = read_option<double>(solid_node, "rho");
                    }

                    YAML::Node geometry_node = solid_node["geometry"];
                    if (geometry_node) {
                        string geom_case;
                        if (geometry_node["case"]) {
                            geom_case = geometry_node["case"].as<string>();
                            // consider placing try-catch blocks between each if to give more context for errors
                            if (geom_case == "circle") {
                                auto radius = read_option<double>(geometry_node, "radius");
                                auto n_nodes = read_option<double>(geometry_node, "n_nodes");
                                auto x_center = read_option<double>(geometry_node, "x_center");
                                auto y_center = read_option<double>(geometry_node, "y_center");
                                boundary = solid::generate_circle(radius, n_nodes, x_center, y_center);
                            } else if (geom_case == "wedge") {
                                auto l = read_option<double>(geometry_node, "l");
                                auto half_angle_deg = read_option<double>(geometry_node, "half_angle_deg");
                                auto x_tip = read_option<double>(geometry_node, "x_tip");
                                auto y_tip = read_option<double>(geometry_node, "y_tip");
                                boundary = solid::generate_wedge(l, half_angle_deg, x_tip, y_tip);
                                if (type == solid::SolidBodyType::Dynamic) {
                                    throw runtime_error("Dynamic properties not yet set for wedge geometry!\n");
                                }
                            } else if (geom_case == "diamond_wedge") {
                                auto l = read_option<double>(geometry_node, "l");
                                auto half_angle_deg = read_option<double>(geometry_node, "half_angle_deg");
                                auto x_center = read_option<double>(geometry_node, "x_center");
                                auto y_center = read_option<double>(geometry_node, "y_center");
                                boundary = solid::generate_diamond_wedge(l, half_angle_deg, x_center, y_center);
                                if (type == solid::SolidBodyType::Dynamic) {
                                    double h = l * sin(half_angle_deg * M_PI / 180);
                                    M = 2 * h * l * rho;
                                    I = rho / 3 * h * l *
                                        (sqr(h) + sqr(l)); // Calculated by solving a double integral, might be wrong
                                    CM = {x_center, y_center};
                                }
                            } else if (geom_case == "rectangle") {
                                auto W = read_option<double>(geometry_node, "width");
                                auto H = read_option<double>(geometry_node, "height");
                                auto rotation_angle_deg = read_option<double>(geometry_node, "rotation_angle_deg");
                                auto x_center = read_option<double>(geometry_node, "x_center");
                                auto y_center = read_option<double>(geometry_node, "y_center");
                                boundary = solid::generate_rectangle(W, H, rotation_angle_deg, x_center, y_center);
                                if (type == solid::SolidBodyType::Dynamic) {
                                    M = W * H * rho;
                                    I = M / 12 * (sqr(W) + sqr(H));
                                    CM = {x_center, y_center};
                                }
                            } else if (geom_case == "random_polygon") {
                                YAML::Node points_node = geometry_node["points"];
                                if (points_node) {
                                    try {
                                        for (const auto &p_node : points_node) {
                                            solid::Point p;
                                            p.x = p_node[0].as<double>();
                                            p.y = p_node[1].as<double>();
                                            boundary.push_back(p);
                                        }
                                    } catch (exception &e) {
                                        throw runtime_error("Failed reading list of points. Specify on the format\n "
                                                            "points:\n  - [px1, py1]\n  - [px2, py2]\n  - ...");
                                    }
                                } else {
                                    throw runtime_error("\'points\' must be specified for a \'random_polygon\'");
                                }
                                if (type == solid::SolidBodyType::Dynamic) {
                                    throw runtime_error("Dynamic solid not available for random polygon geometry");
                                }
                            } else {
                                throw runtime_error("Invalid geometry case");
                            }
                        }
                    } else {
                        throw runtime_error("option \'geometry\' not specified ");
                    }

                    shared_ptr<solid::SolidBody> solid;
                    if (type == solid::SolidBodyType::Static) {
                        solid = make_shared<solid::SolidBody>(*fvm, move(boundary), type);
                    } else {
                        solid = make_shared<solid::DynamicRigid>(*fvm, move(boundary), CM, M, I);
                    }
                    fsi->add_solid(move(solid));
                } catch (exception &e) {
                    throw runtime_error("Failed to parse entries with base name \'" + root_name_solids + "\'\n" +
                                        "Failed to parse solid geometry object number " + to_string(counter) + "\n" +
                                        string(e.what()));
                }
            }
        }
    } catch (exception &e) {
        throw runtime_error("Failed to parse input file:\n" + string(e.what()));
    }
}

string InputParser::option_not_specified_msg(string root_name, string option_name, string extra_msg) const {
    string msg = "Option \"" + option_name + "\" not specified in the input file \"" + input_filename + "\"";
    msg += " with base name \"" + root_name + "\"";
    if (!extra_msg.empty())
        msg += "\n" + extra_msg;
    return msg;
}
