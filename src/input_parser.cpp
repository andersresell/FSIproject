#include "input_parser.hpp"
#include "fluid/fvm_solver.hpp"

InputParser::InputParser(const string &input_filename) : input_filename{input_filename} {
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
    base_dir = std::filesystem::current_path();
    if (base_dir[base_dir.size() - 1] != '/')
        base_dir += '/';
    cout << "Simulation directory: \n" << base_dir;
}

void InputParser::create_solvers(unique_ptr<fluid::FVM_Solver> &fvm, unique_ptr<FSI_Solver> &fsi,
                                 unique_ptr<fluid::ExternalBCs> &bcs) {
    const string root_name_setup = "setup";
    const string root_name_ext_bcs = "external_bcs";
    const string root_name_initial_cond = "initial_condition";
    const string root_name_solids = "solids";
    try {

        /*--------------------------------------------------------------------
        Setup
        --------------------------------------------------------------------*/
        const int ni = read_required_option<int>(root_name_setup, "nx");
        const int nj = read_required_option<int>(root_name_setup, "ny");
        const double L_x = read_required_option<double>(root_name_setup, "L_x");
        const double L_y = read_required_option<double>(root_name_setup, "L_y");
        const double CFL = read_required_option<double>(root_name_setup, "CFL");
        const fluid::FluxScheme flux_scheme = read_optional_enum_option<fluid::FluxScheme>(
            root_name_setup, "flux_scheme", flux_scheme_from_string, fluid::FluxScheme::HLLC);
        const fluid::OdeScheme ode_scheme = read_optional_enum_option<fluid::OdeScheme>(
            root_name_setup, "ode_scheme", ode_scheme_from_string, fluid::OdeScheme::TVD_RK3);
        const int fvm_write_stride = read_required_option<int>(root_name_setup, "fvm_write_stride");

        const StoppingCrit stopping_crit_type = read_required_enum_option<StoppingCrit>(
            root_name_setup, "stopping_criterion_type", stopping_crit_from_string);
        const double stopping_crit_value = read_required_option<double>(root_name_setup, "stopping_criterion_value");

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

        bcs = make_unique<fluid::ExternalBCs>(ni, nj, west_bc, east_bc, south_bc, north_bc, M_inf, p_inf, rho_inf, "");

        const string output_dir = base_dir + "/output";
        fvm = make_unique<fluid::FVM_Solver>(ni, nj, L_x, L_y, CFL, ode_scheme, flux_scheme, *bcs, output_dir);

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
                    double M, rho, I;
                    solid::Point CM;
                    vector<solid::Point> boundary;
                    solid::SolidBodyType type = solid::SolidBodyType::Static;
                    if (solid_node["is_static"] && solid_node["is_static"].as<bool>() == false) {
                        type = solid::SolidBodyType::Dynamic;
                        rho = solid_node["rho"].as<double>();
                    }

                    YAML::Node geometry_node = solid_node["geometry"];
                    if (geometry_node) {
                        string geom_case;
                        if (geometry_node["case"]) {
                            geom_case = geometry_node["case"].as<string>();

                            if (geom_case == "circle") {
                                auto radius = geometry_node["radius"].as<double>();
                                auto n_nodes = geometry_node["n_nodes"].as<double>();
                                auto x_center = geometry_node["x_center"].as<double>();
                                auto y_center = geometry_node["y_center"].as<double>();
                                boundary = solid::generate_circle(radius, n_nodes, x_center, y_center);
                            } else if (geom_case == "wedge") {
                                auto l = geometry_node["l"].as<double>();
                                auto half_angle_deg = geometry_node["half_angle_deg"].as<double>();
                                auto x_tip = geometry_node["x_tip"].as<double>();
                                auto y_tip = geometry_node["y_tip"].as<double>();
                                boundary = solid::generate_wedge(l, half_angle_deg, x_tip, y_tip);
                                if (type == solid::SolidBodyType::Dynamic) {
                                    throw runtime_error("Dynamic properties not yet set for wedge geometry!\n");
                                }
                            } else if (geom_case == "diamond_wedge") {
                                auto l = geometry_node["l"].as<double>();
                                auto half_angle_deg = geometry_node["half_angle_deg"].as<double>();
                                auto x_center = geometry_node["x_center"].as<double>();
                                auto y_center = geometry_node["y_center"].as<double>();
                                boundary = solid::generate_diamond_wedge(l, half_angle_deg, x_center, y_center);
                                if (type == solid::SolidBodyType::Dynamic) {
                                    double h = l * sin(half_angle_deg * M_PI / 180);
                                    M = 2 * h * l * rho;
                                    I = rho / 3 * h * l *
                                        (sqr(h) + sqr(l)); // Calculated by solving a double integral, might be wrong
                                    CM = {x_center, y_center};
                                }
                            } else if (geom_case == "rectangle") {
                                auto W = geometry_node["width"].as<double>();
                                auto H = geometry_node["height"].as<double>();
                                auto rotation_angle_deg = geometry_node["rotation_angle_deg"].as<double>();
                                auto x_center = geometry_node["x_center"].as<double>();
                                auto y_center = geometry_node["y_center"].as<double>();
                                boundary = solid::generate_rectangle(W, H, rotation_angle_deg, x_center, y_center);
                                if (type == solid::SolidBodyType::Dynamic) {
                                    M = W * H * rho;
                                    I = M / 12 * (sqr(W) + sqr(H));
                                    CM = {x_center, y_center};
                                }
                            } else if (geom_case == "random_polygon") {
                                YAML::Node points_node = geometry_node["points"];
                                if (points_node) {
                                    for (const auto &p_node : points_node) {
                                        solid::Point p;
                                        p.x = p_node[0].as<double>();
                                        p.y = p_node[1].as<double>();
                                        boundary.push_back(p);
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
                    throw runtime_error("Failed to parse solid geometry object number " + to_string(counter) + "\n" +
                                        string(e.what()));
                }
            }
        }
    } catch (exception &e) {
        throw runtime_error("Failed to parse input file:\n" + string(e.what()));
    }
}

string InputParser::option_not_specified_msg(string root_name, string option_name, string extra_msg) const {
    string msg = "\"" + option_name + "\" not specified in the input file \"" + input_filename + "\"\n";
    msg += "\"with base name \"" + root_name + "\"";
    if (!extra_msg.empty())
        msg += "\n" + extra_msg;
    return msg;
}

void InputParser::check_that_root_name_is_valid(string root_name) const {

    auto it = std::find(available_root_nodes.begin(), available_root_nodes.end(), root_name);
    if (it == available_root_nodes.end()) {
        throw runtime_error("Illegal root name \'" + root_name + "\' encountered");
    }
}

void InputParser::add_parsed_option(string root_name, string option_name) {
    parsed_options.emplace_back(root_name, option_name);
}

static bool yaml_file_option_exist(const vector<pair<string, string>> &parsed_options,
                                   pair<string, string> option_pair) {
    auto it = std::find(parsed_options.begin(), parsed_options.end(), option_pair);
    return it != parsed_options.end();
}

void InputParser::report_invalid_options() const {
    for (const pair<string, string> &parsed_option : parsed_options) {
        const string &root_name = parsed_option.first;
        const string &option_name = parsed_option.second;
        if (!root_node[root_name][option_name]) {
            cout << "Warning: Input option with name \'" + option_name + "\' and root name \'" + root_name +
                        "\' is not specified in the input file\n";
        }
    }
    for (const string &root_name : available_root_nodes) {
        for (YAML::const_iterator it = root_node[root_name].begin(); it != root_node[root_name].end(); it++) {
            const string option_name = it->first.as<string>();
            const pair<string, string> option_pair = {root_name, option_name};
            if (!yaml_file_option_exist(parsed_options, option_pair)) {
                cout << "Warning: Input option with name \'" + option_name + "\' and root name \'" + root_name +
                            "\' specified in the input file is invalid\n";
            }
        }
    }
}

InputParser::~InputParser() {
    report_invalid_options();
}