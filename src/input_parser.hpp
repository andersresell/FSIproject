#pragma once

#include "fsi/fsi_solver.hpp"
#include "includes.hpp"
#include <yaml-cpp/yaml.h>

class InputParser {

  public:
    InputParser(const string &input_filename);

    void create_solvers(unique_ptr<fluid::FVM_Solver> &fvm, unique_ptr<FSI_Solver> &fsi);

  private:
    string option_not_specified_msg(string root_name, string option_name, string extra_msg = "") const;

    template <typename T> static T read_option(const YAML::Node &node, string option_name) {
        if (node[option_name]) {
            return node[option_name].as<T>();
        } else {
            throw runtime_error("Failed reading option \'" + option_name + "\'");
        }
    }

    template <typename T> T read_required_option(string root_name, string option_name, string extra_msg = "") {

        if (root_node[root_name][option_name]) {
            return root_node[root_name][option_name].as<T>();
        } else {
            throw std::runtime_error(option_not_specified_msg(root_name, option_name, extra_msg));
        }
    }

    template <typename T> T read_optional_option(string root_name, string option_name, T default_value) {
        if (root_node[root_name][option_name]) {
            return root_node[root_name][option_name].as<T>();
        } else {
            return default_value;
        }
    }

    template <typename EnumType>
    static EnumType lookup_enum_option_map(const std::map<string, EnumType> &map, const string &key, string option_name,
                                           const string &option_parent_name = "") {
        if (map.count(key) == 1) {
            return map.at(key);
        } else {
            string keys;
            for (const auto &pair : map)
                keys += "'" + pair.first + "'\n";
            if (option_parent_name.size() > 0)
                option_name = option_parent_name + ": " + option_name;
            throw std::runtime_error("Illegal value '" + key + "' specified for setting '" + option_name +
                                     "'. Legal values are:\n" + keys);
        }
    }

    template <typename EnumType>
    EnumType read_required_enum_option(string root_name, string option_name,
                                       const std::map<string, EnumType> &enum_map) {
        if (root_node[root_name][option_name]) {
            return lookup_enum_option_map(enum_map, root_node[root_name][option_name].as<string>(), option_name);
        } else {
            throw std::runtime_error(option_not_specified_msg(root_name, option_name));
        }
    }

    template <typename EnumType>
    EnumType read_optional_enum_option(string root_name, string option_name, const std::map<string, EnumType> &enum_map,
                                       EnumType default_value) {
        if (root_node[root_name][option_name]) {
            return lookup_enum_option_map(enum_map, root_node[root_name][option_name].as<string>(), option_name);
        } else {
            return default_value;
        }
    }
    string base_dir;
    string output_dir;
    YAML::Node root_node;
    string input_filename;

    static inline const map<string, fluid::BC_Type> bc_type_from_string{
        {"inviscid_wall", fluid::BC_Type::InvicidWall},
        {"nonreflecting_outflow", fluid::BC_Type::NonreflectingOutflow},
        {"supersonic_inflow", fluid::BC_Type::SupersonicInflow},
        {"time_history", fluid::BC_Type::TimeHistory}};

    static inline const map<string, fluid::FluxScheme> flux_scheme_from_string{{"hllc", fluid::FluxScheme::HLLC},
                                                                               {"rusanov", fluid::FluxScheme::Rusanov}};

    static inline const map<string, fluid::OdeScheme> ode_scheme_from_string{
        {"explicit_euler", fluid::OdeScheme::ExplicitEuler}, {"tvd_rk3", fluid::OdeScheme::TVD_RK3}};

    static inline const map<string, fluid::Limiter> limiter_from_string{{"minmod", fluid::Limiter::Minmod},
                                                                        {"mc", fluid::Limiter::MC}};

    static inline const map<string, StoppingCrit> stopping_crit_from_string{{"time", StoppingCrit::Time},
                                                                            {"timesteps", StoppingCrit::Timesteps},
                                                                            {"convergence", StoppingCrit::Convergence}};
};