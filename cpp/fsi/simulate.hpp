//
// Created by anders on 2/23/22.
//

#ifndef FSIPROJECT_SIMULATE_HPP
#define FSIPROJECT_SIMULATE_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "../../includes.hpp"
#include "fsi_solver.hpp"

struct Simulate{
    Simulate(const std::string& input_file);

    fluid::BC_Type bc_type_from_str(std::string&& bc_type);

};



#endif //FSIPROJECT_SIMULATE_HPP
