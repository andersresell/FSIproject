//
// Created by anders on 12/25/21.
//
#include "fvm_utilities.hpp"

void ExternalBCs::all_walls(){
    west = BC::Wall;
    east = BC::Wall;
    south = BC::Wall;
    north = BC::Wall;
}