//
// Created by anders on 2/24/22.
//
#include "simulate.hpp"
#include "input_parser.hpp"

void simulate(const std::string &input_file) {
    unique_ptr<fluid::FVM_Solver> fvm;
    unique_ptr<FSI_Solver> fsi;

    {
        InputParser input_parser{input_file};
        input_parser.create_solvers(fvm, fsi);
#pragma omp parallel
        { printf("Hello from thread %i\n", omp_get_thread_num()); }

        printf("---------------------- Setup complete ----------------------\n\n");
    }
    fsi->solve();
}
