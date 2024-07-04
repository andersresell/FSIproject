
#include "fluid/fvm_utilities.hpp"
#include "simulate.hpp"

int main(int argc, char *argv[]) {
    try {
        cout << "\n\n////////////////////////////////////////////\n"
             << "////////////// FSI Solver //////////////\n"
             << "////////////////////////////////////////////\n\n";
        if (argc != 2) {
            throw runtime_error("Specify the yaml input file\n");
        }

        const string input_file = argv[1];

        simulate(input_file);
    }

    catch (exception &e) {
        cerr << "Exception caught:\n" << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}