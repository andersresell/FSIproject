//
// Created by anders on 12/25/21.
//
#include "fvm_utilities.hpp"

namespace fluid {

    ExternalBCs::ExternalBCs(int ni, int nj, BC_Type west, BC_Type east, BC_Type south, BC_Type north,
                             std::string history_output_west_folder,
                             double M_inf, double p_inf, double rho_inf)
            : ni{ni}, nj{nj}, west{west}, east{east}, south{south}, north{north}{
        if (east == BC_Type::SupersonicInflow || north == BC_Type::SupersonicInflow || south == BC_Type::SupersonicInflow) {
            std::cerr << "Supersonic inflow only permitted at western wall\n";
            exit(1);
        }
        double u_inf = M_inf*sqrt(Gamma*p_inf/rho_inf);
        U_inf = {rho_inf,rho_inf*u_inf, 0, p_inf/(Gamma-1) + 0.5*rho_inf*u_inf*u_inf};
        if (west == BC_Type::TimeHistory){
            load_history_output_west(std::move(history_output_west_folder));
        }
    }

    void ExternalBCs::set_BCs(vec4 *U_in, double t) {
        switch (west) {
            case BC_Type::InvicidWall : {
                for (int j{0}; j < nj + 4; j++) {
                    U_in[IX(0, j)] = set_vertical_invicid_wall(U_in[IX(3, j)]);
                    U_in[IX(1, j)] = set_vertical_invicid_wall(U_in[IX(2, j)]);
                }
                break;
            }
            case BC_Type::SupersonicInflow : {
                for (int j{0}; j < nj + 4; j++) {
                    U_in[IX(0, j)] = U_inf;
                    U_in[IX(1, j)] = U_inf;
                }
                break;
            }
            case BC_Type::NonreflectingOutflow : {
                for (int j{ 0 }; j < nj + 4; j++) {
                    U_in[IX(0, j)] = U_in[IX(3, j)];
                    U_in[IX(1, j)] = U_in[IX(2, j)];
                }
                break;
            }
            case BC_Type::TimeHistory : {
                assert(!time_history_west.empty());
                auto it = std::find_if(time_history_west.begin(), time_history_west.end(),
                                       [t](TimeHistory x) { return x.t >= t; });
                vec4 U0, U1;
                if (it == time_history_west.begin()){
                    U0 = it->U_outer_GP;
                    U1 = it->U_inner_GP;
                }else{ //Linear interpolation of time values
                    double epsilon = (t-(it-1)->t)/(it->t - (it-1)->t);
                    U0 = epsilon*(it->U_outer_GP) + (1 - epsilon)*((it-1)->U_outer_GP);
                    U1 = epsilon*(it->U_inner_GP) + (1 - epsilon)*((it-1)->U_inner_GP);
                }
                for (int j{ 0 }; j < nj + 4; j++) {
                    U_in[IX(0, j)] = U0;
                    U_in[IX(1, j)] = U1;
                }

                break;
            }
        }
        switch (east) {
            case BC_Type::InvicidWall : {
                for (int j{ 0 }; j < nj + 4; j++) {
                   U_in[IX(ni + 3, j)] = set_vertical_invicid_wall(U_in[IX(ni, j)]);
                   U_in[IX(ni + 2, j)] = set_vertical_invicid_wall(U_in[IX(ni + 1, j)]);
                }
                break;
            }
            case BC_Type::NonreflectingOutflow : {
                for (int j{ 0 }; j < nj + 4; j++) {
                    U_in[IX(ni + 3, j)] = U_in[IX(ni + 1, j)];
                    U_in[IX(ni + 2, j)] = U_in[IX(ni + 1, j)];
                }
                break;
            }

        }
        switch (south) {
            case BC_Type::InvicidWall : {
                for (int i{0}; i < ni + 4; i++) {
                    U_in[IX(i, 0)] = set_horizontal_invicid_wall(U_in[IX(i, 3)]);
                    U_in[IX(i, 1)] = set_horizontal_invicid_wall(U_in[IX(i, 2)]);
                }
                break;
            }
            case BC_Type::NonreflectingOutflow : {
                for (int i{ 0 }; i < ni + 4; i++) {
                    U_in[IX(i, 0)] = U_in[IX(i, 2)];
                    U_in[IX(i, 1)] = U_in[IX(i, 2)];
                }
                break;
            }
        }
        switch (north) {
            case BC_Type::InvicidWall : {
                for (int i{0}; i < ni + 4; i++) {
                    U_in[IX(i, nj + 3)] = set_horizontal_invicid_wall(U_in[IX(i, nj)]);
                    U_in[IX(i, nj + 2)] = set_horizontal_invicid_wall(U_in[IX(i, nj + 1)]);
                }
                break;
            }
            case BC_Type::NonreflectingOutflow : {
                for (int i{ 0 }; i < ni + 4; i++) {
                    U_in[IX(i, nj + 3)] = U_in[IX(i, nj + 1)];
                    U_in[IX(i, nj + 2)] = U_in[IX(i, nj + 1)];
                }
                break;
            }
        }
    }


    vec4 ExternalBCs::set_vertical_invicid_wall(const vec4& U_in){
        //Enforced by setting x velocity component at ghost cell to the nagative value of internal cell. This can be done
        //by switching the sign of the x momentum.
        return {U_in.u1,
        -U_in.u2,
        U_in.u3,
        U_in.u4};
    }

    vec4 ExternalBCs::set_horizontal_invicid_wall(const vec4& U_in){
        //Enforced by setting y velocity component at ghost cell to the nagative value of internal cell. This can be done
        //by switching the sign of the y momentum.
        return {U_in.u1,
                U_in.u2,
                -U_in.u3,
                U_in.u4};
    }

    void ExternalBCs::load_history_output_west(std::string&& history_output_west_folder){
        std::ifstream ist{"python/output_folders/" + history_output_west_folder +  "/fvm_history_output_west.csv"};
        if (!ist) {
            std::cerr << "Error: couldn't open fvm history output west file\n";
            exit(1);
        }
        std::string line;
        double t;
        vec4 U0;
        vec4 U1;
        //while (std::getline(ist,line)){
          //  std::stringstream ss{line};
            //while (getline(ss))
        //}
        std::string first_line;
        getline(ist,first_line); //skipping first line
        while (ist >> t >> U0.u1 >> U0.u2 >> U0.u3 >> U0.u4 >> U1.u1 >> U1.u2 >> U1.u3 >> U1.u4){
            time_history_west.push_back({t,U1,U0});
        }

        std::cout << "Fvm time history at western boundary successfully loaded\n";

    }

}
