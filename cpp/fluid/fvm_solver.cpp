//
// Created by anders on 12/26/21.
//

#include "fvm_solver.hpp"

namespace fluid {
    using namespace std;

    FVM_Solver::FVM_Solver(int ni, int nj, double L_x, double L_y, double CFL, OdeScheme ode_scheme,
                           FluxScheme flux_scheme, ExternalBCs external_bcs)
            : ni{ni}, nj{nj}, L_x{L_x}, L_y{L_y}, dx{L_x / ni}, dy{L_y / nj}, CFL{CFL}, ode_scheme{ode_scheme},
            flux_scheme{flux_scheme}, external_bcs{external_bcs}
            {
        U = new vec4[(ni + 4) * (nj + 4)];
        V = new vec4[(ni + 4) * (nj + 4)];
        U_left = new vec4[(ni + 2) * nj];
        U_right = new vec4[(ni + 2) * nj];
        U_down = new vec4[ni * (nj + 2)];
        U_up = new vec4[ni * (nj + 2)];
        F_f = new vec4[(ni + 1) * nj];
        G_f = new vec4[ni * (nj + 1)];
        Res = new vec4[ni * nj];
        if (ode_scheme == OdeScheme::TVD_RK3) {
            U_tmp = new vec4[(ni + 4) * (nj + 4)];
        } else {
            U_tmp = nullptr;
        }
    }

    void FVM_Solver::write_fvm_csv_out_file(const std::string& output_folder, int n) {
        //Writing one output data file for every timestep
        conserved2primitive(U);
        std::ofstream ost{"../python/" + output_folder +  "/fvm_out_t" + std::to_string(n) + ".csv"};
        if (!ost) {
            std::cerr << "Error: couldn't open fvm csv output file\n";
            exit(1);
        }
        ost << "#rho,u,v,p\n";
        for (int i{2}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                ost << V[IX(i, j)].u1 << "," << V[IX(i, j)].u2 << "," << V[IX(i, j)].u3 << ","
                    << V[IX(i, j)].u4 << "\n";
            }
        }
    }
    void FVM_Solver::write_fvm_csv_header_file(const std::string& output_folder, int n_last, int write_stride) const{
        //Writing header file containg information about the simulation. This info is used by the python plotter
        std::ofstream ost{"../python/" + output_folder + "/header.csv"};
        if (!ost) {
            std::cerr << "Error: couldn't open fvm csv header file\n";
            exit(1);
        }
        ost << "#ni,nj,L_x,L_y,n_last,write_stride\n";
        ost << ni << ',' << nj << ',' << L_x << ',' << L_y << ',' << n_last << ',' << write_stride << std::endl;

    }


    double FVM_Solver::ode_step() {
        external_bcs.set_BCs(U); //External bc's are only applied once per timestep regardless of ode scheme for now
        double dt = calc_timestep();
        switch (ode_scheme) {
            case OdeScheme::ExplicitEuler: {
                eval_RHS(U);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        U[IX(i, j)] = U[IX(i, j)] + dt * Res[IXR(i, j)];
                    }
                }
                break;
            }
            case OdeScheme::TVD_RK3: {
                eval_RHS(U);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        U_tmp[IX(i, j)] = U[IX(i, j)] + dt * Res[IXR(i, j)];
                    }
                }
                external_bcs.set_BCs(U_tmp);
                eval_RHS(U_tmp);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        U_tmp[IX(i, j)] = 3.0 / 4 * U[IX(i, j)] + 1.0 / 4 * U_tmp[IX(i, j)] + dt / 4 * Res[IXR(i, j)];
                    }
                }
                external_bcs.set_BCs(U_tmp);
                eval_RHS(U_tmp);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        U[IX(i, j)] = 1.0 / 3 * U[IX(i, j)] + 2.0 / 3 * U_tmp[IX(i, j)] + 2 * dt / 3 * Res[IXR(i, j)];
                    }
                }
                break;
            }
        }
        return dt;
    }

    void FVM_Solver::eval_RHS(vec4 *U_in) {
        MUSCL_extrapolate(U_in); //eval U_left, U_right etc

        switch (flux_scheme) {
            case FluxScheme::Rusanov: {
                rusanov();
                break;
            }
            case FluxScheme::HLLC: {
                HLLC();
                break;
            }
        }

        for (int i{2}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                Res[IXR(i, j)] =
                        -1 / dx * (F_f[IXF(i, j)] - F_f[IXF(i - 1, j)]) -
                        1 / dy * (G_f[IXG(i, j)] - G_f[IXG(i, j - 1)]);
            }
        }
    }

    void FVM_Solver::MUSCL_extrapolate(vec4 *U_in) {
        vec4 Delta;
        conserved2primitive(U_in);

        for (int i{1}; i < ni + 3; i++) {
            for (int j{2}; j < nj + 2; j++) {
                Delta = minmod(V[IX(i, j)] - V[IX(i - 1, j)], V[IX(i + 1, j)] - V[IX(i, j)]);
                U_left[IXH(i, j)] = primitive2conserved(V[IX(i, j)] - 0.5 * Delta);
                U_right[IXH(i, j)] = primitive2conserved(V[IX(i, j)] + 0.5 * Delta);
            }
        }
        for (int i{2}; i < ni + 2; i++) {
            for (int j{1}; j < nj + 3; j++) {
                Delta = minmod(V[IX(i, j)] - V[IX(i, j - 1)], V[IX(i, j + 1)] - V[IX(i, j)]);
                U_down[IXV(i, j)] = primitive2conserved(V[IX(i, j)] - 0.5 * Delta);
                U_up[IXV(i, j)] = primitive2conserved(V[IX(i, j)] + 0.5 * Delta);
            }
        }
    }

    vec4 FVM_Solver::minmod(const vec4 &a, const vec4 &b) {
        return {sgn(a.u1) * std::max(0.0, std::min(std::abs(a.u1), sgn(a.u1) * b.u1)),
                sgn(a.u2) * std::max(0.0, std::min(std::abs(a.u2), sgn(a.u2) * b.u2)),
                sgn(a.u3) * std::max(0.0, std::min(std::abs(a.u3), sgn(a.u3) * b.u3)),
                sgn(a.u4) * std::max(0.0, std::min(std::abs(a.u4), sgn(a.u4) * b.u4))};
    }

    void FVM_Solver::conserved2primitive(vec4 *U_in) {
        for (int i{0}; i < (ni + 4) * (nj + 4); i++) {
            V[i].u1 = U_in[i].u1;
            V[i].u2 = U_in[i].u2 / U_in[i].u1;
            V[i].u3 = U_in[i].u3 / U_in[i].u1;
            V[i].u4 = calc_P(U_in[i]);
        }
    }

    void FVM_Solver::rusanov() {
        unsigned int ind, ind_p;
        for (int i{1}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                ind = IXH(i, j);
                ind_p = IXH(i + 1, j);
                F_f[IXF(i, j)] = 0.5 * (calc_F(U_right[ind]) + calc_F(U_left[ind_p])
                                        - std::max(calc_sprad_x(U_right[ind]), calc_sprad_x(U_left[ind_p])) *
                                          (U_left[ind_p] - U_right[ind]));
            }
        }
        for (int i{2}; i < ni + 2; i++) {
            for (int j{1}; j < nj + 2; j++) {
                ind = IXV(i, j);
                ind_p = IXV(i, j + 1);
                G_f[IXG(i, j)] = 0.5 * (calc_G(U_up[ind]) + calc_G(U_down[ind_p])
                                        - std::max(calc_sprad_y(U_up[ind]), calc_sprad_y(U_down[ind_p])) *
                                          (U_down[ind_p] - U_up[ind]));
            }
        }
    }


    double FVM_Solver::calc_timestep() const {
        double maxval{0};
        for (int i{0}; i < (ni + 4) * (nj + 4); i++) {
            maxval = std::max(maxval, calc_sprad_x(U[i]) / dx + calc_sprad_y(U[i]) / dy);
        }
        return CFL / maxval;
    }

    FVM_Solver::~FVM_Solver() {
        delete[] U;
        delete[] V;
        delete[] U_left;
        delete[] U_right;
        delete[] U_down;
        delete[] U_up;
        delete[] F_f;
        delete[] G_f;
        delete[] Res;

        if (ode_scheme == OdeScheme::TVD_RK3) {
            delete[] U_tmp;
        }
    }
}