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
        U_left = new vec4[(ni + 1) * nj];
        U_right = new vec4[(ni + 1) * nj];
        U_down = new vec4[ni * (nj + 1)];
        U_up = new vec4[ni * (nj + 1)];
        F_f = new vec4[(ni + 1) * nj];
        G_f = new vec4[ni * (nj + 1)];
        Res = new vec4[ni * nj];
        if (ode_scheme == OdeScheme::TVD_RK3) {
            U_tmp = new vec4[(ni + 4) * (nj + 4)];
        } else {
            U_tmp = nullptr;
        }
    }

    void FVM_Solver::write_simple_fvm_csv_file(std::string filename) {
        std::ofstream ost{"../python/fvm_output/" + filename};
        if (!ost) {
            std::cerr << "Error: couldn't open file\n";
            exit(1);
        }
        //Header
        ost << "#ni,nj,L_x,L_y\n";
        ost << ni << "," << nj << "," << L_x << "," << L_y << std::endl;
        ost << "#u1,u2,u3,u4\n";
        for (int i{2}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                ost << U[IX(i, j)].u1 << "," << U[IX(i, j)].u2 << "," << U[IX(i, j)].u3 << ","
                    << U[IX(i, j)].u4 << "\n";
            }
        }
    }


    double FVM_Solver::ode_step() {
        external_bcs.set_BCs(U); //External bc's are only applied once per timestep regardless of ode scheme for now
        double dt = calc_timestep(CFL);
        std::cout << "dt = " << dt << std::endl;
        switch (ode_scheme) {
            case OdeScheme::ExplicitEuler: {
                eval_RHS(U);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        U[IX(i, j)] = U[IX(i, j)] + dt * Res[IXR(i - 2, j - 2)];
                    }
                }
                break;
            }
            case OdeScheme::TVD_RK3: {
                //IMPLEMENT
                break;
            }
        }
        return dt;
    }

    void FVM_Solver::eval_RHS(vec4 *U_in) {
        //eval U_left, U_right etc
        MUSCL_extrapolate(U_in);
        switch (flux_scheme) {
            case FluxScheme::Rusanov: {
                rusanov();
                break;
            }
        }

        for (int i{0}; i < ni; i++) {
            for (int j{0}; j < nj; j++) {
                Res[IXR(i, j)] =
                        -1 / dx * (F_f[IXH(i + 1, j)] - F_f[IXH(i, j)]) -
                        1 / dy * (G_f[IXV(i, j + 1)] - G_f[IXV(i, j)]);
            }
        }
    }

    void FVM_Solver::MUSCL_extrapolate(vec4 *U_in) {
        vec4 tmpV;
        conserved2primitive(U_in);

        for (int i{1}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                tmpV = V[IX(i, j)] + 0.5 * minmod(V[IX(i, j)] - V[IX(i - 1, j)], V[IX(i + 1, j)] - V[IX(i, j)]);
                U_left[IXH(i - 1, j - 2)] = primitive2conserved(tmpV);

                tmpV = V[IX(i + 1, j)] - 0.5 * minmod(V[IX(i + 2, j)] - V[IX(i + 1, j)], V[IX(i + 1, j)] - V[IX(i, j)]);
                U_right[IXH(i - 1, j - 2)] = primitive2conserved(tmpV);
            }
        }
        for (int i{2}; i < ni + 2; i++) {
            for (int j{1}; j < nj + 2; j++) {
                tmpV = V[IX(i, j)] + 0.5 * minmod(V[IX(i, j)] - V[IX(i, j - 1)], V[IX(i, j + 1)] - V[IX(i, j)]);
                U_down[IXV(i - 2, j - 1)] = primitive2conserved(tmpV);

                tmpV = V[IX(i, j + 1)] - 0.5 * minmod(V[IX(i, j + 2)] - V[IX(i, j + 1)], V[IX(i, j + 1)] - V[IX(i, j)]);
                U_up[IXV(i - 2, j - 1)] = primitive2conserved(tmpV);
            }
        }
    }

    inline vec4 FVM_Solver::minmod(const vec4 &a, const vec4 &b) const {
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

    vec4 FVM_Solver::conserved2primitive(const vec4 &U_in) {
        return {U_in.u1,
                U_in.u2 / U_in.u1,
                U_in.u3 / U_in.u1,
                calc_P(U_in)};
    }

    double FVM_Solver::calc_P(const vec4 &U_in) {
        return (Gamma - 1) * (U_in.u4 - 0.5 * (U_in.u2 * U_in.u2 + U_in.u3 * U_in.u3) / U_in.u1);
    }

    vec4 FVM_Solver::primitive2conserved(const vec4 &V_in) {
        return {V_in.u1,
                V_in.u2 * V_in.u1,
                V_in.u3 * V_in.u1,
                V_in.u4 / (Gamma - 1) + 0.5 * V_in.u1 * (V_in.u2 * V_in.u2 + V_in.u3 * V_in.u3)};
    }

    inline vec4 FVM_Solver::calc_F(const vec4 &U_in) const {
        double P = calc_P(U_in);
        return {U_in.u2,
                U_in.u2 * U_in.u2 / U_in.u1 + P,
                U_in.u2 * U_in.u3 / U_in.u1,
                (U_in.u4 + P) * U_in.u2 / U_in.u1};
    }

    inline vec4 FVM_Solver::calc_G(const vec4 &U_in) const {
        double P = calc_P(U_in);
        return {U_in.u3,
                U_in.u2 * U_in.u3 / U_in.u1,
                U_in.u3 * U_in.u3 / U_in.u1 + P,
                (U_in.u4 + P) * U_in.u3 / U_in.u1};
    }

    void FVM_Solver::rusanov() {
        for (int i{0}; i < ni + 1; i++) {
            for (int j{0}; j < nj; j++) {
                int ind = IXH(i, j);
                F_f[ind] = 0.5 * (calc_F(U_left[ind]) + calc_F(U_right[ind])
                                  - std::max(calc_sprad_x(U_left[ind]), calc_sprad_x(U_right[ind])) *
                                    (U_right[ind] - U_left[ind]));
            }
        }
        for (int i{0}; i < ni; i++) {
            for (int j{0}; j < nj + 1; j++) {
                int ind = IXV(i, j);
                G_f[ind] = 0.5 * (calc_G(U_down[ind]) + calc_G(U_up[ind])
                                  - std::max(calc_sprad_y(U_down[ind]), calc_sprad_y(U_up[ind])) *
                                    (U_up[ind] - U_down[ind]));
            }
        }
    }

    double FVM_Solver::calc_sprad_x(const vec4 &U_in) const {
        return std::abs(U_in.u2 / U_in.u1) + calc_sound_speed(U_in);
    }

    double FVM_Solver::calc_sprad_y(const vec4 &U_in) const {
        return std::abs(U_in.u3 / U_in.u1) + calc_sound_speed(U_in);
    }


    inline double FVM_Solver::calc_sound_speed(const vec4 &U_in) const {
        return sqrt(Gamma / U_in.u1 * calc_P(U_in));
    }


    double FVM_Solver::calc_timestep(double CFL) const {
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