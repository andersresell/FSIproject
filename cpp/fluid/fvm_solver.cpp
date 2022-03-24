//
// Created by anders on 12/26/21.
//

#include "fvm_solver.hpp"
#include "../solid/solid_body.hpp"

namespace fluid {
    using namespace std;

    FVM_Solver::FVM_Solver(int ni, int nj, double L_x, double L_y, double CFL, OdeScheme ode_scheme,
                           FluxScheme flux_scheme, const ExternalBCs& external_bcs,std::string output_folder)
            : ni{ni}, nj{nj}, L_x{L_x}, L_y{L_y}, dx{L_x / ni}, dy{L_y / nj}, CFL{CFL}, ode_scheme{ode_scheme},
              flux_scheme{flux_scheme}, external_bcs{external_bcs}, output_folder{std::move(output_folder)}{
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
        cell_status = new CellStatus[(ni + 4) * (nj + 4)];
        is_static = new bool[(ni + 4) * (nj + 4)];
        initialize_solids();
    }

    void FVM_Solver::write_fvm_output(const std::string& output_folder, int n) {
        //Writing one output data file for every timestep
        conserved2primitive(U);
        std::ofstream ost{"python/output_folders/" + output_folder +  "/fvm_output_t" + std::to_string(n) + ".csv"};
        if (!ost) {
            std::cerr << "Error: couldn't open fvm csv output file\n";
            exit(1);
        }
        ost << "#rho,u,v,p\n";
        int ind;
        for (int i{2}; i < ni + 2; i++) {
            for (int j{2}; j < nj + 2; j++) {
                ind = IX(i,j);
                if (cell_status[ind] == CellStatus::Solid){
                    ost << "nan,nan,nan,1e5\n";
                }
                else{ //Writing the ghost points as well
                    ost << V[ind].u1 << "," << V[ind].u2 << "," << V[ind].u3 << ","
                        << V[ind].u4 << "\n";
                }
            }
        }
    }
    void FVM_Solver::write_fvm_header(const std::string& output_folder, int write_stride, int n_last, double t_end) const{
        //Writing header file containg information about the simulation. This info is used by the python plotter
        std::ofstream ost{"python/output_folders/" + output_folder + "/fvm_header.csv"};
        if (!ost) {
            std::cerr << "Error: couldn't open fvm csv header file\n";
            exit(1);
        }
        ost << "#ni,nj,L_x,L_y,write_stride,n_last,t_end\n";
        ost << ni << ',' << nj << ',' << L_x << ',' << L_y << ',' << write_stride << ',' << n_last  << ',' << t_end << std::endl;

    }

    void FVM_Solver::initialize_solids() {
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                if (i == 0 || i == 1 || i == ni + 2 || i == ni + 3 || j == 0 || j == 1 || j == nj + 2 || j == nj + 3) {
                    cell_status[IX(i,j)] = CellStatus::Ghost; //setting all values to fluid, including the ones at the exterior
                } else {
                    cell_status[IX(i, j)] = CellStatus::Fluid;
                }
                is_static[IX(i, j)] = false;
            }
        }
    /*
        for (auto &sb_ptr: solid_bodies) {
            if (sb_ptr->type == solid::SolidBodyType::Static) {
                sb_ptr->find_solid_cells();
                sb_ptr->flag_static();
                sb_ptr->find_ghost_cells();
                sb_ptr->find_intercepts();
            }
        }
        solids_initialized = true;
        */
    }

    void FVM_Solver::step_solids(vec4* U_in, double dt, bool update_solid_pos){
        for (auto& sb_ptr : solid_bodies){
            sb_ptr->step(U_in,dt,update_solid_pos);
        }
    }

    double FVM_Solver::calc_solid_timestep(){
        double v_norm_max{0};
        for (auto &sb_ptr: solid_bodies) {
            v_norm_max = std::max(v_norm_max, sb_ptr->max_boundary_speed());
        }
        return std::min(dx,dy)/v_norm_max; //Could perhaps use the CFL condition here additionally
    }

    double FVM_Solver::ode_step(double dt_old) {
        external_bcs.set_BCs(U);
        step_solids(U,dt_old,true);
        double dt = std::min(calc_timestep(),calc_solid_timestep());
        switch (ode_scheme) {
            case OdeScheme::ExplicitEuler: {
                eval_RHS(U);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        if (cell_status[IX(i,j)] == CellStatus::Fluid) U[IX(i, j)] = U[IX(i, j)] + dt * Res[IXR(i, j)];
                    }
                }
                break;
            }
            case OdeScheme::TVD_RK3: {
                eval_RHS(U);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        if (cell_status[IX(i,j)] == CellStatus::Fluid) U_tmp[IX(i, j)] = U[IX(i, j)] + dt * Res[IXR(i, j)];

                    }
                }

                external_bcs.set_BCs(U_tmp);
                step_solids(U_tmp,dt, false);
                eval_RHS(U_tmp);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        if (cell_status[IX(i,j)] == CellStatus::Fluid)
                            U_tmp[IX(i, j)] = 3.0 / 4 * U[IX(i, j)] + 1.0 / 4 * U_tmp[IX(i, j)] + dt / 4 * Res[IXR(i, j)];
                    }
                }

                external_bcs.set_BCs(U_tmp);
                step_solids(U_tmp,dt, false);
                eval_RHS(U_tmp);
                for (int i{2}; i < ni + 2; i++) {
                    for (int j{2}; j < nj + 2; j++) {
                        if (cell_status[IX(i,j)] == CellStatus::Fluid)
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
                if (cell_status[IX(i,j)] == CellStatus::Fluid) {
                    //assert(!Res[IXR(i, j)].isnan());
                    Res[IXR(i, j)] =
                            -1 / dx * (F_f[IXF(i, j)] - F_f[IXF(i - 1, j)]) -
                            1 / dy * (G_f[IXG(i, j)] - G_f[IXG(i, j - 1)]);
                }
            }
        }
    }

    void FVM_Solver::MUSCL_extrapolate(vec4 *U_in) {
        vec4 Delta{};

        conserved2primitive(U_in);

        for (int i{1}; i < ni + 3; i++) {
            for (int j{2}; j < nj + 2; j++) {
                if (cell_status[IX(i,j)] != CellStatus::Solid) {
                    Delta = minmod(V[IX(i, j)] - V[IX(i - 1, j)], V[IX(i + 1, j)] - V[IX(i, j)]);
                    U_left[IXH(i, j)] = primitive2conserved(V[IX(i, j)] - 0.5 * Delta);
                    U_right[IXH(i, j)] = primitive2conserved(V[IX(i, j)] + 0.5 * Delta);
              }
            }
        }
        for (int i{2}; i < ni + 2; i++) {
            for (int j{1}; j < nj + 3; j++) {
                if (cell_status[IX(i,j)] != CellStatus::Solid) {
                    Delta = minmod(V[IX(i, j)] - V[IX(i, j - 1)], V[IX(i, j + 1)] - V[IX(i, j)]);
                    U_down[IXV(i, j)] = primitive2conserved(V[IX(i, j)] - 0.5 * Delta);
                    U_up[IXV(i, j)] = primitive2conserved(V[IX(i, j)] + 0.5 * Delta);

                }
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
                if (cell_status[IX(i,j)] != CellStatus::Solid) {
                    ind = IXH(i, j);
                    ind_p = IXH(i + 1, j);
                    F_f[IXF(i, j)] = 0.5 * (calc_F(U_right[ind]) + calc_F(U_left[ind_p])
                                            - std::max(calc_sprad_x(U_right[ind]), calc_sprad_x(U_left[ind_p])) *
                                              (U_left[ind_p] - U_right[ind]));
                }
            }
        }
        for (int i{2}; i < ni + 2; i++) {
            for (int j{1}; j < nj + 2; j++) {
                if (cell_status[IX(i,j)] != CellStatus::Solid) {
                    ind = IXV(i, j);
                    ind_p = IXV(i, j + 1);
                    G_f[IXG(i, j)] = 0.5 * (calc_G(U_up[ind]) + calc_G(U_down[ind_p])
                                            - std::max(calc_sprad_y(U_up[ind]), calc_sprad_y(U_down[ind_p])) *
                                              (U_down[ind_p] - U_up[ind]));
                }
            }
        }
    }


    double FVM_Solver::calc_timestep() const {
        double maxval{0};
        int ind;
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                ind = IX(i,j);
                if (cell_status[ind] != CellStatus::Solid)
                    maxval = std::max(maxval, calc_sprad_x(U[ind]) / dx + calc_sprad_y(U[ind]) / dy);

            }
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
        delete[] cell_status;
        delete[] is_static;

        if (ode_scheme == OdeScheme::TVD_RK3) {
            delete[] U_tmp;
        }
    }


}