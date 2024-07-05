#pragma once
//
// Created by anders on 12/26/21.
//
// A carthesian FVM solver employing a semi discrete approach, using an ODE solver. It uses MUSCL extrapolation
// with minmod limiter to get second order accuracy and an approximate Riemann solver to compute the fluxes.
// The boundaries are handled by ghost points

#include "../includes.hpp"
#include "fvm_bcs.hpp"
#include "fvm_utilities.hpp"

namespace solid {
// Forward declaration
class SolidBody;
} // namespace solid

namespace fluid {
// typedef vec4 (*bc_func_ptr)(const vec4& U_in); //In order to set appropriate BC's functions are passed as arguents
// and return types

enum class CellStatus {
    Fluid = 0,
    Ghost,
    Solid
};

class FVM_Solver {
    // Computes one timestep of the fluid domain with the explicit euler method
  public:
    const int ni, nj;
    const double L_x, L_y;
    const double dx, dy;
    const double CFL;
    const OdeScheme ode_scheme;
    const FluxScheme flux_scheme;
    const Limiter limiter;
    ExternalBCs external_bcs;
    std::vector<std::shared_ptr<solid::SolidBody>> solid_bodies;
    CellStatus *cell_status;
    bool *is_static;
    const std::string output_folder;
    // bool first_timestep;
    // bool solids_initialized;
    vec4 *U;

  private:
    vec4 *U_tmp, *V, *Res, *U_left, *U_right, *U_down, *U_up, *F_f, *G_f;

  public:
    FVM_Solver(int ni, int nj, double L_x, double L_y, double CFL, OdeScheme ode_scheme, FluxScheme flux_scheme,
               Limiter limiter, const ExternalBCs &external_bcs, std::string output_folder);

    void write_fvm_output(const std::string &output_folder, int n, double t);
    void write_fvm_header(const std::string &output_folder, int write_stride, int n_last, double t_end) const;

    double ode_step(double dt_old, double t_old);

    void initialize_solids();

  private:
    void step_solids(vec4 *U_in, double dt, bool update_solid_pos); // Enforces fluid bc's and updates solid positions

    double calc_solid_timestep();

    double calc_timestep() const;

    void eval_RHS(vec4 *U_in);

    void MUSCL_extrapolate(vec4 *U_in);

    void rusanov();

    void HLLC();

    static vec4 F_f_HLLC(const vec4 &U_L, const vec4 &U_R); // Computes the HLLC flux F_{i+1/2,j}

    static vec4 G_f_HLLC(const vec4 &U_D, const vec4 &U_U); // Computes the HLLC flux G_{i,j+1/2}

  public:
    ~FVM_Solver();
};

} // namespace fluid
