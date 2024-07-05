#pragma once
#include "fvm_utilities.hpp"

namespace fluid {

struct BlastLoadBC {
    /*Refer to "A coupling of empirical explosive blast loads to ALE air domains in LS-DYNA"

    the paramter t_bo which is a shift in time is assumed to be zero

    */
    double p_s;           // peak overpressure
    double D_p;           // overpressure duration phase
    double t_a;           // blast arrival time
    double zeta;          //?
    double phi;           //?
    Vector2 blast_center; // origin of the blast
    double alpha;
    double rho0;
    double p0;
    double c0;

    double calc_tau(double t) const { return (t - t_a) / D_p; }
    double calc_p_inc(double t) const {
        double tau = calc_tau(t);
        return p_s * (1 - tau) * exp(-alpha * tau);
    }
    double calc_rho1() const; //{}
    double calc_p1() const;   //{}
    double calc_p0() const;   //{}

    vec4 calc_bc_val_conservative(Vector2 coord, double t);
};

inline vec4 BlastLoadBC::calc_bc_val_conservative(Vector2 coord, double t) {

    assert(false); // need to account for what happens before the blast arrives.

    const Vector2 r = coord - blast_center;
    const double R = r.norm(); // needed?
    const double p1 = calc_p1();
    const double p_inc = calc_p_inc(t);
    const double rho1 = rho0 * (6 * p1 / p0 + 1) / (p1 / p0 + 6);
    const double rho = rho1 * pow((p_inc + p0) / p1, 1.0 / Gamma);
    const double u_p_sqr = sqr(c0) * 25 * sqr(p1 / p0 - 1) / (42 * p1 / p0 + 7);
    const double q_s = 0.5 * rho1 * u_p_sqr;
    const double q = q_s * (1 - zeta) * exp(-phi * zeta);
    const double u_mag = sqrt(2 * q / rho);
    const Vector2 u = u_mag * r.normalized();

    const vec4 V = {rho, u.x(), u.y(), p_inc};
    return primitive2conserved(V);
}

struct ExternalBCs {
    BlastLoadBC blast_bc; // used in case the type is blast load
    int ni, nj;
    double dx, dy;
    BC_Type west, east, south, north;
    vec4 U_inf;                            // Used in case of supersonic inflow
    vector<TimeHistory> time_history_west; // Used for time history at western boundary

    ExternalBCs(int ni, int nj, double dx, double dy, BC_Type west, BC_Type east, BC_Type south, BC_Type north,
                double M_inf, double p_inf, double rho_inf, string history_output_west_folder);

    // Used in case BlastLoad is used
    ExternalBCs() = default;

    Vector2 get_coord(int i, int j) { return {dx * (i - 1.5), dy * (j - 1.5)}; }
    void set_BCs(vec4 *U_in, double t);
    static vec4 set_vertical_invicid_wall(const vec4 &U_in);
    static vec4 set_horizontal_invicid_wall(const vec4 &U_in);

    void load_history_output_west(string history_output_west_folder);
};

inline ExternalBCs create_external_bcs_blast_load(int ni, int nj, double Lx, double Ly, Vector2 blast_center,
                                                  double rho0, double p0) {

    BlastLoadBC bl;
    bl.rho0 = rho0;
    bl.p0 = p0;
    bl.c0 = sqrt(Gamma * p0 / rho0);
    bl.blast_center = blast_center;

    ExternalBCs bcs;
    bcs.ni = ni;
    bcs.nj = nj;

    if (blast_center.x() >= 0 && blast_center.x() <= Lx && blast_center.x() >= 0 && blast_center.y() <= Ly) {
        string err_msg = "Blast center location (x,y) = (" + to_string(blast_center.x()) + ", " +
                         to_string(blast_center.y()) + ") is contained inside the simulation domain";
        throw runtime_error(err_msg);
    }

    /*Each face that is "facing" the blast is assigned blast condition, otherwise nonreflecting outflow is assigned*/
    bcs.west = (blast_center.x() < 0.0) ? BC_Type::BlastLoad : BC_Type::NonreflectingOutflow;
    bcs.east = (blast_center.x() > Lx) ? BC_Type::BlastLoad : BC_Type::NonreflectingOutflow;
    bcs.south = (blast_center.y() < 0.0) ? BC_Type::BlastLoad : BC_Type::NonreflectingOutflow;
    bcs.north = (blast_center.y() > Ly) ? BC_Type::BlastLoad : BC_Type::NonreflectingOutflow;

    return bcs;
}

} // namespace fluid
