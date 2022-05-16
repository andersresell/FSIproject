//
// Created by anders on 1/4/22.
//

#include "setup_cases.hpp"

namespace fluid {


    void set_initial_cond_quiecent(vec4* U, int ni, int nj, double rho, double p){
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                U[IX(i, j)] = fluid::FVM_Solver::primitive2conserved({rho,0,0,p});
            }
        }
        std::cout << "Initial condition quiecent set\n";
    }

    void set_initial_cond1(vec4 *U, int ni, int nj) {
        //Sets the pressure in the bottom left part high and the pressure in the rest low
        double rho_low = 1.2;
        double p_low = 1e5;
        double p_high = 10 * p_low;
        double rho_high = rho_low*p_high/p_low;

        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                if (i <= (ni + 3) / 2 && j <= (nj + 3) / 2) {
                    U[IX(i, j)] = fluid::FVM_Solver::primitive2conserved({rho_high,0,0,p_high});
                } else {
                    U[IX(i, j)] = fluid::FVM_Solver::primitive2conserved({rho_low,0,0,p_low});
                }
            }
        }
        std::cout << "Initial condition 1 set\n";
    }

    void set_initial_cond2(vec4* U, int ni, int nj){
        //Sets the pressure in the left part high and the pressure in the rest low
        double rho_low = 1.2;
        double p_low = 1e5;
        double p_high = 10 * p_low;
        double rho_high = rho_low*p_high/p_low;
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                if (i <= ni / 10 + 2) {
                    U[IX(i, j)] = fluid::FVM_Solver::primitive2conserved({rho_high,0,0,p_high});
                } else {
                    U[IX(i, j)] = fluid::FVM_Solver::primitive2conserved({rho_low,0,0,p_low});
                }
            }
        }
        std::cout << "Initial condition 2 set\n";
    }

    void set_initial_cond_shock_tube_experiment(vec4* U, int ni, int nj, double L_x, vec4 V_l, vec4 V_r){
        //The left part is the "driver" section, and the right part is the "driven" section
        double dx = L_x/ni;
        double driver_length = 0.27;
        for (int i{0}; i < ni + 4; i++) {
            double x = (i-1.5)*dx;
            for (int j{0}; j < nj + 4; j++) {
                if (x <= driver_length) {
                    U[IX(i, j)] = fluid::FVM_Solver::primitive2conserved(V_l);
                } else {
                    U[IX(i, j)] = fluid::FVM_Solver::primitive2conserved(V_r);
                }
            }
        }
        std::cout << "Shock tube experiment set\n";
    }

    void set_constant_horizontal_flow_cond(vec4* U, int ni, int nj, double M, double rho, double p){
        double u = M*sqrt(Gamma*p/rho);
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                U[IX(i, j)].u1 = rho;
                U[IX(i, j)].u2 = rho * u;
                U[IX(i, j)].u3 = 0;
                U[IX(i, j)].u4 = p / (fluid::Gamma - 1) + 0.5 * rho * u * u;
            }
        }
        std::cout << "Initial condition horizontal constant flow, Mach "+std::to_string(M)+" set\n";
    }

    void set_initial_cond_riemann(vec4* U, int ni, int nj, const vec4& V_l, const vec4& V_r){
        if (ni%2 != 0) {
            std::cerr << "ni must be even for the riemann test\n";
            exit(1);
        }
        vec4 U_l = FVM_Solver::primitive2conserved(V_l);
        vec4 U_r = FVM_Solver::primitive2conserved(V_r);
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                if (i < ni/2+2){
                    U[IX(i,j)] = U_l;
                }
                else{
                    U[IX(i,j)] = U_r;
                }
            }
        }
        std::cout << "Initial condition Riemann set\n";
    }

    void set_initial_cond_pressure_bubble(vec4* U, int ni, int nj, double L_x, double L_y, double x_c,
                                          double y_c, double radius){
        double rho_low = 1.2;
        double p_low = 1e5;
        double p_high = 30 * p_low;
        double T_high = 273.15+200;
        double rho_high = p_high/(287*T_high);
        double u{0};
        double v{0};
        double R = L_x/10;
        double x, y;
        double dx = L_x/ni;
        double dy = L_y/nj;
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                x = dx*(i-1.5);
                y = dy*(j-1.5);
                U[IX(i, j)].u2 = 0;
                U[IX(i, j)].u3 = 0;
                if (sqr(x-x_c) + sqr(y-y_c) <= sqr(R)) {
                    U[IX(i, j)].u1 = rho_high;
                    U[IX(i, j)].u4 = p_high / (fluid::Gamma - 1) + 0.5 * rho_high * (u * u + v * v);
                } else {
                    U[IX(i, j)].u1 = rho_low;
                    U[IX(i, j)].u4 = p_low / (fluid::Gamma - 1) + 0.5 * rho_high * (u * u + v * v);
                }
            }
        }
        std::cout << "Initial condition 2 set\n";
    }

    void set_initial_cond_constant_data(vec4* U, int ni, int nj, vec4 V){
        vec4 U_const = FVM_Solver::primitive2conserved(V);
        for (int i{0}; i < ni + 4; i++) {
            for (int j{0}; j < nj + 4; j++) {
                U[IX(i, j)] = U_const;
            }
        }
        std::cout << "Initial condition constant data set\n";
    }




}

namespace solid {
    std::vector<solid::Point>
    generate_wedge(double l, double half_angle_deg, double x_tip, double y_tip){
        std::vector<solid::Point> wedge;
        double h = l * tan(half_angle_deg * M_PI / 180);
        wedge.push_back(solid::Point{x_tip,y_tip});
        wedge.push_back(solid::Point{x_tip+l,y_tip+h});
        wedge.push_back(solid::Point{x_tip+l,y_tip-h});
        return wedge;
    }

    std::vector<solid::Point>
    generate_diamond_wedge(double l, double half_angle_deg, double x_center, double y_center) {
        std::vector<solid::Point> diamond_wedge;
        double h = l * tan(half_angle_deg * M_PI / 180);
        diamond_wedge.push_back(solid::Point{x_center - l, y_center});
        diamond_wedge.push_back(solid::Point{x_center, y_center - h});
        diamond_wedge.push_back(solid::Point{x_center + l, y_center});
        diamond_wedge.push_back(solid::Point{x_center, y_center + h});
        return diamond_wedge;
    }

    std::vector<solid::Point> generate_circle(double R, int n_nodes, double x_center, double y_center) {
        std::vector<Point> circle;
        for (int i{0}; i < n_nodes; i++) {
            double theta = 2 * M_PI * i / n_nodes;
            circle.push_back(solid::Point{x_center + R * cos(theta), y_center + R * sin(theta)});
        }
        return circle;
    }

    std::vector<solid::Point>
    generate_rectangle(double W, double H, double rotation_angle_deg, double x_center, double y_center) {
        std::vector<solid::Point> rectangle;
        double theta = rotation_angle_deg * M_PI / 180;
        double c = cos(theta);
        double s = sin(theta);
        rectangle.push_back({-W / 2, -H / 2});
        rectangle.push_back({W / 2, -H / 2});
        rectangle.push_back({W / 2, H / 2});
        rectangle.push_back({-W / 2, H / 2});
        for (auto &p: rectangle) { //rotating all points and tranlating by the center
            double tmp_x = c * p.x - s * p.y;
            double tmp_y = s * p.x + c * p.y;
            p = {tmp_x + x_center, tmp_y + y_center};
        }
        return rectangle;
    }


}
