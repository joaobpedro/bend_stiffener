#include "bend_stiffener.h"
#include "bs_physics.h"
#include "integration_helper.h"
#include <cmath>
#include <vector>

typedef std::vector<double> State;

bend_stiffener::bend_stiffener(double root_dia, double length, double tip_dia,
                               double inner_dia)
    : m_root_dia(root_dia), m_length(length), m_tip_dia(tip_dia),
      m_inner_dia(inner_dia) {};

bend_stiffener::~bend_stiffener() {};

double bend_stiffener::calculate_volume() {
    double volume_outer;
    double volume_inner;
    double volume;

    const double pi = 3.141592653599;

    volume_outer = 1.0 / 3.0 * pi * m_length *
                   (pow(m_root_dia / 2, 2) + pow(m_tip_dia / 2, 2) +
                    m_root_dia / 2 + m_tip_dia / 2);
    volume_inner = pi * pow(m_inner_dia / 2, 2) * m_length;
    volume = volume_outer - volume_inner;

    return volume;
};

double bend_stiffener::get_dia(double x) {
    return m_root_dia - (m_root_dia - m_tip_dia) * (x / m_length);
};

double bend_stiffener::get_Inertia(double x) {
    const double pi = 3.141592653599;
    // calculate outer diameter
    // this assumes the conic shape
    double outer_dia = m_root_dia - (m_root_dia - m_tip_dia) * (x / m_length);

    // inertia is given by  I = pi/64 *(Do^4 - Di^4)

    double inertia = (pi / 64) * (pow(outer_dia, 4) - pow(m_inner_dia, 4));
    return inertia;
}

std::vector<State> bend_stiffener::solve_equations(double tension,
                                                   double angle) {

    // boundary conditions based on cantelever beam model
    double y0 = 0;
    double theta0 = 0;
    double target_ML = 0;

    int steps =
        1000; // this always discretized the bend stiffener into 1000 pieces

    // calculate the moment and shear force at the root
    std::pair<double, double> result1;
    result1 = Integrator::solve_tapered_bvp(
        m_length, steps, y0, theta0, target_ML, angle, bs_physics::equations,
        bend_stiffener::get_Inertia(0.0));

    double M0 = result1.first;
    double V0 = result1.second;

    // integrate over the length

    double h = m_length / (double)steps;
    double x = 0.0; // start at the root

    State y = {y0, theta0, M0, V0}; // conditions at the root
    std::vector<State> Results;

    double inertia;
    for (int i = 0; i < steps; i++) {
        inertia = bend_stiffener::get_Inertia(x);
        y = (Integrator::RK4(x, y, h, bs_physics::equations, inertia));
        Results.push_back(y);
        x += h;
    };

    return Results;
}

State bend_stiffener::calculate_strain(std::vector<State> Deformations) {

    // curvature is d_theta/d_s
    double step = m_length / 1000; // HARDCODED: discretized
    double x = 0.0;
    //
    State strain;

    for (int i = 0; i < Deformations.size(); i++) {
        double EI = bs_physics::get_EI(bend_stiffener::get_Inertia(x));
        strain.push_back((bend_stiffener::get_dia(x) * Deformations[i][2]) /
                         (EI));
        x += step;
    }

    return strain;
}
