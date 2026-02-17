#include "bend_stiffener.h"
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

double get_IE(double x) {
    double EI_base = 10000; // TODO: change this to the correct value
    double EI_tip = 1000;

    return EI_base -
           (EI_base - EI_tip) *
               (x /
                6.0); // TODO: length hardcoded for now but need to enter here
}

// define the exact equations for the bending of the bend stiffener

std::vector<State> bend_stiffener::calculate_strain(double tension,
                                                    double angle) {

    // boundary conditions based on cantelever beam modenl
    double y0 = 0;
    double theta0 = 0;
    double target_ML = 0;

    int steps = 1000;

    // calculate the moment and shear force at the root
    std::pair<double, double> result1;
    result1 = Integrator::solve_tapered_bvp(m_length, steps, y0, theta0,
                                            target_ML, angle);

    double M0 = result1.first;
    double V0 = result1.second;

    // integrate over the length

    double h = m_length / (double)steps;
    double x = 0.0; // start at the root

    State y = {y0, theta0, M0, V0}; // conditions at the root
    std::vector<State> Results;

    for (int i = 0; i < steps; i++) {
        y = (RK4(x, y, h, Integrator::equations));
        Results.push_back(y);
        x += h;
    };

    return Results;
}
