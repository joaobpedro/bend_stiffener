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

double bend_stiffener::get_EI(double x) {
    double EI_base = 10000; // TODO: change this to the correct value
    double EI_tip = 1000;

    return EI_base -
           (EI_base - EI_tip) *
               (x /
                6.0); // TODO: length hardcoded for now but need to enter here
}

// define the exact equations for the bending of the bend stiffener
State bend_stiffener::equations(double x, const State &y) {

    State dydx(4);
    double EI = get_EI(x);

    // y[0] = y, y[1] = theta, y[2] = moment, y[3] = shear

    dydx[0] = std::tan(y[1]);
    dydx[1] = y[2] / (EI * std::cos(y[1]));
    dydx[2] = y[3];
    dydx[3] = 0; // point load, depends on the loading

    return dydx;
}

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
                                            target_ML, angle, equations);

    double M0 = result1.first;
    double V0 = result1.second;

    // integrate over the length

    double h = m_length / (double)steps;
    double x = 0.0; // start at the root

    State y = {y0, theta0, M0, V0}; // conditions at the root
    std::vector<State> Results;

    for (int i = 0; i < steps; i++) {
        y = (Integrator::RK4(x, y, h, equations));
        Results.push_back(y);
        x += h;
    };

    return Results;
}
