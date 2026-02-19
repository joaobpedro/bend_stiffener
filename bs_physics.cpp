#include "bs_physics.h"
#include <cmath>
#include <vector>

typedef std::vector<double> State;

double bs_physics::get_EI(double inertia) { return m_E * inertia; }

// define the exact equations for the bending of the bend stiffener
State bs_physics::equations(double x, const State &y, double inertia) {

    State dydx(4);
    double EI = bs_physics::get_EI(inertia);

    // y[0] = y, y[1] = theta, y[2] = moment, y[3] = shear

    dydx[0] = std::tan(y[1]);
    dydx[1] = y[2] / (EI * std::cos(y[1]));
    dydx[2] = y[3];
    dydx[3] = 0; // point load, depends on the loading

    return dydx;
}
