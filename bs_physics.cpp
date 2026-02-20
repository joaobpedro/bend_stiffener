#include "bs_physics.h"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

typedef std::vector<double> State;

double bs_physics::get_EI(double inertia, double strain) {

    if (strain == 0) {
        return m_E * inertia;
    } else {
        return bs_physics::get_non_linear_E(strain) * inertia;
    }
};

// define the exact equations for the bending of the bend stiffener
State bs_physics::equations(double x, const State &y, double inertia,
                            double strain) {

    State dydx(4);
    double EI = bs_physics::get_EI(inertia, strain);

    // y[0] = y, y[1] = theta, y[2] = moment, y[3] = shear

    dydx[0] = std::tan(y[1]);
    dydx[1] = y[2] / (EI * std::cos(y[1]));
    dydx[2] = y[3];
    dydx[3] = 0; // point load, depends on the loading

    return dydx;
}

double bs_physics::get_non_linear_E(double strain) {
    std::map<double, double> non_linear_E = {
        {0, 215800},   {0.01, 215800}, {0.02, 179500}, {0.03, 138100},
        {0.04, 94300}, {0.05, 59500},  {0.06, 39700},  {0.07, 28200},
        {0.08, 21100}, {0.09, 17000},  {0.1, 14300},   {0.11, 12700},
        {0.12, 11000}, {0.13, 10600},  {0.14, 10000},  {0.15, 9700},
        {0.16, 9600},  {0.17, 8600},   {0.18, 9400},   {0.19, 8600},
        {0.2, 8500}

    };

    // if strain is larger than the curve raise an error

    // finds the key that its larger or equal to my strain
    auto key_up = non_linear_E.lower_bound(strain);

    if (key_up == non_linear_E.begin()) {
        std::cout
            << "Strain is out of the lower bound, ZERO strain will be used"
            << std::endl;
        return key_up->second;
    }

    if (key_up == non_linear_E.end()) {
        std::cout << "Strain is out of the upper bound, MAX strain will be used"
                  << std::endl;
        return key_up->second;
    }

    if (key_up->first == strain) {
        return key_up->second;
    }

    // get the previous key from the key_up
    auto key_low = std::prev(key_up);

    double x0 = key_low->first;
    double x1 = key_up->first;
    double y0 = key_low->second;
    double y1 = key_up->second;
    double result = y1 + (strain - x0) * (y0 - y1) / (x1 - x0);
    return result * 1000.0 * 1000.0; // convert to MPa
}
