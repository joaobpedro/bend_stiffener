#ifndef BS_ENGINE_H_
#define BS_ENGINE_H_

#include <map>
#include <vector>

typedef std::vector<double> State;

struct Dimensions {
    double length;
    double root_dia;
    double tip_dia;
    double inner_dia;
};

 
double get_non_linear_E(double strain) {
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
        // std::cout
        //     << "Strain is out of the lower bound, ZERO strain will be used"
        //     << std::endl;
        return key_up->second * 1000.0;
    }

    if (key_up == non_linear_E.end()) {
        // std::cout << "Strain is out of the upper bound, MAX strain will be
        // used"
        //           << std::endl;
        return std::prev(key_up)->second* 1000.0;
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
    return result * 1000.0; // convert to MPa
};

double get_dia(Dimensions& dimensions, double x) {

    double outer_dia = dimensions.root_dia - (dimensions.root_dia - dimensions.tip_dia) * 
        (x/ dimensions.length);
    return outer_dia;
}

double get_Inertia(Dimensions& dimensions, double x) {

    const double pi = 3.141592653599;
    // calculate outer diameter
    // this assumes the conic shape
    double outer_dia = dimensions.root_dia - (dimensions.root_dia - dimensions.tip_dia) * 
        (x/ dimensions.length);

    // inertia is given by  I = pi/64 *(Do^4 - Di^4)

    double inertia = (pi / 64) * (pow(outer_dia, 4) - pow(dimensions.inner_dia, 4));
    return inertia;
}

double get_EI(Dimensions& dimensions, double strain, double x) {
    
    const double m_E = 1100000;
    double inertia = get_Inertia(dimensions, x);
    if (strain == 0) {
        return m_E * inertia;
    } else {
        return get_non_linear_E(strain) * inertia;
    }
};


State equations(double x, const State &y, Dimensions& dimensions,
                            double strain) {

    State dydx(4);
    // double inertia_x = bend_stiffener::get_Inertia(x);
    double EI = get_EI(dimensions, strain, x);

    // y[0] = y, y[1] = theta, y[2] = moment, y[3] = shear

    dydx[0] = std::tan(y[1]);
    dydx[1] = y[2] / (EI * std::cos(y[1]));
    dydx[2] = y[3];
    dydx[3] = 0; // point load, depends on the loading

    return dydx;
}


State RK4(double x, const State &y, double h, Dimensions& dimensions, double strain) {
    size_t n = y.size(); // NOTE: changed type from int to size_t because of
                         // compiler complains
    State k1 = equations(x, y, dimensions, strain);
    State y_temp(n);

    for (int i = 0; i < n; i++) {
        y_temp[i] = y[i] + 0.5 * h * k1[i];
    }
    State k2 = equations(x + 0.5 * h, y_temp, dimensions, strain);

    for (int i = 0; i < n; i++) {
        y_temp[i] = y[i] + 0.5 * h * k2[i];
    }
    State k3 = equations(x + 0.5 * h, y_temp, dimensions, strain);

    for (int i = 0; i < n; i++) {
        y_temp[i] = y[i] + h * k3[i];
    }
    State k4 = equations(x + h, y_temp, dimensions, strain);

    State y_next(n);

    for (int i = 0; i < n; i++) {
        y_next[i] =
            y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    return y_next;
}


State shoot(Dimensions& dimensions, int steps, double y0, double theta0, double guessed_M0, double guessed_V0, std::vector<double>& strain) {

    double h = dimensions.length / (double)steps;
    double x = 0.0; // starting at the root

    State y = {y0, theta0, guessed_M0,
               guessed_V0}; // inital boundary conditions

    for (int i = 0; i < steps; i++) {
        double strain_d = strain[i];
        y = RK4(x, y, h, dimensions, strain_d);
        x += h;
    }

    return y;
}


double solve_V0(Dimensions& dimensions, int steps, double y0, double theta0,
                            double curr_guessed_M0, double target_theta_L,
                            std::vector<double>& strain) {

    double v0 = 100.0;
    double v1 = -100.0; // random initial values

    double f0 = shoot(dimensions, steps, y0, theta0, curr_guessed_M0, v0, strain)[1] -
                target_theta_L;
    double f1 = shoot(dimensions, steps, y0, theta0, curr_guessed_M0, v1,
                      strain)[1] -
                target_theta_L;
    // for degbug
    State temp_state;
    for (int i = 0; i < 100; i++) {
        if (std::abs(f1) < 1e-6)
            return v1;
        if (std::abs(f1 - f0) < 1e-12)
            break;

        double v2 = v1 - f1 * (v1 - v0) / (f1 - f0);
        v0 = v1;
        f0 = f1;
        v1 = v2;
        temp_state = shoot(dimensions, steps, y0, theta0, curr_guessed_M0, v1,
                   strain);
        f1 = temp_state[1] - target_theta_L;
    }

    return v1;
}

std::pair<double, double>
solve_tapered_bvp(Dimensions& dimensions, int steps, double y0,
                              double theta0, double target_ML,
                              double target_thetaL,
                               std::vector<double>& strain) {

    double u0 = -100.0;
    double u1 = 100.0;

    double correct_v0_for_u0 =
        solve_V0(dimensions, steps, y0, theta0, u0, target_thetaL, strain);
    double f0 = shoot(dimensions, steps, y0, theta0, u0, correct_v0_for_u0,
                      strain)[2] -
                target_ML;

    double correct_v0_for_u1 =
        solve_V0(dimensions, steps, y0, theta0, u1, target_thetaL, strain); 
    double f1 = shoot(dimensions, steps, y0, theta0, u1, correct_v0_for_u1,
                      strain)[2] -
                target_ML;

    for (int i = 0; i < 100; i++) {
        if (std::abs(f1) < 1e-6)
            return {u1, correct_v0_for_u1};
        if (std::abs(f1 - f0) < 1e-12)
            break;

        double u2 = u1 - f1 * (u1 - u0) / (f1 - f0);
        u0 = u1;
        f0 = f1;
        u1 = u2;

        correct_v0_for_u1 = solve_V0(dimensions, steps, y0, theta0, u1,
                                     target_thetaL, strain);
        f1 = shoot(dimensions, steps, y0, theta0, u1, correct_v0_for_u1, strain)[2] -
             target_ML;
    }
    return {u1, correct_v0_for_u1};
}

void calculate_strain(size_t steps, Dimensions& dimensions, State& strain, std::vector<State> Deformations) {

    // curvature is d_theta/d_s
    double h = dimensions.length / (double) steps; // HARDCODED: discretized
    double x = 0.0;

    for (int i = 0; i < Deformations.size(); i++) {
        double EI = get_EI(dimensions, strain[i], x); 
        strain[i] = (get_dia(dimensions, x) * Deformations[i][2]) / (EI);
        x += h;
    }
}
#endif
