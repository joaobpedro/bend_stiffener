#include "integration_helper.h"
#include <cmath>
#include <functional>
#include <vector>

typedef std::vector<double> State;
typedef std::function<State(double, const State &, double, double)> Equations;

// runge kutta solver, because its what I know to implement
State Integrator::RK4(double x, const State &y, double h, const Equations &f,
                      double inertia, double strain) {
    size_t n = y.size(); // NOTE: changed type from int to size_t because of
                         // compiler complains
    State k1 = f(x, y, inertia, strain);
    State y_temp(n);

    for (int i = 0; i < n; i++) {
        y_temp[i] = y[i] + 0.5 * h * k1[i];
    }
    State k2 = f(x + 0.5 * h, y_temp, inertia, strain);

    for (int i = 0; i < n; i++) {
        y_temp[i] = y[i] + 0.5 * h * k2[i];
    }
    State k3 = f(x + 0.5 * h, y_temp, inertia, strain);

    for (int i = 0; i < n; i++) {
        y_temp[i] = y[i] + h * k3[i];
    }
    State k4 = f(x + h, y_temp, inertia, strain);

    State y_next(n);

    for (int i = 0; i < n; i++) {
        y_next[i] =
            y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    return y_next;
}

// the shooting function
State Integrator::shoot(double length, int steps, double y0, double theta0,
                        double guessed_M0, double guessed_V0,
                        const Equations &equations, double inertia,
                        double strain) {

    double h = length / (double)steps;
    double x = 0.0; // starting at the root

    State y = {y0, theta0, guessed_M0,
               guessed_V0}; // inital boundary conditions
    //
    for (int i = 0; i < steps; i++) {
        y = RK4(x, y, h, equations, inertia, strain);
        x += h;
    }

    return y;
}

// inner secant: finds the shear[0] to hit the target theta(L)

double Integrator::solve_V0(double length, int steps, double y0, double theta0,
                            double curr_guessed_M0, double target_theta_L,
                            const Equations &equations, double inertia,
                            double strain) {

    double v0 = 100000.0;
    double v1 = -100000.0; // random initial values

    double f0 = shoot(length, steps, y0, theta0, curr_guessed_M0, v0, equations,
                      inertia, strain)[1] -
                target_theta_L;
    double f1 = shoot(length, steps, y0, theta0, curr_guessed_M0, v1, equations,
                      inertia, strain)[1] -
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
        temp_state = shoot(length, steps, y0, theta0, curr_guessed_M0, v1, equations,
                   inertia, strain);
        f1 = temp_state[1] - target_theta_L;
    }

    return v1;
}

// finds M(0) to it the M(L) target, M(L) = 0, cantelever beam

std::pair<double, double>
Integrator::solve_tapered_bvp(double length, int steps, double y0,
                              double theta0, double target_ML,
                              double target_thetaL, const Equations &equations,
                              double inertia, double strain) {

    double u0 = -10000.0;
    double u1 = 10000.0;

    double correct_v0_for_u0 =
        solve_V0(length, steps, y0, theta0, u0, target_thetaL, equations,
                 inertia, strain);
    double f0 = shoot(length, steps, y0, theta0, u0, correct_v0_for_u0,
                      equations, inertia, strain)[2] -
                target_ML;

    double correct_v0_for_u1 =
        solve_V0(length, steps, y0, theta0, u1, target_thetaL, equations,
                 inertia, strain);
    double f1 = shoot(length, steps, y0, theta0, u1, correct_v0_for_u1,
                      equations, inertia, strain)[2] -
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

        correct_v0_for_u1 = solve_V0(length, steps, y0, theta0, u1,
                                     target_thetaL, equations, inertia, strain);
        f1 = shoot(length, steps, y0, theta0, u1, correct_v0_for_u1, equations,
                   inertia, strain)[2] -
             target_ML;
    }
    return {u1, correct_v0_for_u1};
}
