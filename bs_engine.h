#ifndef BS_ENGINE_H_
#define BS_ENGINE_H_

#include <map>
#include <vector>
#include <math.h>

#define DISCRETIZATION 1000
#define NUMBER_STATES 4

typedef std::vector<double> State;

typedef struct Dimensions {
    double length;
    double root_dia;
    double tip_dia;
    double inner_dia;
} bend_stiffner;

// if I am doing this super performant I need to re-write
// 1. maps - I will just get two arrays and make them the same size
// 2. vectors

double get_non_linear_E(double strain) {
    
    const int material_data_size = 21;
    
    double strain_vector[material_data_size] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20};
    double Emod_vector[material_data_size] = {215800, 215800, 179500, 138100, 94300, 59500, 39700, 28200, 21100, 17000, 14300, 12700, 11000, 10600, 10000, 9700, 9600, 8600, 9400, 8600, 8500};
    
    double distance = 1000;
    int min_index = 0;
    int max_index = 0;
    
    // need to use the absolute values here
    for (int i = 0; i < sizeof(strain); i++)
    {
        double prev_distance = distance;
        distance = strain_vector[i] - strain;
        if (fabs(distance < prev_distance)) 
        {
            if (distance < 0)
            {
                min_index = i;
                max_index = i + 1;
            };
            if (distance > 0)
            {
                min_index = i - 1;
                max_index = 1;
            }
        };
    }
    
    // need to deal with the edge cases
    if (min_index < 0) {
        min_index = 1;
        max_index = min_index + 1;
    };
    
    if (max_index > sizeof(strain_vector))
    {
        max_index = sizeof(strain_vector);
        min_index = max_index - 1;
    }
    
    double x0 = strain_vector[min_index];
    double x1 = strain_vector[max_index];
    double y0 = Emod_vector[min_index];
    double y1 = Emod_vector[max_index];
    double result = y1 + (strain - x0) * (y0 - y1) / (x1 - x0);
    return result * 1000.0; // convert to MPa
};

double get_dia(const Dimensions& dimensions, double x) {

    double outer_dia = dimensions.root_dia - (dimensions.root_dia - dimensions.tip_dia) * 
        (x/ dimensions.length);
    return outer_dia;
}

double get_Inertia(const Dimensions& dimensions, double x) {

    const double pi = 3.141592653599;
    // calculate outer diameter
    // this assumes the conic shape
    double outer_dia = dimensions.root_dia - (dimensions.root_dia - dimensions.tip_dia) * 
        (x/ dimensions.length);

    // inertia is given by  I = pi/64 *(Do^4 - Di^4)

    double inertia = (pi / 64) * (pow(outer_dia, 4) - pow(dimensions.inner_dia, 4));
    return inertia;
}

double get_EI(const Dimensions& dimensions, double strain, double x) {
    
    const double m_E = 1100000;
    double inertia = get_Inertia(dimensions, x);
    if (strain == 0) {
        return m_E * inertia;
    } else {
        return get_non_linear_E(strain) * inertia;
    }
};


State equations(double x, const State &y, const Dimensions& dimensions,
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


State RK4(double x, const State &y, double h, const Dimensions& dimensions, double strain) {
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


State shoot(const Dimensions& dimensions, size_t steps, double y0, double theta0, double guessed_M0, double guessed_V0, const std::vector<double>& strain) {

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


double solve_V0(const Dimensions& dimensions, size_t steps, double y0, double theta0,
                            double curr_guessed_M0, double target_theta_L, double V0_guess,
                            const std::vector<double>& strain) {

    double v0 = -V0_guess;
    double v1 = V0_guess; // random initial values

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
solve_tapered_bvp(const Dimensions& dimensions, size_t steps, double y0,
                              double theta0, double target_ML,
                              double target_thetaL, double M0_guess, double V0_guess,
                              const std::vector<double>& strain) {

    double u0 = -M0_guess;
    double u1 = M0_guess;

    double correct_v0_for_u0 =
        solve_V0(dimensions, steps, y0, theta0, u0, target_thetaL, V0_guess, strain);
    double f0 = shoot(dimensions, steps, y0, theta0, u0, correct_v0_for_u0,
                      strain)[2] -
                target_ML;

    double correct_v0_for_u1 =
        solve_V0(dimensions, steps, y0, theta0, u1, target_thetaL, V0_guess, strain); 
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
        V0_guess = correct_v0_for_u1;
        correct_v0_for_u1 = solve_V0(dimensions, steps, y0, theta0, u1,
                                     target_thetaL, V0_guess, strain);
        f1 = shoot(dimensions, steps, y0, theta0, u1, correct_v0_for_u1, strain)[2] -
             target_ML;
    }
    return {u1, correct_v0_for_u1};
}

void calculate_strain(size_t steps, const Dimensions& dimensions, const std::vector<State>& Deformations, State& strain)  {

    // curvature is d_theta/d_s
    double h = dimensions.length / (double) steps; 
    double x = 0.0;

    for (int i = 0; i < Deformations.size(); i++) {
        double EI = get_EI(dimensions, strain[i], x); 
        strain[i] = (get_dia(dimensions, x) * Deformations[i][2]) / (EI);
        x += h;
    }
}
#endif
