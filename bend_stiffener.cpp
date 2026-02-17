#include "bend_stiffener.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>
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

  return EI_base - (EI_base - EI_tip) * (x / 6.0); // TODO: length hardcoded for now but need to enter here
}

// define the exact equations for the bending of the bend stiffener

State equations(double x, const State &y) {

  State dydx(4);
  double EI = get_IE(x);

  // y[0] = y, y[1] = theta, y[2] = moment, y[3] = shear

  dydx[0] = std::tan(y[1]);
  dydx[1] = y[2] / (EI * std::cos(y[1]));
  dydx[2] = y[3];
  dydx[3] = 0; // point load, depends on the loading

  return dydx;
}

// runge kutta solver, because its what I know to implement
State RK4(double x, const State &y, double h,
          const std::function<State(double, const State &)> &f) {
  size_t n = y.size(); // TODO: changed type from int to size_t because of compiler complains
  State k1 = f(x, y);
  State y_temp(n);

  for (int i = 0; i < n; i++) {
    y_temp[i] = y[i] + 0.5 * h * k1[i];
  }
  State k2 = f(x + 0.5 * h, y_temp);

  for (int i = 0; i < n; i++) {
    y_temp[i] = y[i] + 0.5 * h * k2[i];
  }
  State k3 = f(x + 0.5 * h, y_temp);

  for (int i = 0; i < n; i++) {
    y_temp[i] = y[i] + h * k3[i];
  }
  State k4 = f(x + h, y_temp);

  State y_next(n);

  for (int i = 0; i < n; i++) {
    y_next[i] = y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }

  return y_next;
}

// the shooting function
State shoot(double length, int steps, double y0, double theta0,
            double guessed_M0, double guessed_V0) {

  double h = length / (double)steps;
  double x = 0.0; // starting at the root

  State y = {y0, theta0, guessed_M0, guessed_V0}; // inital boundary conditions
  //
  for (int i = 0; i < steps; i++) {
    y = RK4(x, y, h, equations);
    x += h;
  }

  return y;
}

// inner secant: finds the shear[0] to hit the target theta(L)

double solve_V0(double length, int steps, double y0, double theta0,
                double curr_guessed_M0, double target_theta_L) {

  double v0 = 0;
  double v1 = -10.0; // random initial values

  double f0 =
      shoot(length, steps, y0, theta0, curr_guessed_M0, v0)[1] - target_theta_L;
  double f1 =
      shoot(length, steps, y0, theta0, curr_guessed_M0, v1)[1] - target_theta_L;

  for (int i = 0; i < 100; i++) {
    if (std::abs(f1) < 1e-6)
      return v1;
    if (std::abs(f1 - f0) < 1e-12)
      break;

    double v2 = v1 - f1 * (v1 - v0) / (f1 - f0);
    v0 = v1;
    f0 = f1;
    v1 = v2;
    f1 = shoot(length, steps, y0, theta0, curr_guessed_M0, v1)[1] -
         target_theta_L;
  }

  return v1;
}

// finds M(0) to it the M(L) target, M(L) = 0, cantelever beam

std::pair<double, double> solve_tapered_bvp(double length, int steps, double y0,
                                            double theta0, double target_ML,
                                            double target_thetaL) {

  double u0 = 0.0;
  double u1 = -100.0;

  double correct_v0_for_u0 =
      solve_V0(length, steps, y0, theta0, u0, target_thetaL);
  double f0 =
      shoot(length, steps, y0, theta0, u0, correct_v0_for_u0)[2] - target_ML;

  double correct_v0_for_u1 =
      solve_V0(length, steps, y0, theta0, u1, target_thetaL);
  double f1 =
      shoot(length, steps, y0, theta0, u1, correct_v0_for_u1)[2] - target_ML;

  for (int i = 0; i < 100; i++) {
    if (std::abs(f1) < 1e-6)
      return {u1, correct_v0_for_u1};
    if (std::abs(f1 - f0) < 1e-12)
      break;

    double u2 = u1 - f1 * (u1 - u0) / (f1 - f0);
    u0 = u1;
    f0 = f1;
    u1 = u2;

    correct_v0_for_u1 = solve_V0(length, steps, y0, theta0, u1, target_thetaL);
    f1 = shoot(length, steps, y0, theta0, u1, correct_v0_for_u1)[2] - target_ML;
  }
  return {u1, correct_v0_for_u1};
}

std::vector<State> bend_stiffener::calculate_strain(double tension, double angle) {

  // boundary conditions based on cantelever beam modenl
  double y0 = 0;
  double theta0 = 0;
  double target_ML = 0;

  int steps = 1000;

  // calculate the moment and shear force at the root
  std::pair<double, double> result1;
  result1 =
      solve_tapered_bvp(m_length, steps, y0, theta0, target_ML, angle);

  double M0 = result1.first;
  double V0 = result1.second;

  // integrate over the length

  double h = m_length / (double)steps;
  double x = 0.0; // start at the root

  State y = {y0, theta0, M0, V0}; // conditions at the root
  std::vector<State> Results;

  for (int i = 0; i < steps; i++) {
    y = (RK4(x, y, h, equations));
    Results.push_back(y);
    x += h;
  };

  return Results;
}
