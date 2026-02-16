#include "bend_stiffener.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

typedef std::tuple<std::vector<double>, std::vector<double>,
                   std::vector<double>>
    solution;
typedef std::vector<double> State;

// functions prototypes

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

State derivs(double x, const State &y) {
  State dydx(4);
  dydx[0] = y[1];
  dydx[1] = y[2];
  dydx[2] = y[3];
  dydx[3] = -10; // this is a dummy uniform loadind.

  return dydx;
};

State rk4_step(double x, const State &Y, double h,
               const std::function<State(double, const State &)> &f) {
  int n = Y.size();
  State k1 = f(x, Y);

  State Y_temp(n);
  for (int i = 0; i < n; ++i)
    Y_temp[i] = Y[i] + 0.5 * h * k1[i];
  State k2 = f(x + 0.5 * h, Y_temp);

  for (int i = 0; i < n; ++i)
    Y_temp[i] = Y[i] + 0.5 * h * k2[i];
  State k3 = f(x + 0.5 * h, Y_temp);

  for (int i = 0; i < n; ++i)
    Y_temp[i] = Y[i] + h * k3[i];
  State k4 = f(x + h, Y_temp);

  State Y_next(n);
  for (int i = 0; i < n; ++i) {
    Y_next[i] = Y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }
  return Y_next;
};

double shoot(double L, int steps, double y0, double yp0, double ypp0,
             double guessed_yppp0) {
  double h = L / steps;
  double x = 0.0;

  // Initial state: [y(0), y'(0), y''(0), guessed y'''(0)]
  State Y = {y0, yp0, ypp0, guessed_yppp0};

  for (int i = 0; i < steps; ++i) {
    Y = rk4_step(x, Y, h, derivs);
    x += h;
  }

  // We are trying to hit a specific SLOPE at x = L, so we return Y[1]
  return Y[1];
}

double solve_bvp(double L, double y0, double yp0, double ypp0,
                 double target_ypL, int steps) {
  // Two initial guesses for the unknown initial 3rd derivative y'''(0)
  double u0 = 0.0;
  double u1 = 1.0;

  double tol = 1e-6;
  int max_iter = 100;

  // Evaluate the error (resulting final slope minus target final slope)
  double f0 = shoot(L, steps, y0, yp0, ypp0, u0) - target_ypL;
  double f1 = shoot(L, steps, y0, yp0, ypp0, u1) - target_ypL;

  for (int iter = 0; iter < max_iter; ++iter) {
    if (std::abs(f1) < tol) {
      std::cout << "Converged in " << iter << " iterations." << std::endl;
      return u1; // We found the correct y'''(0)
    }

    if (std::abs(f1 - f0) < 1e-12) {
      std::cerr << "Secant method failed: division by zero." << std::endl;
      break;
    }

    // Secant calculation
    double u2 = u1 - f1 * (u1 - u0) / (f1 - f0);

    u0 = u1;
    f0 = f1;
    u1 = u2;
    f1 = shoot(L, steps, y0, yp0, ypp0, u1) - target_ypL;
  }

  std::cerr << "Warning: Maximum iterations reached without convergence."
            << std::endl;
  return u1;
}

void bend_stiffener::calculate_strain(double tension, double angle) {

  // dedfine integration params
  double L = 10.0;          // Domain length
  double y0 = 0.0;          // Known value at x = 0
  double yp0 = 0.0;         // Known slope at x = 0
  double ypp0 = 5.0;        // Known curvature at x = 0
  double target_ypL = -2.0; // Known target SLOPE at x = L
  int steps = 1000;

  double correct_yppp0 = solve_bvp(L, y0, yp0, ypp0, target_ypL, steps);

  // Step 2: Integrate one last time to output the true trajectory
  double h = L / steps;
  double x = 0.0;
  State Y = {y0, yp0, ypp0, correct_yppp0};

  std::cout << "Final Trajectory (sampled):" << std::endl;
  for (int i = 0; i <= steps; ++i) {
    if (i % 250 == 0) {
      std::cout << "x = " << x << " | y = " << Y[0] << " | y' = " << Y[1]
                << " | y'' = " << Y[2] << " | y''' = " << Y[3] << std::endl;
    }
    if (i < steps) {
      Y = rk4_step(x, Y, h, derivs);
      x += h;
    }
  }
};
