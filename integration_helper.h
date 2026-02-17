#ifndef INTEGRATION_HELPER_H
#define INTEGRATION_HELPER_H

#include <functional>
#include <vector>

typedef std::vector<double> State;

/*
 * This is a collection of methods for integration
 * it is a very specific problem:
 *  - Cantelever beam integration, this is an assumption
 *  for the bend stiffner condition
 * Zero displacment and curvature at the root
 * Zero moment and shear force are the tip
 */

class Integrator {
public:
  // // construction
  // Integrator();
  // ~Integrator();

  // define the equations for the bend stiffner bending - high angle bending
  State equations(double x, const State &y);

  // Runge Kutta integration method, this is the one I know who to do
  State RK4(double x, const State &y, double h,
            const std::function<State(double, const State &)> &f);

  // the shooting method, basically tries different values of inputs to guess M0
  // and V0
  State shoot(double length, int steps, double y0, double theta0,
              double guessed_M0, double guessed_V0);

  // solving the secant method to find the V0 from the M0
  double solve_V0(double length, int steps, double y0, double theta0,
                  double curr_guessed_M0, double target_theta_L);

  // solves the boundary value problem
  std::pair<double, double> solve_tapered_bvp(double length, int steps,
                                              double y0, double theta0,
                                              double target_ML,
                                              double target_thetaL);

private:
};

#endif
