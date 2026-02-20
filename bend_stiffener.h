#ifndef BEND_STIFFENER_H
#define BEND_STIFFENER_H

#include <tuple>
#include <vector>

typedef std::vector<double> State;

class bend_stiffener {
  public:
    // constructor
    bend_stiffener(double root_dia, double length, double tip_dia,
                   double inner_dia);

    ~bend_stiffener();

    // return the variable diameter for the bend stiffner
    double get_dia(double x);

    // return the variable EI for the bend stiffner
    double get_Inertia(double x);

    // calculate volume
    double calculate_volume();

    // calculate strain
    std::vector<State> solve_equations(double tension, double angle);

    // calculate strain from the solved deformations
    State calculate_strain(std::vector<State> Deformations);

    // optmize method here

  private:
    // the bend stiffener is defined by these four variables
    // it is assumed that bend stiffener is a perfect cone
    double m_root_dia;  // meters
    double m_length;    // meters
    double m_tip_dia;   // meters
    double m_inner_dia; // meters

    // this stores the solution state from one angle iteration to the other
    State m_strain = State(1000, 0);
    // State m_initial_conditions = State(3, 0);
};

#endif // BEND_STIFFNER_H
