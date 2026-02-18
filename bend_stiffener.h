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

    // return the variable EI for the bend stiffner
    static double get_EI(double x);

    // calculate volume
    double calculate_volume();

    // double get_IE(double length);

    // Eqautions defining the behavior
    static State equations(double x, const State &y);

    // calculate strain
    std::vector<State> calculate_strain(double tension, double angle);

    // optmize
    // optmize method here

  private:
    // the bend stiffener is defined by these four variables
    // it is assumed that bend stiffener is a perfect cone
    double m_root_dia;  // meters
    double m_length;    // meters
    double m_tip_dia;   // meters
    double m_inner_dia; // meters
};

#endif // BEND_STIFFNER_H
