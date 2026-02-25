#ifndef BS_PHYSICS_H
#define BS_PHYSICS_H

#include <tuple>
#include <vector>

typedef std::vector<double> State;

struct Dimensions {
    double length;
    double root_diameter;
    double tip_diameter;
    double inner_diameter;
};

class bs_physics {
  public:
    // return the I, based on diameter, length
    static double get_Inertia(Dimensions& dimensions, double x);

    // return the variable EI for the bend stiffner
    static double get_EI(Dimensions& dimensions, double strain, double x);

    // Eqautions defining the behavior
    static State equations(double x, const State &y, Dimensions& dimensions,
                           double strain);

    // returns the non-linear value for the E-mod depending on the strain
    static double get_non_linear_E(double strain);

  private:
    // material properties
    static constexpr double m_E =
        1000e6; // 1000MPa for 60D PU constant PU material for now
};

#endif // BEND_STIFFNER_H
