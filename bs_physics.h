#ifndef BS_PHYSICS_H
#define BS_PHYSICS_H

#include <tuple>
#include <vector>

typedef std::vector<double> State;

class bs_physics {
  public:
    // constructor

    // return the variable EI for the bend stiffner
    static double get_EI(double inertia);

    // Eqautions defining the behavior
    static State equations(double x, const State &y, double inertia);

  private:
    // material properties
    static constexpr double m_E =
        1000e6; // 1000MPa for 60D PU constant PU material for now
};

#endif // BEND_STIFFNER_H
