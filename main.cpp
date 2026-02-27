#include <iostream>
using namespace std;

#include "bs_engine.h"

typedef std::vector<double> State;


int main(int argc, char **argv) {
    // TODO: insert code
    // l = 6m
    // RD = 0.65
    // TD = 0.232
    // ID = 0.218

    Dimensions bs_jotun;

    bs_jotun.length = 6.0; // m
    bs_jotun.root_dia = 0.65; // m
    bs_jotun.tip_dia = 0.232; // m
    bs_jotun.inner_dia = 0.218; // m

    // initial conditions

    double thetaL = 0.1;
    double theta = 0.0;
    double theta0 = 0.0;
    double y0 = 0;
    double mL = 0.0;

    std::pair<double, double> results0;

    size_t steps = 1000;
    size_t iterations = 100;

    std::vector<double> strain(steps, 0); // initilize the strain with zeros
    std::vector<State> Results(steps, State(4));

    while (theta < thetaL) {
        for (int I = 0; I < iterations; I++) {

            std::vector<State> temp_results;
            results0 = solve_tapered_bvp(bs_jotun, steps, y0, theta0, mL, theta, strain); // initial strain is zero.

            double M0 = results0.first;
            double V0 = results0.second;
            State y = {y0, theta0, M0, V0};
            double x = 0.0;
            double h = bs_jotun.length / (double) steps;

            for (int J = 0; J < steps; J++) {
                y = RK4(x, y, h, bs_jotun, strain[J]);
                temp_results.push_back(y);
                x += h;
            }
            // double error = temp_results[0][2] - Results[0][2]; // error measured in the theta

            //update strain
            std::vector<double> prev_strain = strain;
            calculate_strain(steps, bs_jotun, temp_results, strain);
            // store the Results
            Results = temp_results;

            double error = strain[0] - prev_strain[0];
            std::cout << "Iteration no:" << I << ". Error: " << error << std::endl;
            // stops if error less than 1e-6
            if (std::abs(error) < 1e-6) {
                break;
            }

        }
        std::cout << "Applied theta: " << theta << std::endl;
        theta += 0.01;
    }

    std::cout << "Length Deformation Theta Moment Shear" << std::endl;
    double h = bs_jotun.length / steps;
    for (auto item :  Results) {
        std::cout << h << " " <<  item[0] << " " << item[1] << " " << item[2] << " " << item[3] << std::endl;
        h += bs_jotun.length/steps;
    }

    std::cout << "Strain" << std::endl;
    for (auto item : strain) {
        std::cout << item << std::endl;
    }

    return 0;
}
