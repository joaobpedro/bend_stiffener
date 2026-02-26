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
//     std::vector<State> to_plot;
//
//     // this was the test
//     bend_stiffener Jotun(0.65, 6.0, 0.232, 0.218);
//     // to_plot =
//     //     Jotun.solve_equations(0.2, 0.1); // right now these inputs are dummy
//     //
//     // std::cout << "Deformations Theta Moment Shear" << std::endl;
//     // for (int i = 0; i < to_plot.size(); i++) {
//     //     std::cout << to_plot[i][0] << ' ' << to_plot[i][1] << ' '
//     //               << to_plot[i][2] << ' ' << to_plot[i][3] << std::endl;
//     // }
//     //
//     // std::vector<double> strain = Jotun.calculate_strain(to_plot);
//     //
//     // this will apply the angles in a progressive way
//
//     double max_angle = 0.1; // this is 0.2 rad
//     // right now we do not do anything with the tension
//
//     std::vector<double> Strain(1000, 0);
//     std::vector<double> Prev_Strain(1000, 0);
//     std::vector<State> deformations;
//
//     // solve for maximum angle first
//     double error = 10; // initialise error high
// //    while (error > 1e-6) {
// //        deformations = Jotun.solve_equations(1, max_angle);
// //        Prev_Strain = Strain; // store last strain
// //        // update the current strain
// //        Strain = Jotun.calculate_strain(deformations);
// //
// //        double error1 = std::abs(Strain[0] - Prev_Strain[0]);
// //        double error2 = std::abs(Strain[500] - Prev_Strain[500]);
// //        double error3 = std::abs(Strain[999] - Prev_Strain[999]);
// //
// //        error = std::max({error1, error2, error3});
// //
// //        std::cout << "Convergence Error: " << error << std::endl;
// //    }
//
//     deformations = Jotun.solve_equations(1, max_angle);
//     Prev_Strain = Strain; // store last strain
//     // update the current strain
//     Strain = Jotun.calculate_strain(deformations);
//
//     double error1 = std::abs(Strain[0] - Prev_Strain[0]);
//     double error2 = std::abs(Strain[500] - Prev_Strain[500]);
//     double error3 = std::abs(Strain[999] - Prev_Strain[999]);
//
//     error = std::max({error1, error2, error3});
//
//     std::cout << "Convergence Error: " << error << std::endl;
//
//     std::cout << "Analysis Finished with Convergence Error: " << error
//               << std::endl;
//
//     for (int i = 0; i < deformations.size(); i++) {
//         std::cout << deformations[i][0] << ' ' << deformations[i][1] << ' '
//                   << deformations[i][2] << ' ' << deformations[i][3]
//                   << std::endl;
//     }
//
//     std::cout << "Strain Results:" << std::endl;
//     for (int i = 0; i < Strain.size(); i++) {
//         std::cout << Strain[i] << std::endl;
//     }

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
            calculate_strain(steps, bs_jotun, strain, temp_results);
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
