#include <iostream>
using namespace std;

#include <wx/app.h>
#include <wx/defs.h>
#include <wx/platinfo.h> // platform info

#include <wx/cmdline.h> // command line parser
#include <wx/file.h>    // wxFile
#include <wx/filefn.h>  // File Functions
#include <wx/string.h>  // wxString

#include "bs_engine.h"

typedef map<wxString, wxString> CmdLineMap; // CmdLineMap
typedef std::vector<double> State;

/*
   switch 	  This is a boolean option which can be given or not, but which
   doesn't have any value. We use the word switch to distinguish such boolean
   options from more generic options like those described below. For example, -v
   might be a switch meaning "enable verbose mode".

   option 	  Option for us here is something which comes with a value
   unlike a switch. For example, -o:filename might be an option which allows to
   specify the name of the output file.

   parameter  This is a required program argument.

   More info at:
   http://docs.wxwidgets.org/2.8/wx_wxcmdlineparser.html#wxcmdlineparser

*/

static const wxCmdLineEntryDesc cmdLineDesc[] = {
    //   kind           shortName      longName         description
    //   parameterType          flag(s)
    {wxCMD_LINE_SWITCH, wxT_2("h"), wxT_2("help"),
     wxT_2("command line parameter help"), wxCMD_LINE_VAL_NONE,
     wxCMD_LINE_OPTION_HELP},
    {wxCMD_LINE_SWITCH, wxT_2("v"), wxT_2("verbose"), wxT_2("Dump parameters"),
     wxCMD_LINE_VAL_NONE, wxCMD_LINE_PARAM_OPTIONAL},
    //  { wxCMD_LINE_OPTION, wxT_2("p"),    wxT_2("project"), wxT_2("Project
    //  name"),                wxCMD_LINE_VAL_STRING,
    //  wxCMD_LINE_OPTION_MANDATORY  },
    //  { wxCMD_LINE_OPTION, wxT_2("t"),    wxT_2("target"),  wxT_2("Target
    //  name"),                 wxCMD_LINE_VAL_STRING,
    //  wxCMD_LINE_OPTION_MANDATORY  },
    //  { wxCMD_LINE_OPTION, wxT_2("tdir"), wxT_2("tdir"),    wxT_2("Target
    //  output directory"),     wxCMD_LINE_VAL_STRING,
    //  wxCMD_LINE_OPTION_MANDATORY  },
    {wxCMD_LINE_NONE, wxT_2(""), wxT_2(""), wxT_2(""), wxCMD_LINE_VAL_NONE,
     wxCMD_LINE_PARAM_OPTIONAL}};

// sample command line params: --project=cpde_export  --target=W32_Debug --tdir=w32\bin\Debug\

void ParserToMap(size_t argc, wxCmdLineParser &parser, CmdLineMap &cmdMap) {
    size_t pcount = sizeof(cmdLineDesc) / sizeof(wxCmdLineEntryDesc) - 1;
    if (argc < pcount)
        pcount = argc;
    for (size_t i = 0; i < pcount; i++) {
        wxString pname = cmdLineDesc[i].longName;
        if (cmdLineDesc[i].kind == wxCMD_LINE_PARAM) {
            wxString pvalue = parser.GetParam(i - 1);
            cmdMap.insert(make_pair(pname, pvalue));
        } else {
            // switch or option, mush check if present
            if (parser.Found(pname)) {
                wxString pvalue;
                if (cmdLineDesc[i].type == wxCMD_LINE_VAL_STRING) {
                    parser.Found(pname, &pvalue);
                } else if (cmdLineDesc[i].type == wxCMD_LINE_VAL_NUMBER) {
                    long lvalue = 0;
                    parser.Found(pname, &lvalue);
                    pvalue.Printf(wxT("%i"), lvalue);
                }
                cmdMap.insert(make_pair(pname, pvalue));
            }
        }
    }
}

int main(int argc, char **argv) {
    // initialise wxWidgets library
    wxInitializer initializer(argc, argv);

    // parse command line
    wxCmdLineParser parser(cmdLineDesc);
    parser.SetSwitchChars(wxT("-"));
    parser.SetCmdLine(argc, argv);
    if (parser.Parse() != 0) {
        // command line parameter error
        return 1;
    }

    // parser success
    // convert parameters to map
    CmdLineMap cmdMap;
    ParserToMap(argc, parser, cmdMap);

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

    double tethaL = 0.1;
    double theta = 0.0;
    double theta0 = 0.0;
    double y0 = 0;
    double mL = 0.0;

    std::pair<double, double> results0;

    size_t steps = 1000;
    size_t iterations = 100;

    std::vector<double> strain(steps, 0); // initilize the strain with zeros
    std::vector<State> Results(steps, State(4));

    while (theta < tethaL) {
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
                x =+ h;
            }
            // double error = temp_results[0][2] - Results[0][2]; // error measured in the tetha

            //update strain
            std::vector<double> prev_strain = strain;
            calculate_strain(steps, bs_jotun, strain, temp_results);
            // store the Results
            Results = temp_results;

            double error = strain[0] - prev_strain[0];
            std::cout << "Iteration error: " << error << std::endl;
            // stops if error less than 1e-6
            if (std::abs(error) < 1e-6) {
                break;
            }

        }
        std::cout << "Applied tetha: " << theta << std::endl;
        theta += 0.01;
    }

    std::cout << "Deformation Theta Moment Shear" << std::endl;
    for (auto item :  Results) {
        std::cout << item[0] << " " << item[1] << " " << item[2] << " " << item[3] << std::endl;
    }

    std::cout << "Strain" << std::endl;
    for (auto item : strain) {
        std::cout << item << std::endl;
    }

    return 0;
}
