#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>

// Boost Includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace ublas = boost::numeric::ublas;
using matrix_t = ublas::matrix<double>;
using vector_t = ublas::vector<double>;

const double PI = 3.14159265358979323846;

// ==========================================
// 1. Problem Definition Structures
// ==========================================

// Defines a single equation:
// A*u'' + B*u' + C*u + D*v'' + E*v' + F*v = RHS
struct EquationCoeffs {
    std::function<double(double)> coeff_u_d2; // coeff for u''
    std::function<double(double)> coeff_u_d1; // coeff for u'
    std::function<double(double)> coeff_u;    // coeff for u

    std::function<double(double)> coeff_v_d2; // coeff for v''
    std::function<double(double)> coeff_v_d1; // coeff for v'
    std::function<double(double)> coeff_v;    // coeff for v

    std::function<double(double)> rhs;        // Right Hand Side function
};

// Defines Boundary Condition: a*y + b*y' = c
struct BoundaryCondition {
    double a; // coeff for function value
    double b; // coeff for derivative
    double c; // target value
};

// ==========================================
// 2. The General Solver Class
// ==========================================

class CoupledSolver {
private:
    int N;
    vector_t x;
    matrix_t D, D2;

    // Helper: Generate Chebyshev Points
    void initGrid() {
        x.resize(N + 1);
        for (int i = 0; i <= N; ++i) x[i] = std::cos(PI * i / N);
    }

    // Helper: Generate Differentiation Matrices
    void initDiffMatrices() {
        D = matrix_t(N + 1, N + 1);
        std::vector<double> c(N + 1);
        c[0] = 2.0; c[N] = 2.0;
        for (int i = 1; i < N; ++i) c[i] = 1.0;

        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                if (i != j) {
                    D(i, j) = (c[i] * std::pow(-1.0, i + j)) / (c[j] * (x[i] - x[j]));
                }
            }
        }
        for (int i = 1; i < N; ++i) D(i, i) = -x[i] / (2.0 * (1.0 - x[i] * x[i]));
        D(0, 0) = (2.0 * N * N + 1.0) / 6.0;
        D(N, N) = -D(0, 0);

        D2 = ublas::prod(D, D);
    }

public:
    CoupledSolver(int n_points) : N(n_points) {
        initGrid();
        initDiffMatrices();
    }

    // Helper to get defaults (zeros)
    static double Zero(double) { return 0.0; }

    // Main Solve Function
    // Returns a pair of vectors {u, v}
    std::pair<vector_t, vector_t> solve(
        const EquationCoeffs& eq1,
        const EquationCoeffs& eq2,
        const BoundaryCondition& u_left, const BoundaryCondition& u_right,
        const BoundaryCondition& v_left, const BoundaryCondition& v_right
    ) {
        int num_pts = N + 1;
        int sys_size = 2 * num_pts;

        matrix_t A(sys_size, sys_size, 0.0);
        vector_t b(sys_size, 0.0);

        // Fill Matrix
        for (int i = 0; i < num_pts; ++i) {
            double xi = x[i];
            int row_u = i;
            int row_v = i + num_pts;

            // --- Apply Boundary Conditions ---
            // Note: In Chebyshev, i=0 is x=1 (Right), i=N is x=-1 (Left)

            // RIGHT Boundary (x=1, i=0)
            if (i == 0) {
                // BC for u at Right
                for(int j=0; j<num_pts; ++j) {
                    A(row_u, j) = u_right.a * (i==j?1.0:0.0) + u_right.b * D(i, j);
                }
                b(row_u) = u_right.c;

                // BC for v at Right
                for(int j=0; j<num_pts; ++j) {
                    A(row_v, j + num_pts) = v_right.a * (i==j?1.0:0.0) + v_right.b * D(i, j);
                }
                b(row_v) = v_right.c;
                continue;
            }

            // LEFT Boundary (x=-1, i=N)
            if (i == N) {
                // BC for u at Left
                for(int j=0; j<num_pts; ++j) {
                    A(row_u, j) = u_left.a * (i==j?1.0:0.0) + u_left.b * D(i, j);
                }
                b(row_u) = u_left.c;

                // BC for v at Left
                for(int j=0; j<num_pts; ++j) {
                    A(row_v, j + num_pts) = v_left.a * (i==j?1.0:0.0) + v_left.b * D(i, j);
                }
                b(row_v) = v_left.c;
                continue;
            }

            // --- Interior Points (Differential Equations) ---

            // Equation 1 (Applied to row_u)
            // L_u(u) + L_v(v) = f1
            for (int j = 0; j < num_pts; ++j) {
                // u-terms block (Top-Left)
                double term_u = eq1.coeff_u_d2(xi) * D2(i, j) +
                                eq1.coeff_u_d1(xi) * D(i, j) +
                                eq1.coeff_u(xi) * (i == j ? 1.0 : 0.0);
                A(row_u, j) += term_u;

                // v-terms block (Top-Right coupling)
                double term_v = eq1.coeff_v_d2(xi) * D2(i, j) +
                                eq1.coeff_v_d1(xi) * D(i, j) +
                                eq1.coeff_v(xi) * (i == j ? 1.0 : 0.0);
                A(row_u, j + num_pts) += term_v;
            }
            b(row_u) = eq1.rhs(xi);

            // Equation 2 (Applied to row_v)
            // L_u(u) + L_v(v) = f2
            for (int j = 0; j < num_pts; ++j) {
                // u-terms block (Bottom-Left coupling)
                double term_u = eq2.coeff_u_d2(xi) * D2(i, j) +
                                eq2.coeff_u_d1(xi) * D(i, j) +
                                eq2.coeff_u(xi) * (i == j ? 1.0 : 0.0);
                A(row_v, j) += term_u;

                // v-terms block (Bottom-Right)
                double term_v = eq2.coeff_v_d2(xi) * D2(i, j) +
                                eq2.coeff_v_d1(xi) * D(i, j) +
                                eq2.coeff_v(xi) * (i == j ? 1.0 : 0.0);
                A(row_v, j + num_pts) += term_v;
            }
            b(row_v) = eq2.rhs(xi);
        }

        // Solve Linear System
        ublas::permutation_matrix<std::size_t> pm(sys_size);
        if (ublas::lu_factorize(A, pm) != 0) {
            throw std::runtime_error("Singular Matrix - Solution Failed");
        }
        ublas::lu_substitute(A, pm, b);

        // Extract u and v
        vector_t u(num_pts), v(num_pts);
        for(int i=0; i<num_pts; ++i) {
            u(i) = b(i);
            v(i) = b(i + num_pts);
        }
        return {u, v};
    }

    // Accessor for grid
    const vector_t& getGrid() const { return x; }
};

// ==========================================
// 3. User Implementation Example
// ==========================================

int run() {
    int N = 30;
    CoupledSolver solver(N);

    // Eq 1: u'' - 4u + v = sin(x)
    EquationCoeffs eq1;
    eq1.coeff_u_d2 = CoupledSolver::Zero;
    eq1.coeff_u_d1 = [](double x) {return 1.0;};
    eq1.coeff_u    = CoupledSolver::Zero;
    eq1.coeff_v_d2 = CoupledSolver::Zero;
    eq1.coeff_v_d1 = CoupledSolver::Zero;
    eq1.coeff_v    = [](double x){ return -1.0; };  // + 1 * v
    eq1.rhs        = CoupledSolver::Zero;

    // Eq 2: v'' - 4v + u = cos(x)
    EquationCoeffs eq2;
    eq2.coeff_u_d2 = CoupledSolver::Zero;
    eq2.coeff_u_d1 = CoupledSolver::Zero;
    eq2.coeff_u    = [](double x){ return 1.0; };  // + 1 * u
    eq2.coeff_v_d2 = [](double x){ return 1.0; };  // 1 * v''
    eq2.coeff_v_d1 = CoupledSolver::Zero;
    eq2.coeff_v    = [](double x){ return -4.0; }; // -4 * v
    eq2.rhs        = [](double x){ return std::cos(x); };

    // 3. Define Boundary Conditions (Dirichlet)
    // Structure: {a, b, c} -> a*y + b*y' = c
    // u(-1)=0, u(1)=0 -> {1, 0, 0}
    BoundaryCondition bc_dirichlet_zero = {1.0, 0.0, 0.0};

    // Example: Neumann on right side for v? v'(1) = 0 -> {0, 1, 0}
    // Let's stick to Dirichlet for now as per prompt example

    try {
        auto result = solver.solve(
            eq1, eq2,
            bc_dirichlet_zero, bc_dirichlet_zero, // u left/right
            bc_dirichlet_zero, bc_dirichlet_zero  // v left/right
        );

        // Output
        vector_t x = solver.getGrid();
        vector_t u = result.first;
        vector_t v = result.second;

        std::cout << std::setw(10) << "x" << std::setw(15) << "u" << std::setw(15) << "v" << "\n";
        std::cout << std::string(40, '-') << "\n";
        for(int i=0; i<=N; ++i) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(4) << x(i)
                      << std::setw(15) << u(i)
                      << std::setw(15) << v(i) << "\n";
        }

    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
