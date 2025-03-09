#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include "../../src/gaussian_elimination.hpp"
#include "../../src/matrix.hpp"

int main() {
    // Solving a 3x3 system of linear equations:
    // 2x + y - z = 8
    // -3x - y + 2z = -11
    // -2x + y + 2z = -3

    // Create the coefficient matrix A and right-hand side vector b
    geo::Matrix A{{2, 1, -1},    // First equation: 2x + y - z = 8
                  {-3, -1, 2},   // Second equation: -3x - y + 2z = -11
                  {-2, 1, 2}};   // Third equation: -2x + y + 2z = -3
    std::vector<double> b{8, -11, -3};

    // Print the system of equations
    std::cout << "Solving the system of equations:" << std::endl;
    std::cout << "2x + y - z = 8" << std::endl;
    std::cout << "-3x - y + 2z = -11" << std::endl;
    std::cout << "-2x + y + 2z = -3" << std::endl;
    std::cout << std::endl;

    // Solve the system using Gaussian elimination
    std::vector<double> solution;
    try {
        solution = geo::gauss_partial(A, b);

        // Print the solution
        std::cout << "Solution:" << std::endl;
        std::cout << "x = " << std::fixed << std::setprecision(6) << solution[0] << std::endl;
        std::cout << "y = " << solution[1] << std::endl;
        std::cout << "z = " << solution[2] << std::endl;

        // Verify the solution by substituting back into the equations
        std::cout << "\nVerification:" << std::endl;
        double eq1 = 2*solution[0] + solution[1] - solution[2];
        double eq2 = -3*solution[0] - solution[1] + 2*solution[2];
        double eq3 = -2*solution[0] + solution[1] + 2*solution[2];

        std::cout << "Equation 1: " << eq1 << " = 8" << std::endl;
        std::cout << "Equation 2: " << eq2 << " = -11" << std::endl;
        std::cout << "Equation 3: " << eq3 << " = -3" << std::endl;
    } catch (const std::runtime_error& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    return 0;
}