#include <iostream>
#include <iomanip>
#include <array>
#include "../../src/gaussian_elimination.hpp"

int main() {
    // Solving a 3x3 system of linear equations:
    // 2x + y - z = 8
    // -3x - y + 2z = -11
    // -2x + y + 2z = -3

    // Create the augmented matrix [A|b] where A is the coefficient matrix
    // and b is the right-hand side vector
    std::array<std::array<double, 4>, 3> matrix = {{
        {2, 1, -1, 8},    // First equation: 2x + y - z = 8
        {-3, -1, 2, -11},  // Second equation: -3x - y + 2z = -11
        {-2, 1, 2, -3}     // Third equation: -2x + y + 2z = -3
    }};

    // Print the system of equations
    std::cout << "Solving the system of equations:" << std::endl;
    std::cout << "2x + y - z = 8" << std::endl;
    std::cout << "-3x - y + 2z = -11" << std::endl;
    std::cout << "-2x + y + 2z = -3" << std::endl;
    std::cout << std::endl;

    // Solve the system using Gaussian elimination
    auto solution = gaussian_elimination<3>(matrix);

    // Check if a solution exists
    bool has_solution = true;
    for (const auto& x : solution) {
        if (std::isnan(x)) {
            has_solution = false;
            break;
        }
    }

    // Print the solution
    if (has_solution) {
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
    } else {
        std::cout << "No solution exists or system is inconsistent." << std::endl;
    }

    return 0;
}