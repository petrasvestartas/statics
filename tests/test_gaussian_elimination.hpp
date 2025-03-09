#pragma once
#include <cmath>
#include <array>
#include <limits>
#include "test_methods.hpp"
#include "gaussian_elimination.hpp"


void test_gaussian_elimination() {
    std::cout << "Testing Gaussian Elimination..." << std::endl;

    // Test 1: Simple 2x2 system
    // x + y = 3
    // 2x + y = 4
    // Solution: x = 1, y = 2
    std::array<std::array<double, 3>, 2> matrix1 = {{
        {1, 1, 3},
        {2, 1, 4}
    }};
    auto result1 = gaussian_elimination<2>(matrix1);
    my_assert(std::abs(result1[0] - 1.0) < 1e-10 && std::abs(result1[1] - 2.0) < 1e-10);

    // Test 2: 3x3 system
    // 2x + y - z = 8
    // -3x - y + 2z = -11
    // -2x + y + 2z = -3
    // Solution: x = 2, y = 3, z = -1
    std::array<std::array<double, 4>, 3> matrix2 = {{
        {2, 1, -1, 8},
        {-3, -1, 2, -11},
        {-2, 1, 2, -3}
    }};
    auto result2 = gaussian_elimination<3>(matrix2);
    my_assert(std::abs(result2[0] - 2.0) < 1e-10 && 
              std::abs(result2[1] - 3.0) < 1e-10 && 
              std::abs(result2[2] - (-1.0)) < 1e-10);

    // Test 3: Singular matrix (no unique solution)
    // x + y = 1
    // x + y = 2
    std::array<std::array<double, 3>, 2> matrix3 = {{
        {1, 1, 1},
        {1, 1, 2}
    }};
    auto result3 = gaussian_elimination<2>(matrix3);
    my_assert(std::isnan(result3[0]) && std::isnan(result3[1]));

    // Test 4: Zero pivot test with row swapping
    // 0x + y = 1
    // x + y = 2
    std::array<std::array<double, 3>, 2> matrix4 = {{
        {0, 1, 1},
        {1, 1, 2}
    }};
    auto result4 = gaussian_elimination<2>(matrix4);
    my_assert(std::abs(result4[0] - 1.0) < 1e-10 && std::abs(result4[1] - 1.0) < 1e-10);

    std::cout << "All Gaussian Elimination tests passed!" << std::endl;
}

int test_gaussian_elimination_main() {
    test_gaussian_elimination();
    return 0;
}