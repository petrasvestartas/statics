#pragma once
#include "core.hpp" // Include the header of the class you're testing
#include "limit_analysis.hpp"
#include <iomanip>
#include <cmath>
#include "test_methods.hpp"

void test_limit_analysis() {

    // horizontal force - 16.4 18.4, 26.2 - 29.3 kN
    geo::Shell shell(std::array<geo::Point, 3>{geo::Point{-4, 0, 0}, geo::Point{0, 0.4, 0}, geo::Point{4, 0, 0}}, 0.12, 480, 9145, -0.03, -0.03, 10);

    // my_assert(std::abs(arc.radius - 20.2) < 0.0001);

}

int test_limit_analysis_main() {
    test_limit_analysis();
    return 0;
}