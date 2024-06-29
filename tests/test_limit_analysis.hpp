#pragma once
#include "./src/core.hpp" // Include the header of the class you're testing
#include "limit_analysis.hpp"
#include <iomanip>
#include <cmath>
#include "test_methods.hpp"

void test_limit_analysis() {

    geo::Shell shell(std::array<geo::Point, 3>{geo::Point{-4, 0, 0}, geo::Point{0, 0.4, 0}, geo::Point{4, 0, 0}}, 0.12, 480, 1000, 0.04, 0.04, 10);

    // my_assert(std::abs(arc.radius - 20.2) < 0.0001);


}

int test_limit_analysis_main() {
    test_limit_analysis();
    return 0;
}