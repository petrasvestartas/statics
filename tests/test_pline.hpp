#pragma once
#include "core.hpp" // Include the header of the class you're testing
#include "pline.hpp"
#include <iomanip>
#include <cmath>
#include "test_methods.hpp"

void test_cut() {
    
    std::vector<geo::Point> points = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {0, 0, 0}
    };

    geo::Plane plane(geo::Point(0, 0, 0), geo::Vector(0, 0, 1));
    std::vector<geo::Point> result = geo::Pline::cut(points, plane); 

}


int test_pline_main() {
    test_cut();
    return 0;
}