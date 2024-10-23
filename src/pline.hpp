#pragma once

#include "core.hpp" 


namespace geo {
    class Pline {
    public:
        // Declaration of the cut method
        static std::vector<Point> cut(std::vector<Point>& points, const Plane& plane);
    };
}