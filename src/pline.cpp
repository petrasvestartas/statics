#include "pline.hpp"

#include <algorithm>

namespace geo {
    std::vector<Point> Pline::cut(std::vector<Point>& points, const Plane& plane) {
        std::vector<Point> result;
        Point point = plane.get_point();
        Vector normal = plane.get_normal();
        log(point.to_string());
        auto vector_name = normal.to_string();
        log(vector_name);
        // result.reserve(points.size());

        // for (const auto& point : points) {
        //     if (plane.normal().dot(point - plane.origin()) >= 0) {
        //         result.push_back(point);
        //     }
        // }

        return result;
    }
}