#include "Plane.hpp"
#include <cmath>

namespace geo
{
    Plane::Plane(const Point& point, const Vector& normal)
        : point(point), normal(normal) {}

    Plane::Plane(const double equation[4]) {
        normal = Vector(equation[0], equation[1], equation[2]);
        // Choose a point on the plane. If the normal is not zero, we can choose a point on the x, y or z axis.
        if (normal[0] != 0) {
            point = Point(-equation[3] / equation[0], 0, 0);
        } else if (normal[1] != 0) {
            point = Point(0, -equation[3] / equation[1], 0);
        } else if (normal[2] != 0) {
            point = Point(0, 0, -equation[3] / equation[2]);
        } else {
            // If the normal is zero, the equation represents the entire space, and we can choose any point.
            point = Point(0, 0, 0);
        }
    }

    const Point& Plane::get_point() const {
        return point;
    }

    const Vector& Plane::get_normal() const {
        return normal;
    }

    double Plane::distance_to_point(const Point& point) const {
        double equation[4];
        get_plane_equation(equation);
        double numerator = std::abs(equation[0] * point[0] + equation[1] * point[1] + equation[2] * point[2] + equation[3]);
        double denominator = std::sqrt(equation[0] * equation[0] + equation[1] * equation[1] + equation[2] * equation[2]);
        return numerator / denominator;
    }

    double Plane::squared_distance_to_point(const Point& point) const {
        double equation[4];
        get_plane_equation(equation);
        double numerator = equation[0] * point[0] + equation[1] * point[1] + equation[2] * point[2] + equation[3];
        numerator = numerator * numerator; // Square the numerator
        double denominator = equation[0] * equation[0] + equation[1] * equation[1] + equation[2] * equation[2];
        return numerator / denominator;
    }

    Point Plane::closest_point(const Point& point) const {
        double equation[4];
        get_plane_equation(equation);
        double distance = equation[0] * point[0] + equation[1] * point[1] + equation[2] * point[2] + equation[3];
        return Point(point[0] - equation[0] * distance,
                    point[1] - equation[1] * distance,
                    point[2] - equation[2] * distance);
    }

    void Plane::get_plane_equation(double (&equation)[4]) const {
        equation[0] = normal[0];
        equation[1] = normal[1];
        equation[2] = normal[2];
        equation[3] = -(normal[0] * point[0] + normal[1] * point[1] + normal[2] * point[2]);
    }

}