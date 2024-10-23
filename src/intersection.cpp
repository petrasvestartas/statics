#include "intersection.hpp"

namespace geo {
    // Constructor
    Intersection::Intersection() {
        // Constructor implementation
    }

    bool Intersection::line_plane(const Line &line, const Plane &plane, Point &output, bool is_finite) {

        
        bool result = false;
        double a, b, d, fd, t;

        auto pt0 = line[0];
        auto pt1 = line[1];
        a = plane.value_at(pt0);
        b = plane.value_at(pt1);
        d = a - b;
        if (d == 0.0)
        {
            if (fabs(a) < fabs(b))
                t = 0.0;
            else if (fabs(b) < fabs(a))
                t = 1.0;
            else
                t = 0.5;
        }
        else
        {
            d = 1.0 / d;
            fd = fabs(d);
            if (fd > 1.0 && (fabs(a) >= geo::GLOBALS::DOUBLE_MAX / fd || fabs(b) >= geo::GLOBALS::DOUBLE_MAX / fd))
            {
                // double overflow - line may be (nearly) parallel to plane
                t = 0.5;
            }
            else
            {
                t = a / (a - b); // = a*d;  a/(a-b) has more precision than a*d
                result = true;
            }
        }

        const double s = 1.0 - t;

        output = Point(
            (line[0].x() == line[1].x()) ? line[0].x() : s * line[0].x() + t * line[1].x(),
            (line[0].y() == line[1].y()) ? line[0].y() : s * line[0].y() + t * line[1].y(),
            (line[0].z() == line[1].z()) ? line[0].z() : s * line[0].z() + t * line[1].z());

        if (is_finite && (t < 0 || t > 1))
            return false;

        return result;
    }

}