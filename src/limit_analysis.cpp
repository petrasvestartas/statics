#include "limit_analysis.hpp"

namespace geo
{

    Shell::Shell(std::array<Point, 3> arc_points, double thickness, double density, double horizontal_reaction, double eccentricity0, double eccentricity1, int divisions){

        // Three point arc
        geo::Point p0(-4, 0, 0);
        geo::Point p1(0, 0.4, 0);
        geo::Point p2(4, 0, 0);
        geo::Arc arc(arc_points[0], arc_points[1], arc_points[2]);

        // Divide arc into points
        std::vector<geo::Point> points = arc.divide_arc_into_points(divisions);

        // Offset polyline
        std::vector<geo::Point> offset0 = geo::offset_polyline(points, -thickness*0.5);
        std::vector<geo::Point> offset1 = geo::offset_polyline(points, thickness*0.5);

        // Create rectangular polygons
        std::vector<std::vector<geo::Point>> polygons;
        for (size_t i = 0; i < points.size() - 1; ++i) {
            std::vector<geo::Point> polygon = {offset0[i], offset0[i + 1], offset1[i + 1], offset1[i]};
            polygons.push_back(polygon);
        }

        // Calculate weight as negative vertical vector
        W = Vector(0,0,0);
        for (auto& polygon : polygons) {
            W[1] -= Point::area(polygon)*density;
        }

        // Set one horizontal force and eccentricities
        Hr = Vector(horizontal_reaction, 0, 0);
        Vector el_v = offset0.front()-offset1.front().to_vector();
        Vector er_v = offset0.back()-offset1.back().to_vector();
        el_v.unitize();
        er_v.unitize();
        Vector temp_el_v = el_v * eccentricity0; // Ensure this operation is valid and returns a Vector
        Point newPoint = offset1.front() + temp_el_v; // Use operator+ here
        el = Line(offset1.front(), newPoint);
        el = Line(offset1.front(), offset1.front() + el_v * eccentricity0);
        er = Line(offset1.back(), offset1.back() + er_v * eccentricity1);


        // Move geometry to set the right top eccentricity point as the origin
        Vector translation_vector = Vector(-static_cast<double>(er[1][0]), -static_cast<double>(er[1][1]), -static_cast<double>(er[1][2]));
        for (auto& polygon : polygons) {
            for (auto& point : polygon) {
                point.translate(translation_vector);
            }
        }

        el.translate(translation_vector);
        er.translate(translation_vector);

        // Print location to verify the points of the arc
        for (auto& polygon : polygons) {
            for (auto& point : polygon) {
                std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl;
            }
        }

        // Apply Horizontal, Vertical and Moment equilibrium to find the right reaction force
        // Vl + Vr = W
        // Hr = Hl
        // W[0] * centroid[0] = Vl * eccentricities[0][0][1] + Vr * eccentricities[1][1][0]
        
    }



} // namespace geo