#pragma once

#include "Point.hpp"
#include "Vector.hpp"

namespace geo
{
    /**
     * @class Line
     * @brief A class representing a line in 3D space.
     */
    class Line
    {
    public:
        /**
         * @brief Default constructor. Initializes the line with two points at the origin.
         */
        Line();

        /**
         * @brief Constructor that initializes the line with the given points.
         * @param p1 The first point.
         * @param p2 The second point.
         */
        Line(const Point& p1, const Point& p2);

        /**
         * @brief Subscript operator for non-const access.
         * @param index The index of the point to access (0 for the first point, 1 for the second point).
         * @return A reference to the point at the given index.
         */
        Point& operator[](int index);

        /**
         * @brief Subscript operator for const access.
         * @param index The index of the point to access (0 for the first point, 1 for the second point).
         * @return A const reference to the point at the given index.
         */
        const Point& operator[](int index) const;

        /**
         * @brief Calculates the length of the line.
         * @return The length of the line.
         */
        double length() const;

        /**
         * @brief Calculates the squared length of the line.
         * @return The squared length of the line.
         */
        double squared_length() const;

        /**
         * @brief Converts the line to a vector.
         * @return A vector that points from the first point to the second point of the line.
         */
        Vector to_vector() const;

        /**
         * @brief Calculates the point at a given parameter t.
         * @param t The parameter value.
         * @return The point at the given parameter value.
         */
        Point point_at(const double& t) const;

        /**
         * @brief Scales the line up by a factor of global SCALE.
         */
        void scale_up();

        /**
         * @brief Scales the line down by a factor of global SCALE.
         */
        void scale_down();

    private:
        /**
         * @brief The points of the line.
         */
        Point points[2];
    };
}