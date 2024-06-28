#include "globals.hpp"
#include <cmath>

#pragma once

namespace geo
{
    /**
     * @class Point
     * @brief A class representing a point in 3D space.
     */
    class Point
    {
    public:
        /**
         * @brief Default constructor. Initializes the point at the origin.
         */
        Point();

        /**
         * @brief Constructor that initializes the point at the given _xyz.
         * @param x The x-coordinate.
         * @param y The y-coordinate.
         * @param z The z-coordinate.
         */
        Point(double x, double y, double z=0);

        /**
         * @brief Subscript operator for non-const access.
         * @param index The index of the coordinate to access (0 for x, 1 for y, 2 for z).
         * @return A reference to the coordinate at the given index.
         */
        double& operator[](int index);

        /**
         * @brief Subscript operator for const access.
         * @param index The index of the coordinate to access (0 for x, 1 for y, 2 for z).
         * @return A const reference to the coordinate at the given index.
         */
        const double& operator[](int index) const;

        /**
         * @brief Scale the Point by a factor (multiplication).
         * @param factor The factor to scale the Point by.
         * @return A reference to the Point after scaling.
         */
        Point& operator*=(double factor);

        /**
         * @brief Scale the Point by a factor (division).
         * @param factor The factor to scale the Point by.
         * @return A reference to the Point after scaling.
         */
        Point& operator/=(double factor);

        /**
         * @brief Add another Point to this Point.
         * @param other The Point to add.
         * @return A reference to the Point after addition.
         */
        Point& operator+=(const Point &other);

        /**
         * @brief Subtract another Point from this Point.
         * @param other The Point to subtract.
         * @return A reference to the Point after subtraction.
         */
        Point& operator-=(const Point &other);

        /**
         * @brief Multiply the Point by a factor.
         * @param factor The factor to multiply the Point by.
         * @return A new Point that is the result of the multiplication.
         */
        Point operator*(double factor) const;

        /**
         * @brief Divide the Point by a factor.
         * @param factor The factor to divide the Point by.
         * @return A new Point that is the result of the division.
         */
        Point operator/(double factor) const;

        /**
         * @brief Add another Point to this Point.
         * @param other The Point to add.
         * @return A new Point that is the result of the addition.
         */
        Point operator+(const Point &other) const;

        /**
         * @brief Subtract another Point from this Point.
         * @param other The Point to subtract.
         * @return A new Point that is the result of the subtraction.
         */
        Point operator-(const Point &other) const;

        /**
         * @brief Check if this Point is equal to another Point.
         * @param other The Point to compare with.
         * @return True if the Points are equal, false otherwise.
         */
        bool operator==(const Point &other) const;

        /**
         * @brief Check if this Point is not equal to another Point.
         * @param other The Point to compare with.
         * @return True if the Points are not equal, false otherwise.
         */
        bool operator!=(const Point &other) const;

        /**
         * @brief Negate the Point.
         * @return A new Point that is the result of the negation.
         */
        Point operator-() const;

        /**
         * @brief Retrieves the x-coordinate.
         * @return The x-coordinate.
         */
        double x() const;

        /**
         * @brief Retrieves the y-coordinate.
         * @return The y-coordinate.
         */
        double y() const;

        /**
         * @brief Retrieves the z-coordinate.
         * @return The z-coordinate.
         */
        double z() const;

        /**
         * @brief Scales the point by a given factor.
         * @param factor The factor to scale by.
         */
        void scale(double factor);

        /**
         * @brief Check if the three points of triangle are in a counterclockwise order.
         * @param a The first point.
         * @param b The second point.
         * @param c The third point.
         * @return A positive value if the points are in counterclockwise order, a negative value if they are in clockwise order, and 0 if they are collinear.
         */
        static float Point::ccw(Point &a, Point &b, Point &c);

        /**
         * @brief Computes the midpoint between two points.
         * @param p The other point.
         * @return The midpoint between the two points.
         */
        Point mid_point(Point &p);

        /**
         * @brief Computes the distance between two points.
         * @param p The other point.
         * @return The distance between the two points.
         */
        double distance(Point &p);

        /**
         * @brief Scales the point up by a factor of global SCALE.
         */
        void scale_up();

        /**
         * @brief Scales the point down by a factor of global SCALE.
         */
        void scale_down();

    private:
        /**
         * @brief The coordinates of the point.
         */
        double _xyz[3];
    };
}  // namespace geo