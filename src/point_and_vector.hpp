#pragma once
#include "./core.hpp" 

namespace geo{
    
    /**
     * @brief Adds a vector to a point.
     * 
     * @param p The point to which the vector will be added.
     * @param v The vector to be added to the point.
     * @return Point The result of adding the vector to the point.
     */
    Point operator+(const Point& p, const Vector& v);

    /**
     * @brief Subtracts a vector from a point.
     * 
     * @param p The point from which the vector will be subtracted.
     * @param v The vector to be subtracted from the point.
     * @return Point The result of subtracting the vector from the point.
     */
    Point operator-(const Point& p, const Vector& v);


    /**
    * @brief Given a 2d vector, decompose it into x and y components. Then determine the signs of each vector.
    * Source Pierre Varignon's theorem. Source: Engineering Mechanics Statics, by R. C. Hibbeler, 12th edition, page 13.
    * 
    * ---------> -
    * ---------- x
    * ---------> +
    * 
    * <--------- +
    * ---------- x
    * <--------- -
    * 
    *  | + | - |
    *  |   |   |
    *  |   |   |
    * \|/  |  \|/
    * 
    * /|\  |  /|\
    *  |   |   |
    *  |   |   |
    *  | - | + |
    * 
    * @param [in] p 2d point
    * @param [in] v 2d vector
    * @return : x and y axes sign, positive sign, counter-clockwise, -1: negative sign, clockwise, 0: vector is algined to axis.
    */
    std::array<int, 2> moment_component_signs_varignon(Point& p, Vector &v);

    /**
     * @brief Compute the moment of a force about a point using Varignon's theorem. Origin is considere 0,0.
     * 
     * @param point The point of application of the force.
     * @param force The force vector.
     * @return double The moment of the force about the point.
     */
    double moment_varignon(Point& point, Vector& force);

    /**
     * Compute the sum of moments.
     * As a convension, the moment is positive if it is counterclockwise and negative if it is clockwise.
     * If the sum is positive the object will rotate counterclockwise, if it is negative it will rotate clockwise.
     * 
     * @param [in] origins point of application
     * @param [in] forces vector of forces, the length of the vector is the magnitude of the force
     * @return the sum of moments
     */
    double moments_varignon_sum(std::vector<Point>& origins, std::vector<Vector>& forces);

}