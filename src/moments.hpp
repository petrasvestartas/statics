#pragma once
#include "core.hpp"
#include "point.hpp"
#include "vector.hpp"

namespace geo{

    /**
     * @brief Get the moment sign object.
     * CW = 1, CCW = -1, COLLINEAR = 0
     * 
     * /+\    |    /-\   |      |      |
     *  |     |     |    |      |      |
     *  |     |     |   \-/     |     \+/
     * 
     * —————— + ————>    <————— - —————
     * ——————————————    ——————————————
     * —————— - ————>    <————— + —————
     * 
     * @param p Point
     * @param f Vector
     * @param cw_is_positive bool, change to false if you want CCW to be positive
     * @return std::array<int, 2>
     */
    std::array<int, 2> get_moment_sign(const Point &p, const Vector &f, bool cw_is_positive = true);

    /**
     * @brief Get the moment object.
     * @param p Point
     * @param f Vector
     * @return std::array<double, 2>
     */
    std::array<double, 2> get_moment(const Point &p, const Vector &f);

    /**
     * @brief Get the moment with eccentricity object.
     *
     *   fx  | fy
     *   --> V
     *       *    e
     *        \
     *         *    p
     * 
     * @param p Point
     * @param f Vector force applied at the point
     * @param e Point eccentricity point, this point is used to adjust the force by the distance between p and e
     * @return std::array<double, 2>
     */
    std::array<double, 2> get_moment_with_eccentricity(const Point &p, const Vector &f, const Point &e);
}