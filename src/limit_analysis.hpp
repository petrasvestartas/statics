#pragma once

#include "Point.hpp"
#include "Vector.hpp"
#include "Line.hpp"

namespace geo
{


    /**
     * @class shell_2d
     * @brief Represents a 2D shell element with reactions, eccentricities, weight, and centroid.
     *
     * This class models a 2D shell element, encapsulating its physical and geometric properties,
     * including the reactions at the shell's supports, the eccentricities of these reactions,
     * the weight of the shell, and the calculated centroid based on its geometry.
     */
    class shell_2d{

        public:
            std::array<Vector, 2> reactions; ///< Reactions at the shell's supports.
            std::array<Line, 2> eccentricities; ///< Eccentricities of the reactions.
            double weight; ///< Weight of the shell.
            geo::Point centroid; ///< Calculated centroid of the shell.


            shell_2d(std::array<Point, 3> arc_points, double thickness, double density, double horizontal_reaction, double eccentricity0, double eccentricity1){
                
            }



    };


    class block_2d{

        public:
            double weight;

    };


} // namespace geo