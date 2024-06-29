#pragma once
#include "core.hpp"
#include "arc.hpp"
#include "offset_2d.hpp"

namespace geo
{
    /**
     * @class shell_2d
     * @brief Represents a 2D shell element with reactions, eccentricities, weight, and centroid.
     *
     * This class models a 2D shell element, encapsulating its physical and geometric properties,
     * including the reactions at the shell's supports, the eccentricities of these reactions,
     * the weight of the shell, and the calculated centroid based on its geometry.
     * 
     * Example:
     * @code
     * // Create a 2D shell element
     * geo::shell_2d shell({{-4, 0, 0}, {0, 0.4, 0}, {4, 0, 0}}, 0.12, 480, 1000, 0.04, 0.04);
     */
    class Shell{

        public:

            // UNITS ARE NEWTONS AND METERS

            Vector Vl; ///< Left vertical reaction.
            Vector Hl; ///< Left horizontal reaction.
            Line el; ///< Left eccentricity line.

            Vector Vr; ///< Right vertical reaction.
            Vector Hr; ///< Right horizontal reaction.
            Line er; ///< Right eccentricity line.

            Vector W; ///< Weight of the shell.
            Point c; ///< Calculated centroid of the shell.


            Shell(std::array<Point, 3> arc_points, double thickness, double density, double horizontal_reaction, double eccentricity0, double eccentricity1, int divisions = 10);



    };


    class Block{

        public:
            double weight;

    };


} // namespace geo