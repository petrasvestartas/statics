#pragma once

#include "core.hpp"

namespace geo {
class Pline {
   public:
    // Declaration of the cut method
    // https://github.com/petrasvestartas/compas_timbervaultedfloor/blob/main/src/compas_timbervaultedfloor/floor_geometry/mesh_conic_projection.py
    static bool cut(std::vector<Point>& points, const Plane& plane, std::vector<Point>& result);
};
}  // namespace geo