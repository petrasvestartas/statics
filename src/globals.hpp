#pragma once

namespace geo
{


   struct GLOBALS
   {
        public:
            static bool IS_FINITE(double x) { return 0x7FF0 != (*((unsigned short *)(&x) + 3) & 0x7FF0); }
            static constexpr double DOUBLE_MIN = 2.22507385850720200e-308;
            static constexpr double DOUBLE_MAX = 1.7976931348623158e+308;
            static constexpr double EPSILON = 2.2204460492503131e-16;
            static constexpr double SQRT_EPSILON = 1.490116119385000000e-8;
            
            static constexpr double M_PI = 3.14159265358979323846;
            static constexpr double TO_RAD =   0.01745329251994329576;
            static double SCALE;
            static double ANGLE;

            static constexpr double GIGA = 1e9;
            static constexpr double MEGA = 1e6;
            static constexpr double KILO = 1e3;
            static constexpr double ZERO_TOLERANCE = 2.3283064365386962890625e-10;
            static constexpr double MILLI = 1e-3;
            static constexpr double MICRO = 1e-6;
            static constexpr double NANO = 1e-9;

            static constexpr double FORCE = 1;
            static constexpr double MASS = 1;
            static constexpr double LENGTH = 1;
   };




} // namespace geo