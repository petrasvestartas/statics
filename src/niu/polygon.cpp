            // double get_closest_distance(const IK::Point_3 &point, CGAL_Polyline &polyline, size_t &edge_id)
            // {
            //     edge_id = 0;
            //     IK::Segment_3 segment(polyline[0], polyline[1]);
            //     double closest_distance = DBL_MAX;

            //     for (size_t i = 0; i < polyline.size() - 1; i++)
            //     {
            //         IK::Segment_3 segment_(polyline[i], polyline[i + 1]);

            //         double t;
            //         closest_point_to(point, segment_, t);

            //         double closest_distance_temp = std::abs(CGAL::squared_distance(point, point_at(segment_, t)));
            //         if (closest_distance_temp < closest_distance_temp)
            //         {
            //             closest_distance = closest_distance_temp;
            //             edge_id = i;
            //         }

            //         if (closest_distance < wood::GLOBALS::DISTANCE_SQUARED)
            //             break;
            //     }

            //     return closest_distance;
            // }