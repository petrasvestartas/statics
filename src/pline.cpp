#include "pline.hpp"

#include <algorithm>

namespace geo {
    std::vector<Point> Pline::cut(std::vector<Point>& points, const Plane& plane) {
        std::vector<Point> result;
        Point point = plane.get_point();
        Vector normal = plane.get_normal();
        log(point.to_string());
        auto vector_name = normal.to_string();
        log(vector_name);



        // # Orient the polygon from the cutting frame to the XY frame
        // frame = Frame.from_plane(plane)
        // xform = Transformation.from_frame_to_frame(frame, Frame.worldXY())
        // xformI = Transformation.from_frame_to_frame(Frame.worldXY(), frame)
        
        
        // polygon_xformed = polygon.transformed(xform)
        
        // points = polygon_xformed.points
        
        // flag = True
        // below = []
        
        // # Check if the points are below the plane
        // counter = 0;
        // for i in range(len(points)):
        //     if (points[i][2] < 0):
                
        //         below.append(True)
        //         flag = False
        //         counter = counter + 1
        //     else:
        //         below.append(False)
                
        // if (flag):
        //     return [polygon]
        
        // print(flag)
                    
        // # For faces that coincide with plane
        // polygons_culled = []

        // if (counter != 0 and counter != len(points)):
            
        //     # Get split parameters
        //     polygon_points = []
        //     polygon_points_bool = []
            

        //     counter = 0
        //     for i in range(len(points)):
        //         polygon_points.append(points[i])
        //         polygon_points_bool.append(False)
        //         # points_lists[-1].append(points[i])
        //         line = Line(points[i], points[(i + 1)%len(points)])
        //         result = intersection_segment_plane(line, plane.worldXY())
        //         print("_",result)
                
        //         if result:
        //             intersection_point = Point(*result)
        //             polygon_points.append(intersection_point)
        //             polygon_points_bool.append(True)
                    
        //             if counter == 0:
        //                 counter = len(polygon_points)-1
        //             # points_lists[-1].append(intersection_point)
        //             # points_lists.append([intersection_point])
            
        //     # polygon_points.append(points[-1])
        //     # polygon_points_bool.append(False)
            
        //     print("counter", counter)
            
        //     def shift_right(lst, n):
        //         n = n % len(lst)  # To handle cases where n > len(lst)
        //         return lst[-n:] + lst[:-n]
                    
        //     def shift_left(lst, n):
        //         n = n % len(lst)  # To handle cases where n > len(lst)
        //         return lst[n:] + lst[:n]
            
        //     print(polygon_points_bool)
        //     polygon_points = shift_left(polygon_points, counter)
        //     polygon_points_bool = shift_left(polygon_points_bool, counter)
            
        //     print(polygon_points_bool)
            
        //     # split points into sub lists
        //     points_lists = []
        //     for i in range(len(polygon_points)):
        //         print(polygon_points[i], polygon_points_bool[i])
        //         if polygon_points_bool[i]:
        //             points_lists.append([])
                    
        //         if len(points_lists) > 1 and polygon_points_bool[i]:
        //             points_lists[-2].append(polygon_points[i])

        //         points_lists[-1].append(polygon_points[i])
                
        //     points_lists[-1].append(polygon_points[0])
        //     print(polygon_points[0])
                
            

        //     for list in points_lists:
        //         print(len(list))
        //         for p in list:
        //             print(p)
        
                
        //     from compas import json_dump
        //     json_dump(points_lists, 'C:/brg/2_code/compas_timbervaultedfloor/data/json_dump/debug.json')
                
        //     for i in range(len(points_lists)):
        //         cut_polyline_part = Polyline(points_lists[i])
        //         if (cut_polyline_part.length<tolerance):
        //             continue
        //         print(cut_polyline_part.point_at(0.5, True)[2])
        //         if (cut_polyline_part.point_at(0.5, True)[2] > tolerance):
        //             cut_polyline_part.transform(xformI)

        //             polygons_culled.append(Polygon(cut_polyline_part.points))

        // print("Polygons culled: ", len(polygons_culled), " Points: ", len(polygons_culled[0]))
        // return polygons_culled

        return result;
    }
}