            // CGAL::Aff_transformation_3<IK> plane_to_xy(const IK::Point_3 &origin, const IK::Plane_3 &plane)
            // {
            //     auto x0 = plane.base1();
            //     auto y0 = plane.base2();
            //     auto z0 = plane.orthogonal_vector();
            //     unitize(x0);
            //     unitize(y0);
            //     unitize(z0);

            //     // Move to origin -> T0 translates point P0 to (0,0,0)
            //     CGAL::Aff_transformation_3<IK> t0(CGAL::TRANSLATION, IK::Vector_3(-origin.x(), -origin.y(), -origin.z()));

            //     // Rotate ->
            //     CGAL::Aff_transformation_3<IK> f0(
            //         x0.x(), x0.y(), x0.z(),
            //         y0.x(), y0.y(), y0.z(),
            //         z0.x(), z0.y(), z0.z());

            //     return f0 * t0;
            // }