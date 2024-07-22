#include "vector.hpp"

namespace geo
{

    Vector::Vector() : _xyz{0.0, 0.0, 0.0} {}

    Vector::Vector(double x, double y, double z){
        _xyz[0] = x;
        _xyz[1] = y;
        _xyz[2] = z;
        _has_length = false;
    }

    Vector Vector::XAxis() {
        return Vector(1, 0, 0);
    }

    Vector Vector::YAxis() {
        return Vector(0, 1, 0);
    }

    Vector Vector::ZAxis() {
        return Vector(0, 0, 1);
    }


    Vector Vector::from_start_and_end(const Vector& start, const Vector& end) {
        return Vector(end._xyz[0] - start._xyz[0], 
                    end._xyz[1] - start._xyz[1], 
                    end._xyz[2] - start._xyz[2]);
    }

    double& Vector::operator[](int index) {
        return _xyz[index];
    }

    const double& Vector::operator[](int index) const {
        return _xyz[index];
    }

    Vector& Vector::operator*=(double factor) {
        for (size_t i = 0; i < 3; ++i)
            _xyz[i] *= factor;
        _has_length = false;
        return *this;
    }

    // Friend non-member operator* overload to handle double * Vector
    Vector operator*(double factor, const Vector& v) {
        return v * factor; // Reuse the member operator* overload
    }

    Vector& Vector::operator/=(double factor) {
        for (size_t i = 0; i < 3; ++i)
            _xyz[i] /= factor;
        _has_length = false;
        return *this;
    }

    Vector& Vector::operator+=(const Vector &other) {
        for (size_t i = 0; i < 3; ++i)
            _xyz[i] += other._xyz[i];
        _has_length = false;
        return *this;
    }

    Vector& Vector::operator-=(const Vector &other) {
        for (size_t i = 0; i < 3; ++i)
            _xyz[i] -= other._xyz[i];
        _has_length = false;
        return *this;
    }

    Vector Vector::operator*(double factor) const {
        Vector result = *this;
        result *= factor;
        result._has_length = false;
        return result;
    }

    Vector Vector::operator/(double factor) const {
        Vector result = *this;
        result /= factor;
        result._has_length = false;
        return result;
    }

    Vector Vector::operator+(const Vector &other) const {
        Vector result = *this;
        result += other;
        result._has_length = false;
        return result;
    }

    Vector Vector::operator-(const Vector &other) const {
        Vector result = *this;
        result -= other;
        result._has_length = false;
        return result;
    }

    bool Vector::operator==(const Vector &other) const {
        // Replace x, y, z with your actual member variables
        return _xyz[0] == other[0] && _xyz[1] == other[1] && _xyz[2] == other[2];
    }

    bool Vector::operator!=(const Vector &other) const {
        return !(*this == other);
    }

    void Vector::reverse() {
        for (size_t i = 0; i < 3; ++i)
            _xyz[i] = -_xyz[i];
    }

    double Vector::length()
    {
        if (!_has_length) {


            double x = abs(_xyz[0]);
            double y = abs(_xyz[1]);
            double z = abs(_xyz[2]);
            if (y >= x && y >= z)
            {
                _length = x;
                x = y;
                y = _length;
            }
            else if (z >= x && z >= y)
            {
                _length = x;
                x = z;
                z = _length;
            }

            // For small denormalized doubles (positive but smaller
            // than DOUBLE_MIN), some compilers/FPUs set 1.0/x to +INF.
            // Without the DOUBLE_MIN test we end up with
            // microscopic vectors that have infinite length!
            if (x > GLOBALS::DOUBLE_MIN)
            {
                y /= x;
                z /= x;
                _length = x * sqrt(1.0 + y * y + z * z);
            }
            else if (x > 0.0 && GLOBALS::IS_FINITE(x))
                _length = x;
            else
                _length = 0.0;
            _has_length = true;
        } 
        
        return _length;
    }

    bool Vector::unitize() {
        bool rc = false;
        double d = length();
        if (d > 0.0)
        {
            _xyz[0]/=d;
            _xyz[1]/=d;
            _xyz[2]/=d;
            rc = true;
            _length = 1;
            _has_length = true;
        }
        return rc;
    }

    Vector Vector::unitized() {
        Vector unitized_vector = *this;
        unitized_vector.unitize();
        return unitized_vector;
    }

    Vector Vector::projection(
        Vector& projection_vector,
        double tolerance,
        double* out_projected_vector_length, 
        Vector* out_perpendicular_projected_vector, 
        double* out_perpendicular_projected_vector_length){

        // unitize vector_to_project_on
        double projection_vector_length = projection_vector.length();
        Vector projection_vector_unit = Vector(
            projection_vector[0]/projection_vector_length, 
            projection_vector[1]/projection_vector_length, 
            projection_vector[2]/projection_vector_length);

        // Before unitizing, check if the length is below a tolerance value
        if (projection_vector_length < tolerance) {
            // Handle the case where the projection_vector_length is too small
            // This could involve setting default values or returning an error
            if (out_projected_vector_length) *out_projected_vector_length = 0;
            if (out_perpendicular_projected_vector) *out_perpendicular_projected_vector = Vector(0, 0, 0);
            if (out_perpendicular_projected_vector_length) *out_perpendicular_projected_vector_length = 0;
            return Vector(0, 0, 0); // Return a zero vector or handle as appropriate
        }

        // get the scale factor of the projection vector by dot product
        double projected_vector_length = this->dot(projection_vector_unit);
        if (out_projected_vector_length) {
            *out_projected_vector_length = projected_vector_length;
        }

        // get the projection vector
        Vector out_projection_vector = projection_vector_unit * projected_vector_length;

        // get the perpendicular vector
        if (out_perpendicular_projected_vector) {
            *out_perpendicular_projected_vector = *this - out_projection_vector;
            if (out_perpendicular_projected_vector_length) {
                *out_perpendicular_projected_vector_length = out_perpendicular_projected_vector->length();
            }
        } else if (out_perpendicular_projected_vector_length) {
            // Calculate length only if needed and perpendicular vector is not requested
            Vector temp_perpendicular_vector = *this - out_projection_vector;
            *out_perpendicular_projected_vector_length = temp_perpendicular_vector.length();
        }

        return out_projection_vector;
    }

    int Vector::is_parallel_to(Vector &v) {
        double ll = length() * v.length();
        int result;

        if (ll > 0.0)
        {
            const double cos_angle = (_xyz[0] * v[0] + _xyz[1] * v[1] + _xyz[2] * v[2]) / ll;

            const double angle_in_radians = GLOBALS::ANGLE * (GLOBALS::PI / 180.0); // convert angle from degrees to radians
            const double cos_tol = cos(angle_in_radians);
            if (cos_angle >= cos_tol)
                result = 1; // "Parallel";
            else if (cos_angle <= -cos_tol)
                result = -1; // "Antiparallel"
            else
                result = 0; // "Not parallel"
        }
        else
        {
            result = 0; // "Not parallel"
        }

        return result;
    }

    double Vector::dot(Vector &other) {
        double result = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            result += _xyz[i] * other._xyz[i];
        }

        return result;
    }

    Vector Vector::cross(Vector &other) {
        double x = _xyz[1] * other._xyz[2] - _xyz[2] * other._xyz[1];
        double y = _xyz[2] * other._xyz[0] - _xyz[0] * other._xyz[2];
        double z = _xyz[0] * other._xyz[1] - _xyz[1] * other._xyz[0];
        Vector result = Vector(x, y, z);
        result.unitize();

        return result;
    }

    double Vector::angle(Vector &other, bool degrees, double tolerance) {

        double dot = this->dot(other);
        double len0 = this->length();
        double len1 = other.length();
        
        // Ensure the denominator is not zero
        double denominator = len0 * len1;
        if (denominator < tolerance) {
            return 0; // Or handle the error as appropriate
        }
        
        // Calculate the cosine of the angle
        double cos_angle = dot / denominator;
        
        // Clamp cos_angle to the range [-1, 1] to avoid NaN errors due to floating point inaccuracies
        cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
        
        // Calculate the angle in radians
        double angle = acos(cos_angle);

        // Determine the sign of the angle using the z-component of the cross product
        Vector crossProduct = this->cross(other);
        if (crossProduct._xyz[2] < 0) {
            angle = -angle;
        }

        // Return the angle in degrees if requested
        if (degrees) {
            angle = angle * (180.0 / GLOBALS::PI);
        } else {
            angle = angle;
        }
        
        // double len0 = this->length();
        // double len1 = other.length();
        // Vector sum =  len1 * *this + len0 * other;
        // Vector diff = len1 * *this - len0 * other;
        // double angle =  2.0 * atan(diff.length() / sum.length());

        // // Determine the sign of the angle using the z-component of the cross product
        // Vector crossProduct = this->cross(other);
        // if (crossProduct._xyz[2] < 0)
        //     angle = -angle;

        return angle;
    }



    Vector Vector::get_leveled_vector(double &vertical_height)
    {

        Vector copy = Vector(_xyz[0], _xyz[1], _xyz[2]);

        if (copy.unitize())
        {
            Vector referenceVector(0, 0, 1);
            double angle = copy.angle(referenceVector, true);
            double inclined_offset_by_vertical_distance = vertical_height / std::cos(angle);
            copy *= inclined_offset_by_vertical_distance;
        }

        return copy;
    }

    double Vector::cosine_law(double& triangle_edge_length_a, double& triangle_edge_length_b, double& angle_in_degrees_between_edges){
        double angle_in_radians_between_edges = angle_in_degrees_between_edges * GLOBALS::PI / 180.0;
        double result = sqrt(pow(triangle_edge_length_a, 2) + pow(triangle_edge_length_b, 2) - 2*triangle_edge_length_a*triangle_edge_length_b*cos(angle_in_radians_between_edges));
        return result;
    }

    double Vector::sine_law_angle(double& triangle_edge_length_a, double& angle_in_degrees_in_front_of_a, double& triangle_edge_length_b){
        double angle_in_radians_in_front_of_a = angle_in_degrees_in_front_of_a * GLOBALS::PI / 180.0;
        return asin( (triangle_edge_length_b * sin(angle_in_radians_in_front_of_a)) / triangle_edge_length_a ) * 180.0 / GLOBALS::PI;
    }

    double Vector::sine_law_length(double& triangle_edge_length_a, double& angle_in_degrees_in_front_of_a, double& angle_in_degrees_in_front_of_b){
        double angle_in_radians_in_front_of_a = angle_in_degrees_in_front_of_a * GLOBALS::PI / 180.0;
        double angle_in_radians_in_front_of_b = angle_in_degrees_in_front_of_b * GLOBALS::PI / 180.0;
        return (triangle_edge_length_a * sin(angle_in_radians_in_front_of_b))/sin(angle_in_radians_in_front_of_a);
    }

    double Vector::angle_between_vector_xy_components_degrees(Vector &vector){
        double angle_between_inclined_vector_and_horizontal_component = atan(vector[1]/vector[0]);
        return angle_between_inclined_vector_and_horizontal_component * 180.0 / GLOBALS::PI;
    }

    Vector Vector::sum_of_vectors(std::vector<Vector> &vectors){
        double x = 0, y = 0, z = 0;
        for (Vector vector : vectors){
            x += vector[0];
            y += vector[1];
            z += vector[2];
        }
        return Vector(x, y, z);
    }

    std::array<double, 3> Vector::coordinate_direction_3angles(Vector &v, bool degrees){
        double x = v[0];
        double y = v[1];
        double z = v[2];
        double r = sqrt(x*x + y*y + z*z);

        // unit vector
        // u = F/|F| = F_x/|F| + F_y/|F| + F_z/|F|
        double x_proportion = x/r;
        double y_proportion = y/r;
        double z_proportion = z/r;

        // angles
        double alpha = acos(x_proportion);
        double beta = acos(y_proportion);
        double gamma = acos(z_proportion);

        if (degrees){
            alpha = alpha * 180.0 / GLOBALS::PI;
            beta = beta * 180.0 / GLOBALS::PI;
            gamma = gamma * 180.0 / GLOBALS::PI;
        }

        return std::array<double, 3>{alpha, beta, gamma};


    }

    std::array<double, 2> Vector::coordinate_direction_2angles(Vector &v, bool degrees){
        double x = v[0];
        double y = v[1];
        double z = v[2];
        double r = sqrt(x*x + y*y + z*z);

        // unit vector
        // u = |v|*sin(phi)*cos(theta) i + |v|*sin(phi)*sin(theta) j + |v|*cos(phi) k
        // from the unit vector, lets get angle phi and theta

        double phi = acos(z/r);
        double theta = atan2(y, x);


        if (degrees){
            phi = phi * 180.0 / GLOBALS::PI;
            theta = theta * 180.0 / GLOBALS::PI;
        }

        return std::array<double, 2>{phi, theta};
        
    }

    void Vector::scale(double factor) {
        _xyz[0] *= factor;
        _xyz[1] *= factor;
        _xyz[2] *= factor;
        _has_length = false;
        _has_unit_vector = false;
    }
    
    void Vector::scale_up() {
        scale(GLOBALS::SCALE);
        _has_length = false;
    }
    void Vector::scale_down() {
        scale(1.0 / GLOBALS::SCALE);
        _has_length = false;
    }

    void Vector::rescale(double factor) {
        unitize();
        scale(factor);
        _has_length = false;
        
    }

std::string Vector::to_string() {
    std::ostringstream oss;

    oss << "geo::Vector "
        << _xyz[0] << " "
        << _xyz[1] << " "
        << _xyz[2] << " "
        << this->length();

    return oss.str();
}

} // namespace geo