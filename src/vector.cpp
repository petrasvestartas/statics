#include "Vector.hpp"

namespace geo
{

    Vector::Vector() : _xyz{0.0, 0.0, 0.0} {}

    Vector::Vector(double x, double y, double z){
        _xyz[0] = x;
        _xyz[1] = y;
        _xyz[2] = z;
        _has_length = false;
    }

    std::string Vector::to_string() const {
        std::ostringstream oss;
        oss << "(" << _xyz[0] << ", " << _xyz[1] << ", " << _xyz[2] << ")";
        return oss.str();
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
            // microscopic vectors that have infinte length!
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

    Vector Vector::component(Vector& other) {
        double cosa = this->dot(other);
        Vector component = Vector(other._xyz[0], other._xyz[1], other._xyz[2]);
        double L = component.length();
        component.scale(cosa / L);
        return component;
    }

    int Vector::is_parallel_to(Vector &v) {
        double ll = length() * v.length();
        int result;

        if (ll > 0.0)
        {
            const double cos_angle = (_xyz[0] * v[0] + _xyz[1] * v[1] + _xyz[2] * v[2]) / ll;

            const double angle_in_radians = GLOBALS::ANGLE * (GLOBALS::M_PI / 180.0); // convert angle from degrees to radians
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

    double Vector::angle(Vector &other) {
        double dot_product = this->dot(other);
        double magnitudes = this->length() * other.length();

        // Avoid division by zero
        if (magnitudes == 0) {
            return 0.0;
        }

        double cos_angle = dot_product / magnitudes;

        // Clamp the value to the range [-1, 1] to avoid errors with acos
        cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

        // Compute the angle in radians
        double angle_rad = std::acos(cos_angle);

        // Convert to degrees
        double angle_deg = angle_rad * (180.0 / GLOBALS::M_PI);

        return angle_deg;
    }

    Vector Vector::get_leveled_vector(double &vertical_height)
    {

        Vector copy = Vector(_xyz[0], _xyz[1], _xyz[2]);

        if (copy.unitize())
        {
            double angle = copy.angle(Vector(0, 0, 1));
            double inclined_offset_by_vertical_distance = vertical_height / std::cos(angle);
            copy *= inclined_offset_by_vertical_distance;
        }

        return copy;
    }

    double Vector::cosine_law(double& triangle_edge_length_a, double& triangle_edge_length_b, double& angle_in_degrees_between_edges){
        double angle_in_radians_between_edges = angle_in_degrees_between_edges * GLOBALS::M_PI / 180.0;
        return sqrt(pow(triangle_edge_length_a, 2) + pow(triangle_edge_length_b, 2) - 2*triangle_edge_length_a*triangle_edge_length_b*cos(angle_in_radians_between_edges));
    }

    double Vector::sine_law_angle(double& triangle_edge_length_a, double& angle_in_degrees_in_front_of_a, double& triangle_edge_length_b){
        double angle_in_radians_in_front_of_a = angle_in_degrees_in_front_of_a * GLOBALS::M_PI / 180.0;
        return asin( (triangle_edge_length_b * sin(angle_in_radians_in_front_of_a)) / triangle_edge_length_a ) * 180.0 / GLOBALS::M_PI;
    }

    double Vector::sine_law_length(double& triangle_edge_length_a, double& angle_in_degrees_in_front_of_a, double& angle_in_degrees_in_front_of_b){
        double angle_in_radians_in_front_of_a = angle_in_degrees_in_front_of_a * GLOBALS::M_PI / 180.0;
        double angle_in_radians_in_front_of_b = angle_in_degrees_in_front_of_b * GLOBALS::M_PI / 180.0;
        return (triangle_edge_length_a * sin(angle_in_radians_in_front_of_b))/sin(angle_in_radians_in_front_of_a);
    }

    void Vector::scale(double factor) {
        _xyz[0] *= factor;
        _xyz[1] *= factor;
        _xyz[2] *= factor;
    }
    
    void Vector::scale_up() {
        scale(GLOBALS::SCALE);
    }
    void Vector::scale_down() {
        scale(1.0 / GLOBALS::SCALE);
    }

    void Vector::rescale(double factor) {
        unitize();
        scale(factor);
    }

} // namespace geo