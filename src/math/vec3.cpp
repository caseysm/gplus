#include "math/vec3.h"
#include <cmath>
#include <stdexcept>

namespace gangsta {
namespace math {

Vec3::Vec3() : x(0.0), y(0.0), z(0.0) {}

Vec3::Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

Vec3::Vec3(const Vec3& other) : x(other.x), y(other.y), z(other.z) {}

double Vec3::dist(const Vec3& other) const 
{
    const double dx = x - other.x;
    const double dy = y - other.y;
    const double dz = z - other.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double Vec3::dot(const Vec3& other) const 
{
    return x * other.x + y * other.y + z * other.z;
}

Vec3 Vec3::cross(const Vec3& other) const 
{
    return Vec3(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

double Vec3::length() const 
{
    return std::sqrt(x * x + y * y + z * z);
}

Vec3 Vec3::normalize() const 
{
    const double len = length();
    if (len < 1e-10) return *this;
    return *this / len;
}

double Vec3::operator[](int index) const 
{
    switch (index) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: throw std::out_of_range("Vector index out of range");
    }
}

double& Vec3::operator[](int index) 
{
    switch (index) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: throw std::out_of_range("Vector index out of range");
    }
}

// Operators
Vec3 Vec3::operator+(const Vec3& other) const 
{
    return Vec3(x + other.x, y + other.y, z + other.z);
}

Vec3 Vec3::operator-(const Vec3& other) const 
{
    return Vec3(x - other.x, y - other.y, z - other.z);
}

Vec3 Vec3::operator*(double scalar) const 
{
    return Vec3(x * scalar, y * scalar, z * scalar);
}

Vec3 Vec3::operator/(double scalar) const 
{
    if (std::abs(scalar) < 1e-10) return *this;
    const double invScalar = 1.0 / scalar;
    return Vec3(x * invScalar, y * invScalar, z * invScalar);
}

Vec3& Vec3::operator+=(const Vec3& other) 
{
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vec3& Vec3::operator-=(const Vec3& other) 
{
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vec3& Vec3::operator*=(double scalar) 
{
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

Vec3& Vec3::operator/=(double scalar) 
{
    if (std::abs(scalar) < 1e-10) return *this;
    const double invScalar = 1.0 / scalar;
    x *= invScalar;
    y *= invScalar;
    z *= invScalar;
    return *this;
}

Vec3& Vec3::operator=(const Vec3& other) 
{
    if (this != &other) {
        x = other.x;
        y = other.y;
        z = other.z;
    }
    return *this;
}

bool Vec3::operator==(const Vec3& other) const 
{
    const double epsilon = 1e-10;
    return (std::abs(x - other.x) < epsilon && 
            std::abs(y - other.y) < epsilon && 
            std::abs(z - other.z) < epsilon);
}

bool Vec3::operator!=(const Vec3& other) const 
{
    return !(*this == other);
}

// Non-member operator
Vec3 operator*(double scalar, const Vec3& vec) 
{
    return vec * scalar;
}

} // namespace math
} // namespace gangsta