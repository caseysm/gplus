#include "math/vec4.h"
#include <cmath>

namespace gangsta {
namespace math {

Vec4::Vec4() : w(0.0), x(0.0), y(0.0), z(0.0) {}

Vec4::Vec4(double w, double x, double y, double z) : w(w), x(x), y(y), z(z) {}

Vec4::Vec4(const Vec4& other) : w(other.w), x(other.x), y(other.y), z(other.z) {}

double Vec4::dot(const Vec4& other) const 
{
    return w * other.w + x * other.x + y * other.y + z * other.z;
}

double Vec4::length() const 
{
    return std::sqrt(w * w + x * x + y * y + z * z);
}

Vec4 Vec4::normalize() const 
{
    const double len = length();
    if (len < 1e-10) return *this;
    return *this / len;
}

// Operators
Vec4 Vec4::operator+(const Vec4& other) const 
{
    return Vec4(w + other.w, x + other.x, y + other.y, z + other.z);
}

Vec4 Vec4::operator-(const Vec4& other) const 
{
    return Vec4(w - other.w, x - other.x, y - other.y, z - other.z);
}

Vec4 Vec4::operator*(double scalar) const 
{
    return Vec4(w * scalar, x * scalar, y * scalar, z * scalar);
}

Vec4 Vec4::operator/(double scalar) const 
{
    if (std::abs(scalar) < 1e-10) return *this;
    const double invScalar = 1.0 / scalar;
    return Vec4(w * invScalar, x * invScalar, y * invScalar, z * invScalar);
}

Vec4& Vec4::operator+=(const Vec4& other) 
{
    w += other.w;
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vec4& Vec4::operator-=(const Vec4& other) 
{
    w -= other.w;
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

Vec4& Vec4::operator*=(double scalar) 
{
    w *= scalar;
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

Vec4& Vec4::operator/=(double scalar) 
{
    if (std::abs(scalar) < 1e-10) return *this;
    const double invScalar = 1.0 / scalar;
    w *= invScalar;
    x *= invScalar;
    y *= invScalar;
    z *= invScalar;
    return *this;
}

Vec4& Vec4::operator=(const Vec4& other) 
{
    if (this != &other) {
        w = other.w;
        x = other.x;
        y = other.y;
        z = other.z;
    }
    return *this;
}

bool Vec4::operator==(const Vec4& other) const 
{
    const double epsilon = 1e-10;
    return (std::abs(w - other.w) < epsilon && 
            std::abs(x - other.x) < epsilon && 
            std::abs(y - other.y) < epsilon && 
            std::abs(z - other.z) < epsilon);
}

bool Vec4::operator!=(const Vec4& other) const 
{
    return !(*this == other);
}

// Non-member operator
Vec4 operator*(double scalar, const Vec4& vec) 
{
    return vec * scalar;
}

} // namespace math
} // namespace gangsta