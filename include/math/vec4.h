#ifndef GANGSTA_VEC4_H
#define GANGSTA_VEC4_H

#include <cmath>

namespace gangsta {
namespace math {

/**
 * @brief 4D vector class for mathematical operations, often used for quaternions
 */
class Vec4 
{
public:
    /**
     * @brief Default constructor, initializes to zero
     */
    Vec4();
    
    /**
     * @brief Constructor with initial values
     * @param w W coordinate (typically the real part for quaternions)
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     */
    Vec4(double w, double x, double y, double z);
    
    /**
     * @brief Copy constructor
     */
    Vec4(const Vec4& other);
    
    /**
     * @brief Calculate dot product
     * @param other Another vector
     * @return Dot product value
     */
    double dot(const Vec4& other) const;
    
    /**
     * @brief Calculate vector length
     * @return Length of the vector
     */
    double length() const;
    
    /**
     * @brief Normalize vector to unit length
     * @return Normalized vector
     */
    Vec4 normalize() const;
    
    // Operators
    Vec4 operator+(const Vec4& other) const;
    Vec4 operator-(const Vec4& other) const;
    Vec4 operator*(double scalar) const;
    Vec4 operator/(double scalar) const;
    Vec4& operator+=(const Vec4& other);
    Vec4& operator-=(const Vec4& other);
    Vec4& operator*=(double scalar);
    Vec4& operator/=(double scalar);
    Vec4& operator=(const Vec4& other);
    bool operator==(const Vec4& other) const;
    bool operator!=(const Vec4& other) const;
    
    // Coordinate access
    double w, x, y, z;
};

// Non-member operator
Vec4 operator*(double scalar, const Vec4& vec);

} // namespace math
} // namespace gangsta

#endif // GANGSTA_VEC4_H