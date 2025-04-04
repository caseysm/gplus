#ifndef GANGSTA_VEC3_H
#define GANGSTA_VEC3_H

#include <cmath>

namespace gangsta {
namespace math {

/**
 * @brief 3D vector class for mathematical operations
 */
class Vec3 
{
public:
    /**
     * @brief Default constructor, initializes to zero
     */
    Vec3();
    
    /**
     * @brief Constructor with initial values
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     */
    Vec3(double x, double y, double z);
    
    /**
     * @brief Copy constructor
     */
    Vec3(const Vec3& other);
    
    /**
     * @brief Calculate distance between two vectors
     * @param other Another vector
     * @return Euclidean distance
     */
    double dist(const Vec3& other) const;
    
    /**
     * @brief Calculate dot product
     * @param other Another vector
     * @return Dot product value
     */
    double dot(const Vec3& other) const;
    
    /**
     * @brief Calculate cross product
     * @param other Another vector
     * @return Cross product vector
     */
    Vec3 cross(const Vec3& other) const;
    
    /**
     * @brief Calculate vector length
     * @return Length of the vector
     */
    double length() const;
    
    /**
     * @brief Normalize vector to unit length
     * @return Normalized vector
     */
    Vec3 normalize() const;
    
    /**
     * @brief Access element by index
     * @param index Index of element (0-2)
     * @return Element value
     */
    double operator[](int index) const;
    
    /**
     * @brief Access element by index (mutable)
     * @param index Index of element (0-2)
     * @return Reference to element
     */
    double& operator[](int index);
    
    // Operators
    Vec3 operator+(const Vec3& other) const;
    Vec3 operator-(const Vec3& other) const;
    Vec3 operator*(double scalar) const;
    Vec3 operator/(double scalar) const;
    Vec3& operator+=(const Vec3& other);
    Vec3& operator-=(const Vec3& other);
    Vec3& operator*=(double scalar);
    Vec3& operator/=(double scalar);
    Vec3& operator=(const Vec3& other);
    bool operator==(const Vec3& other) const;
    bool operator!=(const Vec3& other) const;
    
    // Coordinate access
    double x, y, z;
};

// Non-member operator
Vec3 operator*(double scalar, const Vec3& vec);

} // namespace math
} // namespace gangsta

#endif // GANGSTA_VEC3_H