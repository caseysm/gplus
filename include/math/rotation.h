#ifndef GANGSTA_ROTATION_H
#define GANGSTA_ROTATION_H

#include "math/vec3.h"
#include "math/vec4.h"
#include "math/mat4.h"

namespace gangsta {
namespace math {

/**
 * @brief Rotation representation using quaternions
 */
class Rotation 
{
public:
    /**
     * @brief Default constructor (identity rotation)
     */
    Rotation();
    
    /**
     * @brief Construct from a quaternion
     * @param q Quaternion (w, x, y, z) where w is the real part
     */
    explicit Rotation(const Vec4& q);
    
    /**
     * @brief Construct from axis and angle
     * @param axis Rotation axis (must be normalized)
     * @param angle Rotation angle in radians
     */
    Rotation(const Vec3& axis, double angle);
    
    /**
     * @brief Construct from Euler angles
     * @param pitch Rotation around X axis in radians
     * @param yaw Rotation around Y axis in radians
     * @param roll Rotation around Z axis in radians
     */
    Rotation(double pitch, double yaw, double roll);
    
    /**
     * @brief Get quaternion representation
     * @return Quaternion as Vec4
     */
    const Vec4& getQuaternion() const;
    
    /**
     * @brief Get rotation matrix
     * @return 4x4 rotation matrix
     */
    Mat4 getMatrix() const;
    
    /**
     * @brief Rotate a point
     * @param point The point to rotate
     * @return Rotated point
     */
    Vec3 rotate(const Vec3& point) const;
    
    /**
     * @brief Inverse rotation
     * @return New rotation that is the inverse of this one
     */
    Rotation inverse() const;
    
    /**
     * @brief Combine two rotations
     * @param other Another rotation to apply after this one
     * @return Combined rotation
     */
    Rotation operator*(const Rotation& other) const;
    
    /**
     * @brief Linearly interpolate between rotations
     * @param other Target rotation
     * @param t Interpolation factor (0-1)
     * @return Interpolated rotation
     */
    Rotation slerp(const Rotation& other, double t) const;
    
    /**
     * @brief Normalize the rotation quaternion
     */
    void normalize();
    
private:
    Vec4 q; // Quaternion (w, x, y, z)
};

} // namespace math
} // namespace gangsta

#endif // GANGSTA_ROTATION_H