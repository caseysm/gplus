#ifndef GANGSTA_MAT4_H
#define GANGSTA_MAT4_H

#include "math/vec3.h"
#include "math/vec4.h"
#include <array>

namespace gangsta {
namespace math {

class Vec4;

/**
 * @brief 4x4 matrix class for transformations
 */
class Mat4 
{
public:
    /**
     * @brief Default constructor, initializes to identity matrix
     */
    Mat4();
    
    /**
     * @brief Construct with provided array of 16 values
     * @param values The 16 values in row-major order
     */
    Mat4(const std::array<double, 16>& values);
    
    /**
     * @brief Copy constructor
     */
    Mat4(const Mat4& other);
    
    /**
     * @brief Set matrix to identity
     */
    void identity();
    
    /**
     * @brief Set as rotation matrix
     * @param angle Rotation angle in radians
     * @param axis Axis of rotation (must be normalized)
     */
    void setRotation(double angle, const Vec3& axis);
    
    /**
     * @brief Set as translation matrix
     * @param v Translation vector
     */
    void setTranslation(const Vec3& v);
    
    /**
     * @brief Set as scale matrix
     * @param v Scale factors for x, y, z
     */
    void setScale(const Vec3& v);
    
    /**
     * @brief Multiply with another matrix
     * @param other Matrix to multiply with
     * @return Result of multiplication
     */
    Mat4 operator*(const Mat4& other) const;
    
    /**
     * @brief Multiply with vector
     * @param v Vector to transform
     * @return Transformed vector
     */
    Vec4 operator*(const Vec4& v) const;
    
    /**
     * @brief Transform a 3D point (with implied w=1)
     * @param v The point to transform
     * @return Transformed point
     */
    Vec3 transformPoint(const Vec3& v) const;
    
    /**
     * @brief Transform a 3D vector (with implied w=0)
     * @param v The vector to transform
     * @return Transformed vector
     */
    Vec3 transformVector(const Vec3& v) const;
    
    /**
     * @brief Get determinant
     * @return Determinant of the matrix
     */
    double determinant() const;
    
    /**
     * @brief Get inverse matrix
     * @return Inverse matrix
     */
    Mat4 inverse() const;
    
    /**
     * @brief Get transpose
     * @return Transposed matrix
     */
    Mat4 transpose() const;
    
    /**
     * @brief Access element by index
     * @param row Row index (0-3)
     * @param col Column index (0-3)
     * @return Reference to element
     */
    double& at(int row, int col);
    
    /**
     * @brief Access element by index (const version)
     * @param row Row index (0-3)
     * @param col Column index (0-3)
     * @return Element value
     */
    double at(int row, int col) const;
    
    /**
     * @brief Get pointer to internal data
     * @return Pointer to array of 16 doubles
     */
    double* data();
    
    /**
     * @brief Get pointer to internal data (const version)
     * @return Pointer to array of 16 doubles
     */
    const double* data() const;
    
private:
    // Matrix data in row-major order (m00, m01, m02, m03, m10, ..., m33)
    std::array<double, 16> m;
};

} // namespace math
} // namespace gangsta

#endif // GANGSTA_MAT4_H