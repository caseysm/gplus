#include "math/mat4.h"
#include <cmath>
#include <stdexcept>

namespace gangsta {
namespace math {

Mat4::Mat4() 
{
    identity();
}

Mat4::Mat4(const std::array<double, 16>& values) : m(values) 
{
}

Mat4::Mat4(const Mat4& other) : m(other.m) 
{
}

void Mat4::identity() 
{
    m = {1.0, 0.0, 0.0, 0.0,
         0.0, 1.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0};
}

void Mat4::setRotation(double angle, const Vec3& axis) 
{
    identity();
    
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double t = 1.0 - c;
    
    const double x = axis.x;
    const double y = axis.y;
    const double z = axis.z;
    
    // Row 0
    m[0] = t * x * x + c;
    m[1] = t * x * y - s * z;
    m[2] = t * x * z + s * y;
    
    // Row 1
    m[4] = t * x * y + s * z;
    m[5] = t * y * y + c;
    m[6] = t * y * z - s * x;
    
    // Row 2
    m[8] = t * x * z - s * y;
    m[9] = t * y * z + s * x;
    m[10] = t * z * z + c;
}

void Mat4::setTranslation(const Vec3& v) 
{
    identity();
    m[3] = v.x;
    m[7] = v.y;
    m[11] = v.z;
}

void Mat4::setScale(const Vec3& v) 
{
    identity();
    m[0] = v.x;
    m[5] = v.y;
    m[10] = v.z;
}

Mat4 Mat4::operator*(const Mat4& other) const 
{
    std::array<double, 16> result = {0};
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                result[i * 4 + j] += m[i * 4 + k] * other.m[k * 4 + j];
            }
        }
    }
    
    return Mat4(result);
}

Vec4 Mat4::operator*(const Vec4& v) const 
{
    return Vec4(
        m[0] * v.w + m[1] * v.x + m[2] * v.y + m[3] * v.z,
        m[4] * v.w + m[5] * v.x + m[6] * v.y + m[7] * v.z,
        m[8] * v.w + m[9] * v.x + m[10] * v.y + m[11] * v.z,
        m[12] * v.w + m[13] * v.x + m[14] * v.y + m[15] * v.z
    );
}

Vec3 Mat4::transformPoint(const Vec3& v) const 
{
    const double w = m[12] * v.x + m[13] * v.y + m[14] * v.z + m[15];
    
    if (std::abs(w) < 1e-10) {
        return Vec3(0, 0, 0);
    }
    
    const double invW = 1.0 / w;
    
    return Vec3(
        (m[0] * v.x + m[1] * v.y + m[2] * v.z + m[3]) * invW,
        (m[4] * v.x + m[5] * v.y + m[6] * v.z + m[7]) * invW,
        (m[8] * v.x + m[9] * v.y + m[10] * v.z + m[11]) * invW
    );
}

Vec3 Mat4::transformVector(const Vec3& v) const 
{
    return Vec3(
        m[0] * v.x + m[1] * v.y + m[2] * v.z,
        m[4] * v.x + m[5] * v.y + m[6] * v.z,
        m[8] * v.x + m[9] * v.y + m[10] * v.z
    );
}

double Mat4::determinant() const 
{
    return 
        m[12] * m[9] * m[6] * m[3] - m[8] * m[13] * m[6] * m[3] -
        m[12] * m[5] * m[10] * m[3] + m[4] * m[13] * m[10] * m[3] +
        m[8] * m[5] * m[14] * m[3] - m[4] * m[9] * m[14] * m[3] -
        m[12] * m[9] * m[2] * m[7] + m[8] * m[13] * m[2] * m[7] +
        m[12] * m[1] * m[10] * m[7] - m[0] * m[13] * m[10] * m[7] -
        m[8] * m[1] * m[14] * m[7] + m[0] * m[9] * m[14] * m[7] +
        m[12] * m[5] * m[2] * m[11] - m[4] * m[13] * m[2] * m[11] -
        m[12] * m[1] * m[6] * m[11] + m[0] * m[13] * m[6] * m[11] +
        m[4] * m[1] * m[14] * m[11] - m[0] * m[5] * m[14] * m[11] -
        m[8] * m[5] * m[2] * m[15] + m[4] * m[9] * m[2] * m[15] +
        m[8] * m[1] * m[6] * m[15] - m[0] * m[9] * m[6] * m[15] -
        m[4] * m[1] * m[10] * m[15] + m[0] * m[5] * m[10] * m[15];
}

Mat4 Mat4::inverse() const 
{
    const double det = determinant();
    
    if (std::abs(det) < 1e-10) {
        throw std::runtime_error("Matrix is not invertible");
    }
    
    const double invDet = 1.0 / det;
    
    std::array<double, 16> inv = {0};
    
    inv[0] = (m[5] * (m[10] * m[15] - m[14] * m[11]) + 
              m[9] * (m[14] * m[7] - m[6] * m[15]) + 
              m[13] * (m[6] * m[11] - m[10] * m[7])) * invDet;
    
    inv[1] = (m[1] * (m[14] * m[11] - m[10] * m[15]) + 
              m[9] * (m[2] * m[15] - m[14] * m[3]) + 
              m[13] * (m[10] * m[3] - m[2] * m[11])) * invDet;
    
    inv[2] = (m[1] * (m[6] * m[15] - m[14] * m[7]) + 
              m[5] * (m[14] * m[3] - m[2] * m[15]) + 
              m[13] * (m[2] * m[7] - m[6] * m[3])) * invDet;
    
    inv[3] = (m[1] * (m[10] * m[7] - m[6] * m[11]) + 
              m[5] * (m[2] * m[11] - m[10] * m[3]) + 
              m[9] * (m[6] * m[3] - m[2] * m[7])) * invDet;
    
    inv[4] = (m[4] * (m[14] * m[11] - m[10] * m[15]) + 
              m[8] * (m[6] * m[15] - m[14] * m[7]) + 
              m[12] * (m[10] * m[7] - m[6] * m[11])) * invDet;
    
    inv[5] = (m[0] * (m[10] * m[15] - m[14] * m[11]) + 
              m[8] * (m[14] * m[3] - m[2] * m[15]) + 
              m[12] * (m[2] * m[11] - m[10] * m[3])) * invDet;
    
    inv[6] = (m[0] * (m[14] * m[7] - m[6] * m[15]) + 
              m[4] * (m[2] * m[15] - m[14] * m[3]) + 
              m[12] * (m[6] * m[3] - m[2] * m[7])) * invDet;
    
    inv[7] = (m[0] * (m[6] * m[11] - m[10] * m[7]) + 
              m[4] * (m[10] * m[3] - m[2] * m[11]) + 
              m[8] * (m[2] * m[7] - m[6] * m[3])) * invDet;
    
    inv[8] = (m[4] * (m[9] * m[15] - m[13] * m[11]) + 
              m[8] * (m[13] * m[7] - m[5] * m[15]) + 
              m[12] * (m[5] * m[11] - m[9] * m[7])) * invDet;
    
    inv[9] = (m[0] * (m[13] * m[11] - m[9] * m[15]) + 
              m[8] * (m[1] * m[15] - m[13] * m[3]) + 
              m[12] * (m[9] * m[3] - m[1] * m[11])) * invDet;
    
    inv[10] = (m[0] * (m[5] * m[15] - m[13] * m[7]) + 
               m[4] * (m[13] * m[3] - m[1] * m[15]) + 
               m[12] * (m[1] * m[7] - m[5] * m[3])) * invDet;
    
    inv[11] = (m[0] * (m[9] * m[7] - m[5] * m[11]) + 
               m[4] * (m[1] * m[11] - m[9] * m[3]) + 
               m[8] * (m[5] * m[3] - m[1] * m[7])) * invDet;
    
    inv[12] = (m[4] * (m[13] * m[10] - m[9] * m[14]) + 
               m[8] * (m[5] * m[14] - m[13] * m[6]) + 
               m[12] * (m[9] * m[6] - m[5] * m[10])) * invDet;
    
    inv[13] = (m[0] * (m[9] * m[14] - m[13] * m[10]) + 
               m[8] * (m[13] * m[2] - m[1] * m[14]) + 
               m[12] * (m[1] * m[10] - m[9] * m[2])) * invDet;
    
    inv[14] = (m[0] * (m[13] * m[6] - m[5] * m[14]) + 
               m[4] * (m[1] * m[14] - m[13] * m[2]) + 
               m[12] * (m[5] * m[2] - m[1] * m[6])) * invDet;
    
    inv[15] = (m[0] * (m[5] * m[10] - m[9] * m[6]) + 
               m[4] * (m[9] * m[2] - m[1] * m[10]) + 
               m[8] * (m[1] * m[6] - m[5] * m[2])) * invDet;
    
    return Mat4(inv);
}

Mat4 Mat4::transpose() const 
{
    std::array<double, 16> trans = {
        m[0], m[4], m[8], m[12],
        m[1], m[5], m[9], m[13],
        m[2], m[6], m[10], m[14],
        m[3], m[7], m[11], m[15]
    };
    
    return Mat4(trans);
}

double& Mat4::at(int row, int col) 
{
    if (row < 0 || row >= 4 || col < 0 || col >= 4) {
        throw std::out_of_range("Matrix indices out of range");
    }
    return m[row * 4 + col];
}

double Mat4::at(int row, int col) const 
{
    if (row < 0 || row >= 4 || col < 0 || col >= 4) {
        throw std::out_of_range("Matrix indices out of range");
    }
    return m[row * 4 + col];
}

double* Mat4::data() 
{
    return m.data();
}

const double* Mat4::data() const 
{
    return m.data();
}

} // namespace math
} // namespace gangsta