#include "math/rotation.h"
#include <cmath>

namespace gangsta {
namespace math {

Rotation::Rotation() : q(1.0, 0.0, 0.0, 0.0) 
{
}

Rotation::Rotation(const Vec4& q) : q(q) 
{
    normalize();
}

Rotation::Rotation(const Vec3& axis, double angle) 
{
    const double halfAngle = angle * 0.5;
    const double sinHalfAngle = std::sin(halfAngle);
    
    q.w = std::cos(halfAngle);
    q.x = axis.x * sinHalfAngle;
    q.y = axis.y * sinHalfAngle;
    q.z = axis.z * sinHalfAngle;
    
    normalize();
}

Rotation::Rotation(double pitch, double yaw, double roll) 
{
    // Convert Euler angles to quaternion using the ZYX convention
    const double halfPitch = pitch * 0.5;
    const double halfYaw = yaw * 0.5;
    const double halfRoll = roll * 0.5;
    
    const double cosHalfPitch = std::cos(halfPitch);
    const double sinHalfPitch = std::sin(halfPitch);
    const double cosHalfYaw = std::cos(halfYaw);
    const double sinHalfYaw = std::sin(halfYaw);
    const double cosHalfRoll = std::cos(halfRoll);
    const double sinHalfRoll = std::sin(halfRoll);
    
    q.w = cosHalfRoll * cosHalfPitch * cosHalfYaw + sinHalfRoll * sinHalfPitch * sinHalfYaw;
    q.x = cosHalfRoll * sinHalfPitch * cosHalfYaw + sinHalfRoll * cosHalfPitch * sinHalfYaw;
    q.y = cosHalfRoll * cosHalfPitch * sinHalfYaw - sinHalfRoll * sinHalfPitch * cosHalfYaw;
    q.z = sinHalfRoll * cosHalfPitch * cosHalfYaw - cosHalfRoll * sinHalfPitch * sinHalfYaw;
    
    normalize();
}

const Vec4& Rotation::getQuaternion() const 
{
    return q;
}

Mat4 Rotation::getMatrix() const 
{
    const double w = q.w;
    const double x = q.x;
    const double y = q.y;
    const double z = q.z;
    
    const double xx = x * x;
    const double xy = x * y;
    const double xz = x * z;
    const double xw = x * w;
    
    const double yy = y * y;
    const double yz = y * z;
    const double yw = y * w;
    
    const double zz = z * z;
    const double zw = z * w;
    
    std::array<double, 16> m = {
        1.0 - 2.0 * (yy + zz), 2.0 * (xy - zw), 2.0 * (xz + yw), 0.0,
        2.0 * (xy + zw), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - xw), 0.0,
        2.0 * (xz - yw), 2.0 * (yz + xw), 1.0 - 2.0 * (xx + yy), 0.0,
        0.0, 0.0, 0.0, 1.0
    };
    
    return Mat4(m);
}

Vec3 Rotation::rotate(const Vec3& point) const 
{
    // Using quaternion rotation formula: p' = q * p * q^-1
    // Where p is the point as a quaternion (0, x, y, z)
    
    const double w = q.w;
    const double x = q.x;
    const double y = q.y;
    const double z = q.z;
    
    // Vector part of the quaternion
    const Vec3 qv(x, y, z);
    
    // Apply quaternion formula v' = v + 2w(qv x v) + 2(qv x (qv x v))
    const Vec3 temp1 = qv.cross(point);
    const Vec3 temp2 = qv.cross(temp1);
    
    return point + 2.0 * w * temp1 + 2.0 * temp2;
}

Rotation Rotation::inverse() const 
{
    // Quaternion inverse is the conjugate (for unit quaternions)
    return Rotation(Vec4(q.w, -q.x, -q.y, -q.z));
}

Rotation Rotation::operator*(const Rotation& other) const 
{
    // Quaternion multiplication
    const double w1 = q.w;
    const double x1 = q.x;
    const double y1 = q.y;
    const double z1 = q.z;
    
    const double w2 = other.q.w;
    const double x2 = other.q.x;
    const double y2 = other.q.y;
    const double z2 = other.q.z;
    
    Vec4 result(
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    );
    
    return Rotation(result);
}

Rotation Rotation::slerp(const Rotation& other, double t) const 
{
    if (t <= 0.0) return *this;
    if (t >= 1.0) return other;
    
    // Calculate the dot product
    double dot = q.w * other.q.w + q.x * other.q.x + q.y * other.q.y + q.z * other.q.z;
    
    // Ensure shortest path
    Vec4 end = other.q;
    if (dot < 0.0) {
        end = Vec4(-end.w, -end.x, -end.y, -end.z);
        dot = -dot;
    }
    
    // If quaternions are very close, use linear interpolation
    if (dot > 0.9995) {
        Vec4 result = q * (1.0 - t) + end * t;
        return Rotation(result.normalize());
    }
    
    // Use spherical interpolation
    const double theta0 = std::acos(dot);
    const double theta = theta0 * t;
    
    const double sinTheta0 = std::sin(theta0);
    const double sinTheta = std::sin(theta);
    const double sinThetaComp = std::sin(theta0 - theta);
    
    Vec4 result = q * (sinThetaComp / sinTheta0) + end * (sinTheta / sinTheta0);
    return Rotation(result);
}

void Rotation::normalize() 
{
    const double length = q.length();
    if (length > 0.0) {
        q = q / length;
    } else {
        q = Vec4(1.0, 0.0, 0.0, 0.0);
    }
}

} // namespace math
} // namespace gangsta