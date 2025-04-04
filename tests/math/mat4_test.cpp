#include "math/mat4.h"
#include "math/vec3.h"
#include "math/vec4.h"
#include <gtest/gtest.h>
#include <cmath>
#include <array>

using namespace gangsta::math;

class Mat4Test : public ::testing::Test {
protected:
    Mat4 identity;
    Vec3 point;
    Vec3 vec;
    
    void SetUp() override {
        identity = Mat4(); // Default constructor creates identity matrix
        point = Vec3(1, 2, 3);
        vec = Vec3(1, 0, 0);
    }
    
    // Helper function to check if a matrix is identity
    void ExpectIdentity(const Mat4& mat) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) {
                    EXPECT_DOUBLE_EQ(1.0, mat.at(i, j));
                } else {
                    EXPECT_DOUBLE_EQ(0.0, mat.at(i, j));
                }
            }
        }
    }
};

// Constructor tests
TEST_F(Mat4Test, DefaultConstructor) {
    Mat4 mat;
    ExpectIdentity(mat);
}

TEST_F(Mat4Test, ArrayConstructor) {
    std::array<double, 16> values = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };
    
    Mat4 mat(values);
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_DOUBLE_EQ(values[i * 4 + j], mat.at(i, j));
        }
    }
}

TEST_F(Mat4Test, CopyConstructor) {
    std::array<double, 16> values = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };
    
    Mat4 original(values);
    Mat4 copy(original);
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_DOUBLE_EQ(original.at(i, j), copy.at(i, j));
        }
    }
}

// Basic operations tests
TEST_F(Mat4Test, Identity) {
    Mat4 mat;
    mat.at(0, 0) = 5;
    mat.identity();
    ExpectIdentity(mat);
}

TEST_F(Mat4Test, TranslationMatrix) {
    Vec3 translation(5, -3, 2);
    Mat4 mat;
    mat.setTranslation(translation);
    
    // Check that the matrix has the correct form for translation
    EXPECT_DOUBLE_EQ(1.0, mat.at(0, 0));
    EXPECT_DOUBLE_EQ(0.0, mat.at(0, 1));
    EXPECT_DOUBLE_EQ(0.0, mat.at(0, 2));
    EXPECT_DOUBLE_EQ(translation.x, mat.at(0, 3));
    
    EXPECT_DOUBLE_EQ(0.0, mat.at(1, 0));
    EXPECT_DOUBLE_EQ(1.0, mat.at(1, 1));
    EXPECT_DOUBLE_EQ(0.0, mat.at(1, 2));
    EXPECT_DOUBLE_EQ(translation.y, mat.at(1, 3));
    
    EXPECT_DOUBLE_EQ(0.0, mat.at(2, 0));
    EXPECT_DOUBLE_EQ(0.0, mat.at(2, 1));
    EXPECT_DOUBLE_EQ(1.0, mat.at(2, 2));
    EXPECT_DOUBLE_EQ(translation.z, mat.at(2, 3));
    
    EXPECT_DOUBLE_EQ(0.0, mat.at(3, 0));
    EXPECT_DOUBLE_EQ(0.0, mat.at(3, 1));
    EXPECT_DOUBLE_EQ(0.0, mat.at(3, 2));
    EXPECT_DOUBLE_EQ(1.0, mat.at(3, 3));
    
    // Test transforming a point
    Vec3 transformed = mat.transformPoint(point);
    EXPECT_DOUBLE_EQ(point.x + translation.x, transformed.x);
    EXPECT_DOUBLE_EQ(point.y + translation.y, transformed.y);
    EXPECT_DOUBLE_EQ(point.z + translation.z, transformed.z);
}

TEST_F(Mat4Test, ScaleMatrix) {
    Vec3 scale(2, 3, 4);
    Mat4 mat;
    mat.setScale(scale);
    
    // Check that the matrix has the correct form for scaling
    EXPECT_DOUBLE_EQ(scale.x, mat.at(0, 0));
    EXPECT_DOUBLE_EQ(0.0, mat.at(0, 1));
    EXPECT_DOUBLE_EQ(0.0, mat.at(0, 2));
    EXPECT_DOUBLE_EQ(0.0, mat.at(0, 3));
    
    EXPECT_DOUBLE_EQ(0.0, mat.at(1, 0));
    EXPECT_DOUBLE_EQ(scale.y, mat.at(1, 1));
    EXPECT_DOUBLE_EQ(0.0, mat.at(1, 2));
    EXPECT_DOUBLE_EQ(0.0, mat.at(1, 3));
    
    EXPECT_DOUBLE_EQ(0.0, mat.at(2, 0));
    EXPECT_DOUBLE_EQ(0.0, mat.at(2, 1));
    EXPECT_DOUBLE_EQ(scale.z, mat.at(2, 2));
    EXPECT_DOUBLE_EQ(0.0, mat.at(2, 3));
    
    EXPECT_DOUBLE_EQ(0.0, mat.at(3, 0));
    EXPECT_DOUBLE_EQ(0.0, mat.at(3, 1));
    EXPECT_DOUBLE_EQ(0.0, mat.at(3, 2));
    EXPECT_DOUBLE_EQ(1.0, mat.at(3, 3));
    
    // Test transforming a point
    Vec3 transformed = mat.transformPoint(point);
    EXPECT_DOUBLE_EQ(point.x * scale.x, transformed.x);
    EXPECT_DOUBLE_EQ(point.y * scale.y, transformed.y);
    EXPECT_DOUBLE_EQ(point.z * scale.z, transformed.z);
}

TEST_F(Mat4Test, RotationMatrix) {
    // Test rotation around Z axis by 90 degrees
    double angle = M_PI / 2.0; // 90 degrees
    Vec3 axis(0, 0, 1);  // Z axis
    
    Mat4 mat;
    mat.setRotation(angle, axis);
    
    Vec3 point(1, 0, 0);
    Vec3 transformed = mat.transformPoint(point);
    
    // A point (1,0,0) rotated 90 degrees around Z should be approximately (0,1,0)
    EXPECT_NEAR(0.0, transformed.x, 1e-10);
    EXPECT_NEAR(1.0, transformed.y, 1e-10);
    EXPECT_NEAR(0.0, transformed.z, 1e-10);
}

TEST_F(Mat4Test, MatrixMultiplication) {
    // Create translation and scale matrices
    Mat4 translation;
    translation.setTranslation(Vec3(1, 2, 3));
    
    Mat4 scale;
    scale.setScale(Vec3(2, 2, 2));
    
    // Multiply them together
    Mat4 combined = translation * scale;
    
    // The result should first scale, then translate
    Vec3 transformed = combined.transformPoint(point);
    Vec3 expected(point.x * 2 + 1, point.y * 2 + 2, point.z * 2 + 3);
    
    EXPECT_DOUBLE_EQ(expected.x, transformed.x);
    EXPECT_DOUBLE_EQ(expected.y, transformed.y);
    EXPECT_DOUBLE_EQ(expected.z, transformed.z);
}

TEST_F(Mat4Test, VectorTransformation) {
    // Create translation matrix
    Mat4 translation;
    translation.setTranslation(Vec3(1, 2, 3));
    
    // Vectors should not be affected by translation
    Vec3 transformed = translation.transformVector(vec);
    EXPECT_DOUBLE_EQ(vec.x, transformed.x);
    EXPECT_DOUBLE_EQ(vec.y, transformed.y);
    EXPECT_DOUBLE_EQ(vec.z, transformed.z);
    
    // Create rotation matrix - 90 degrees around Z
    Mat4 rotation;
    rotation.setRotation(M_PI / 2.0, Vec3(0, 0, 1));
    
    // Vector (1,0,0) rotated 90 degrees around Z should be approximately (0,1,0)
    transformed = rotation.transformVector(vec);
    EXPECT_NEAR(0.0, transformed.x, 1e-10);
    EXPECT_NEAR(1.0, transformed.y, 1e-10);
    EXPECT_NEAR(0.0, transformed.z, 1e-10);
}

TEST_F(Mat4Test, InverseAndTranspose) {
    // Test inverse of identity
    Mat4 invIdentity = identity.inverse();
    ExpectIdentity(invIdentity);
    
    // Test inverse of translation
    Mat4 translation;
    translation.setTranslation(Vec3(1, 2, 3));
    Mat4 invTranslation = translation.inverse();
    
    // The inverse of a translation by (x,y,z) is a translation by (-x,-y,-z)
    EXPECT_DOUBLE_EQ(-1.0, invTranslation.at(0, 3));
    EXPECT_DOUBLE_EQ(-2.0, invTranslation.at(1, 3));
    EXPECT_DOUBLE_EQ(-3.0, invTranslation.at(2, 3));
    
    // Test inverse of scale
    Mat4 scale;
    scale.setScale(Vec3(2, 3, 4));
    Mat4 invScale = scale.inverse();
    
    // The inverse of scaling by (x,y,z) is scaling by (1/x,1/y,1/z)
    EXPECT_DOUBLE_EQ(1.0/2.0, invScale.at(0, 0));
    EXPECT_DOUBLE_EQ(1.0/3.0, invScale.at(1, 1));
    EXPECT_DOUBLE_EQ(1.0/4.0, invScale.at(2, 2));
    
    // Test transpose
    std::array<double, 16> values = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };
    
    Mat4 mat(values);
    Mat4 transposed = mat.transpose();
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            EXPECT_DOUBLE_EQ(mat.at(i, j), transposed.at(j, i));
        }
    }
}

TEST_F(Mat4Test, Determinant) {
    // Identity determinant should be 1
    EXPECT_DOUBLE_EQ(1.0, identity.determinant());
    
    // Scale determinant should be the product of scale factors
    Mat4 scale;
    scale.setScale(Vec3(2, 3, 4));
    EXPECT_DOUBLE_EQ(2.0 * 3.0 * 4.0, scale.determinant());
}