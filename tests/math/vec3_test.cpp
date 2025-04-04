#include "math/vec3.h"
#include <gtest/gtest.h>
#include <cmath>

using namespace gangsta::math;

class Vec3Test : public ::testing::Test {
protected:
    Vec3 zero;
    Vec3 unit_x;
    Vec3 unit_y;
    Vec3 unit_z;
    Vec3 vec1;
    
    void SetUp() override {
        zero = Vec3(0, 0, 0);
        unit_x = Vec3(1, 0, 0);
        unit_y = Vec3(0, 1, 0);
        unit_z = Vec3(0, 0, 1);
        vec1 = Vec3(1, 2, 3);
    }
};

// Constructor tests
TEST_F(Vec3Test, DefaultConstructor) {
    Vec3 v;
    EXPECT_DOUBLE_EQ(0.0, v.x);
    EXPECT_DOUBLE_EQ(0.0, v.y);
    EXPECT_DOUBLE_EQ(0.0, v.z);
}

TEST_F(Vec3Test, ParameterizedConstructor) {
    Vec3 v(2.5, -1.5, 3.0);
    EXPECT_DOUBLE_EQ(2.5, v.x);
    EXPECT_DOUBLE_EQ(-1.5, v.y);
    EXPECT_DOUBLE_EQ(3.0, v.z);
}

TEST_F(Vec3Test, CopyConstructor) {
    Vec3 original(2.5, -1.5, 3.0);
    Vec3 copy(original);
    EXPECT_DOUBLE_EQ(original.x, copy.x);
    EXPECT_DOUBLE_EQ(original.y, copy.y);
    EXPECT_DOUBLE_EQ(original.z, copy.z);
}

// Basic operations tests
TEST_F(Vec3Test, Distance) {
    EXPECT_DOUBLE_EQ(0.0, zero.dist(zero));
    EXPECT_DOUBLE_EQ(1.0, zero.dist(unit_x));
    EXPECT_DOUBLE_EQ(std::sqrt(14.0), zero.dist(vec1));
    
    // Calculate the distance between (1,0,0) and (1,2,3) -> sqrt((0)^2 + (2)^2 + (3)^2) = sqrt(13)
    EXPECT_DOUBLE_EQ(std::sqrt(13.0), unit_x.dist(vec1));
}

TEST_F(Vec3Test, DotProduct) {
    EXPECT_DOUBLE_EQ(0.0, unit_x.dot(unit_y));
    EXPECT_DOUBLE_EQ(1.0, unit_x.dot(unit_x));
    EXPECT_DOUBLE_EQ(1.0, vec1.dot(unit_x));
    EXPECT_DOUBLE_EQ(2.0, vec1.dot(unit_y));
    EXPECT_DOUBLE_EQ(3.0, vec1.dot(unit_z));
    EXPECT_DOUBLE_EQ(14.0, vec1.dot(vec1));  // Length squared
}

TEST_F(Vec3Test, CrossProduct) {
    Vec3 result = unit_x.cross(unit_y);
    EXPECT_DOUBLE_EQ(0.0, result.x);
    EXPECT_DOUBLE_EQ(0.0, result.y);
    EXPECT_DOUBLE_EQ(1.0, result.z);
    
    result = unit_y.cross(unit_x);
    EXPECT_DOUBLE_EQ(0.0, result.x);
    EXPECT_DOUBLE_EQ(0.0, result.y);
    EXPECT_DOUBLE_EQ(-1.0, result.z);
    
    result = vec1.cross(unit_x);
    EXPECT_DOUBLE_EQ(0.0, result.x);
    EXPECT_DOUBLE_EQ(3.0, result.y);
    EXPECT_DOUBLE_EQ(-2.0, result.z);
}

TEST_F(Vec3Test, Length) {
    EXPECT_DOUBLE_EQ(0.0, zero.length());
    EXPECT_DOUBLE_EQ(1.0, unit_x.length());
    EXPECT_DOUBLE_EQ(std::sqrt(14.0), vec1.length());
}

TEST_F(Vec3Test, Normalize) {
    Vec3 normalized = vec1.normalize();
    double length = vec1.length();
    EXPECT_DOUBLE_EQ(1.0/length, normalized.x/vec1.x);
    EXPECT_DOUBLE_EQ(1.0/length, normalized.y/vec1.y);
    EXPECT_DOUBLE_EQ(1.0/length, normalized.z/vec1.z);
    EXPECT_DOUBLE_EQ(1.0, normalized.length());
}

// Operator tests
TEST_F(Vec3Test, IndexOperator) {
    EXPECT_DOUBLE_EQ(vec1.x, vec1[0]);
    EXPECT_DOUBLE_EQ(vec1.y, vec1[1]);
    EXPECT_DOUBLE_EQ(vec1.z, vec1[2]);
    
    Vec3 v = vec1;
    v[0] = 10.0;
    v[1] = 20.0;
    v[2] = 30.0;
    EXPECT_DOUBLE_EQ(10.0, v.x);
    EXPECT_DOUBLE_EQ(20.0, v.y);
    EXPECT_DOUBLE_EQ(30.0, v.z);
}

TEST_F(Vec3Test, AdditionOperator) {
    Vec3 result = vec1 + unit_x;
    EXPECT_DOUBLE_EQ(2.0, result.x);
    EXPECT_DOUBLE_EQ(2.0, result.y);
    EXPECT_DOUBLE_EQ(3.0, result.z);
}

TEST_F(Vec3Test, SubtractionOperator) {
    Vec3 result = vec1 - unit_x;
    EXPECT_DOUBLE_EQ(0.0, result.x);
    EXPECT_DOUBLE_EQ(2.0, result.y);
    EXPECT_DOUBLE_EQ(3.0, result.z);
}

TEST_F(Vec3Test, MultiplicationOperator) {
    Vec3 result = vec1 * 2.0;
    EXPECT_DOUBLE_EQ(2.0, result.x);
    EXPECT_DOUBLE_EQ(4.0, result.y);
    EXPECT_DOUBLE_EQ(6.0, result.z);
    
    result = 3.0 * vec1;
    EXPECT_DOUBLE_EQ(3.0, result.x);
    EXPECT_DOUBLE_EQ(6.0, result.y);
    EXPECT_DOUBLE_EQ(9.0, result.z);
}

TEST_F(Vec3Test, DivisionOperator) {
    Vec3 result = vec1 / 2.0;
    EXPECT_DOUBLE_EQ(0.5, result.x);
    EXPECT_DOUBLE_EQ(1.0, result.y);
    EXPECT_DOUBLE_EQ(1.5, result.z);
}

TEST_F(Vec3Test, CompoundOperators) {
    Vec3 v = vec1;
    v += unit_x;
    EXPECT_DOUBLE_EQ(2.0, v.x);
    EXPECT_DOUBLE_EQ(2.0, v.y);
    EXPECT_DOUBLE_EQ(3.0, v.z);
    
    v = vec1;
    v -= unit_x;
    EXPECT_DOUBLE_EQ(0.0, v.x);
    EXPECT_DOUBLE_EQ(2.0, v.y);
    EXPECT_DOUBLE_EQ(3.0, v.z);
    
    v = vec1;
    v *= 2.0;
    EXPECT_DOUBLE_EQ(2.0, v.x);
    EXPECT_DOUBLE_EQ(4.0, v.y);
    EXPECT_DOUBLE_EQ(6.0, v.z);
    
    v = vec1;
    v /= 2.0;
    EXPECT_DOUBLE_EQ(0.5, v.x);
    EXPECT_DOUBLE_EQ(1.0, v.y);
    EXPECT_DOUBLE_EQ(1.5, v.z);
}

TEST_F(Vec3Test, EqualityOperators) {
    EXPECT_TRUE(vec1 == vec1);
    EXPECT_FALSE(vec1 == unit_x);
    EXPECT_FALSE(vec1 != vec1);
    EXPECT_TRUE(vec1 != unit_x);
}