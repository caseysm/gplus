#include "algorithm/kabsch.h"
#include "math/rotation.h"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

using namespace gangsta::algorithm;
using namespace gangsta::math;

class KabschTest : public ::testing::Test {
protected:
    Kabsch kabsch;
    std::vector<Vec3> setA;
    std::vector<Vec3> setB;
    std::vector<Vec3> setC;  // Rotated version of setA
    
    void SetUp() override {
        // Create a simple point set
        setA = {
            Vec3(1.0, 0.0, 0.0),
            Vec3(0.0, 1.0, 0.0),
            Vec3(0.0, 0.0, 1.0),
            Vec3(1.0, 1.0, 1.0)
        };
        
        // Create an identical point set with noise
        setB = {
            Vec3(1.1, 0.1, -0.1),
            Vec3(-0.1, 1.0, 0.1),
            Vec3(0.1, -0.1, 0.9),
            Vec3(0.9, 1.1, 1.0)
        };
        
        // Create a rotated version of setA (90 degrees around Z)
        setC = {
            Vec3(0.0, 1.0, 0.0),
            Vec3(-1.0, 0.0, 0.0),
            Vec3(0.0, 0.0, 1.0),
            Vec3(-1.0, 1.0, 1.0)
        };
    }
};

TEST_F(KabschTest, CalculateCentroid) {
    Vec3 centroid = kabsch.calculateCentroid(setA);
    
    // Expected centroid of setA: (0.5, 0.5, 0.5)
    EXPECT_DOUBLE_EQ(0.5, centroid.x);
    EXPECT_DOUBLE_EQ(0.5, centroid.y);
    EXPECT_DOUBLE_EQ(0.5, centroid.z);
}

TEST_F(KabschTest, CalculateRMSD) {
    // RMSD between identical sets should be 0
    double rmsd = kabsch.calculateRMSD(setA, setA);
    EXPECT_NEAR(0.0, rmsd, 1e-10);
    
    // RMSD between setA and setB should be small but non-zero
    rmsd = kabsch.calculateRMSD(setA, setB);
    EXPECT_GT(rmsd, 0.0);
    EXPECT_LT(rmsd, 0.3);  // Small deviation
    
    // RMSD between setA and setC should be larger since rotation is needed
    rmsd = kabsch.calculateRMSD(setA, setC);
    EXPECT_GT(rmsd, 0.5);
}

TEST_F(KabschTest, AlignIdentical) {
    Rotation rotation;
    Vec3 translation;
    
    // Aligning identical sets should give identity rotation and zero translation
    double rmsd = kabsch.align(setA, setA, rotation, translation);
    
    EXPECT_NEAR(0.0, rmsd, 1e-10);
    // Check rotation is identity by verifying quaternion
    Vec4 quat = rotation.getQuaternion();
    EXPECT_NEAR(1.0, quat.w, 1e-10);  // Identity quaternion
    EXPECT_NEAR(0.0, quat.x, 1e-10);
    EXPECT_NEAR(0.0, quat.y, 1e-10);
    EXPECT_NEAR(0.0, quat.z, 1e-10);
    EXPECT_NEAR(0.0, translation.x, 1e-10);
    EXPECT_NEAR(0.0, translation.y, 1e-10);
    EXPECT_NEAR(0.0, translation.z, 1e-10);
}

TEST_F(KabschTest, AlignWithNoise) {
    Rotation rotation;
    Vec3 translation;
    
    // Aligning setA and setB should give near-identity rotation and small translation
    double rmsd = kabsch.align(setA, setB, rotation, translation);
    
    EXPECT_LT(rmsd, 0.3);  // Small RMSD after alignment
    // Check rotation is close to identity
    Vec4 quat = rotation.getQuaternion();
    EXPECT_NEAR(1.0, quat.w, 0.1);  // Near-identity quaternion
    EXPECT_NEAR(0.0, quat.x, 0.1);
    EXPECT_NEAR(0.0, quat.y, 0.1);
    EXPECT_NEAR(0.0, quat.z, 0.1);
    // Small translations possible
    EXPECT_NEAR(0.0, translation.length(), 0.3);
}

TEST_F(KabschTest, AlignRotated) {
    Rotation rotation;
    Vec3 translation;
    
    // Aligning setA and setC should recover the 90-degree rotation around Z
    double rmsd = kabsch.align(setA, setC, rotation, translation);
    
    // Should align perfectly or very close
    EXPECT_NEAR(0.0, rmsd, 0.01);
    
    // Check the resulting rotation is approximately 90 degrees around Z
    // For 90-degree Z rotation, quaternion should be (w=cos(45°), x=0, y=0, z=sin(45°))
    Vec4 quat = rotation.getQuaternion();
    EXPECT_NEAR(std::cos(M_PI/4), std::abs(quat.w), 0.01);
    EXPECT_NEAR(0.0, quat.x, 0.01);
    EXPECT_NEAR(0.0, quat.y, 0.01);
    // The sign might vary depending on implementation
    EXPECT_NEAR(std::sin(M_PI/4), std::abs(quat.z), 0.01);
    
    // Translation should be minimal
    EXPECT_NEAR(0.0, translation.length(), 0.01);
}

TEST_F(KabschTest, AlignTranslated) {
    Rotation rotation;
    Vec3 translation;
    
    // Create a translated version of setA
    std::vector<Vec3> translatedSet;
    Vec3 offset(5.0, -3.0, 2.0);
    for (const Vec3& p : setA) {
        translatedSet.push_back(p + offset);
    }
    
    // Aligning setA and translatedSet should recover the translation
    double rmsd = kabsch.align(setA, translatedSet, rotation, translation);
    
    EXPECT_NEAR(0.0, rmsd, 1e-10);
    // Should be identity rotation
    Vec4 quat = rotation.getQuaternion();
    EXPECT_NEAR(1.0, quat.w, 1e-10);
    EXPECT_NEAR(0.0, quat.x, 1e-10);
    EXPECT_NEAR(0.0, quat.y, 1e-10);
    EXPECT_NEAR(0.0, quat.z, 1e-10);
    // Should recover the translation
    EXPECT_NEAR(offset.x, translation.x, 1e-10);
    EXPECT_NEAR(offset.y, translation.y, 1e-10);
    EXPECT_NEAR(offset.z, translation.z, 1e-10);
}