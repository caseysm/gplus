#include "core/residue.h"
#include <gtest/gtest.h>
#include <string>

using namespace gangsta::core;
using namespace gangsta::math;

class ResidueTest : public ::testing::Test {
protected:
    Residue residue;
    
    void SetUp() override {
        residue = Residue("ALA", 1, 'A');
        
        // Add some atoms to the residue
        residue.addAtom(Atom(Vec3(1.0, 2.0, 3.0), "CA", "ALA", 1, 'A'));
        residue.addAtom(Atom(Vec3(2.0, 3.0, 4.0), "CB", "ALA", 1, 'A'));
        residue.addAtom(Atom(Vec3(3.0, 4.0, 5.0), "N", "ALA", 1, 'A'));
        residue.addAtom(Atom(Vec3(4.0, 5.0, 6.0), "C", "ALA", 1, 'A'));
        residue.addAtom(Atom(Vec3(5.0, 6.0, 7.0), "O", "ALA", 1, 'A'));
    }
};

// Constructor tests
TEST_F(ResidueTest, DefaultConstructor) {
    Residue defaultResidue;
    EXPECT_EQ("UNK", defaultResidue.getName()); // Default name is "UNK" not empty string
    EXPECT_EQ(0, defaultResidue.getNumber());
    EXPECT_EQ(' ', defaultResidue.getChainID());
    EXPECT_EQ(-1, defaultResidue.getSSEIndex());
    EXPECT_TRUE(defaultResidue.getAtoms().empty());
}

TEST_F(ResidueTest, ParameterizedConstructor) {
    EXPECT_EQ("ALA", residue.getName());
    EXPECT_EQ(1, residue.getNumber());
    EXPECT_EQ('A', residue.getChainID());
    EXPECT_EQ(-1, residue.getSSEIndex());
    EXPECT_EQ(5, residue.getAtoms().size());
}

// Atom management tests
TEST_F(ResidueTest, AddAtom) {
    size_t initialCount = residue.getAtoms().size();
    residue.addAtom(Atom(Vec3(6.0, 7.0, 8.0), "CG", "ALA", 1, 'A'));
    EXPECT_EQ(initialCount + 1, residue.getAtoms().size());
}

TEST_F(ResidueTest, GetAtoms) {
    const std::vector<Atom>& atoms = residue.getAtoms();
    EXPECT_EQ(5, atoms.size());
    EXPECT_EQ("CA", atoms[0].getAtomName());
    EXPECT_EQ("CB", atoms[1].getAtomName());
    EXPECT_EQ("N", atoms[2].getAtomName());
    EXPECT_EQ("C", atoms[3].getAtomName());
    EXPECT_EQ("O", atoms[4].getAtomName());
}

TEST_F(ResidueTest, GetAtomByName) {
    // Test getting an existing atom
    Atom* ca = residue.getAtomByName("CA");
    ASSERT_NE(nullptr, ca);
    EXPECT_EQ("CA", ca->getAtomName());
    EXPECT_EQ(Vec3(1.0, 2.0, 3.0), ca->getPosition());
    
    // Test getting a non-existing atom
    Atom* nonExistent = residue.getAtomByName("ZZ");
    EXPECT_EQ(nullptr, nonExistent);
    
    // Test const version
    const Residue& constResidue = residue;
    const Atom* constCa = constResidue.getAtomByName("CA");
    ASSERT_NE(nullptr, constCa);
    EXPECT_EQ("CA", constCa->getAtomName());
}

TEST_F(ResidueTest, GetAlphaCarbon) {
    // Test getting the alpha carbon
    Atom* ca = residue.getAlphaCarbon();
    ASSERT_NE(nullptr, ca);
    EXPECT_EQ("CA", ca->getAtomName());
    
    // Test when there's no alpha carbon
    Residue noAlphaResidue("GLY", 2, 'A');
    noAlphaResidue.addAtom(Atom(Vec3(1.0, 2.0, 3.0), "CB", "GLY", 2, 'A'));
    EXPECT_EQ(nullptr, noAlphaResidue.getAlphaCarbon());
    
    // Test const version
    const Residue& constResidue = residue;
    const Atom* constCa = constResidue.getAlphaCarbon();
    ASSERT_NE(nullptr, constCa);
    EXPECT_EQ("CA", constCa->getAtomName());
}

// Other getter tests
TEST_F(ResidueTest, GetName) {
    EXPECT_EQ("ALA", residue.getName());
}

TEST_F(ResidueTest, GetCode) {
    // Check some common residue codes
    EXPECT_EQ('A', residue.getCode());
    
    Residue gly("GLY", 2, 'A');
    EXPECT_EQ('G', gly.getCode());
    
    Residue cys("CYS", 3, 'A');
    EXPECT_EQ('C', cys.getCode());
    
    Residue unk("UNK", 4, 'A');
    EXPECT_EQ('X', unk.getCode());
}

TEST_F(ResidueTest, GetNumber) {
    EXPECT_EQ(1, residue.getNumber());
}

TEST_F(ResidueTest, GetChainID) {
    EXPECT_EQ('A', residue.getChainID());
}

// SSE-related tests
TEST_F(ResidueTest, SetGetSSEIndex) {
    EXPECT_EQ(-1, residue.getSSEIndex());  // Initial value
    
    int newIndex = 3;
    residue.setSSEIndex(newIndex);
    EXPECT_EQ(newIndex, residue.getSSEIndex());
    
    // Check that all atoms got the SSE index
    for (const Atom& atom : residue.getAtoms()) {
        EXPECT_EQ(newIndex, atom.getSSEIndex());
    }
}

// Center calculation test
TEST_F(ResidueTest, GetCenter) {
    // Center should be the average of all atom positions
    // (1,2,3) + (2,3,4) + (3,4,5) + (4,5,6) + (5,6,7) / 5 = (3,4,5)
    Vec3 expectedCenter(3.0, 4.0, 5.0);
    Vec3 center = residue.getCenter();
    
    EXPECT_DOUBLE_EQ(expectedCenter.x, center.x);
    EXPECT_DOUBLE_EQ(expectedCenter.y, center.y);
    EXPECT_DOUBLE_EQ(expectedCenter.z, center.z);
}