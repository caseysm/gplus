#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <limits>
#include <fstream>

#include "config/config.h"
#include "config/gplus_config.h"
#include "core/molecule.h"
#include "core/atom.h"
#include "core/residue.h"
#include "core/sse.h"
#include "math/vec3.h"
#include "math/rotation.h"
#include "algorithm/alignment.h"
#include "utils/pdb_parser.h"
#include "utils/file.h"
#include "utils/logger.h"
#include "file_utils.h"

using namespace gangsta::config;
using namespace gangsta::core;
using namespace gangsta::algorithm;
using namespace gangsta::math;
using namespace gangsta::utils;

class EdgeCaseTest : public ::testing::Test {
protected:
    Config config;
    GPlusConfig gplusConfig;
    PDBParser parser;
    
    // Output directory for temporary files
    std::string outputDir;
    
    void SetUp() override {
        // Configure default settings
        config.setCoreDistance(4.0);
        config.setResidueDistance(8.0);
        
        gplusConfig.setResultCount(3);
        gplusConfig.setEvaluationDepth(10);
        
        // Set up output directory
        outputDir = "./edge_case_tests/";
        if (!test::utils::directoryExists(outputDir)) {
            test::utils::createDirectory(outputDir);
        }
        
        // Configure logger for testing
        Logger::getInstance().setLevel(Logger::WARNING);
    }
    
    void TearDown() override {
        // Clean up output directory if needed
        // Note: We leave it for inspection
    }
    
    // Helper to create a minimal PDB file
    std::string createMinimalPDB(const std::string& filename, int numResidues = 1) {
        std::string filePath = outputDir + filename;
        std::ofstream file(filePath);
        
        for (int i = 0; i < numResidues; i++) {
            // Write a minimal ATOM record for an alpha carbon
            file << "ATOM  " << std::setw(5) << (i + 1) << " CA  ALA A" << std::setw(4) << (i + 1) 
                 << "    " << std::fixed << std::setprecision(3) 
                 << std::setw(8) << (i * 3.8) << std::setw(8) << 0.000 << std::setw(8) << 0.000 
                 << "  1.00  0.00           C  " << std::endl;
        }
        
        file << "END" << std::endl;
        file.close();
        
        return filePath;
    }
    
    // Helper to create an extremely large PDB file
    std::string createLargePDB(const std::string& filename, int numResidues) {
        std::string filePath = outputDir + filename;
        std::ofstream file(filePath);
        
        for (int i = 0; i < numResidues; i++) {
            // Write a minimal ATOM record for an alpha carbon
            file << "ATOM  " << std::setw(5) << (i + 1) << " CA  ALA A" << std::setw(4) << ((i % 9999) + 1)
                 << "    " << std::fixed << std::setprecision(3) 
                 << std::setw(8) << (i * 3.8) << std::setw(8) << 0.000 << std::setw(8) << 0.000 
                 << "  1.00  0.00           C  " << std::endl;
        }
        
        file << "END" << std::endl;
        file.close();
        
        return filePath;
    }
};

// Test with empty PDB file
TEST_F(EdgeCaseTest, EmptyPDBFile) {
    // Create an empty PDB file
    std::string emptyFilePath = outputDir + "empty.pdb";
    std::ofstream emptyFile(emptyFilePath);
    emptyFile << "END" << std::endl;
    emptyFile.close();
    
    // Attempt to parse the empty file
    Molecule emptyMolecule = parser.parsePDB(emptyFilePath);
    
    // Verify that parsing an empty file results in an empty molecule
    EXPECT_EQ(0, emptyMolecule.getResidueCount());
    EXPECT_EQ(0, emptyMolecule.getSSECount());
    EXPECT_FALSE(emptyMolecule.getName().empty()); // Name should still be set from filename
    
    // Try to align with an empty molecule
    Molecule normalMolecule = parser.parsePDB(createMinimalPDB("normal.pdb", 5));
    
    // Sequential alignment with empty molecule
    std::unique_ptr<Alignment> seqAlignment = Alignment::createSequentialAlignment(config, gplusConfig);
    seqAlignment->setMolecules(normalMolecule, emptyMolecule);
    
    AlignmentResult seqResult = seqAlignment->align();
    
    // Verify that alignment handles empty molecule gracefully
    EXPECT_EQ(0, seqResult.getAlignedCount());
}

// Test with single residue molecules
TEST_F(EdgeCaseTest, SingleResidueMolecule) {
    // Create a PDB file with a single residue
    std::string singleResFilePath = createMinimalPDB("single_residue.pdb", 1);
    
    // Parse the single residue file
    Molecule singleResMolecule = parser.parsePDB(singleResFilePath);
    
    // Verify proper parsing
    EXPECT_EQ(1, singleResMolecule.getResidueCount());
    
    // Try to align with another single residue molecule
    Molecule anotherSingleResMolecule = parser.parsePDB(createMinimalPDB("another_single.pdb", 1));
    
    // Sequential alignment with two single-residue molecules
    std::unique_ptr<Alignment> seqAlignment = Alignment::createSequentialAlignment(config, gplusConfig);
    seqAlignment->setMolecules(singleResMolecule, anotherSingleResMolecule);
    
    AlignmentResult seqResult = seqAlignment->align();
    
    // Verify that alignment works with minimal molecules
    EXPECT_EQ(1, seqResult.getAlignedCount());
    
    // Try non-sequential alignment
    std::unique_ptr<Alignment> nonSeqAlignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
    nonSeqAlignment->setMolecules(singleResMolecule, anotherSingleResMolecule);
    
    AlignmentResult nonSeqResult = nonSeqAlignment->align();
    
    // Verify that alignment works with minimal molecules
    EXPECT_EQ(1, nonSeqResult.getAlignedCount());
}

// Test with very distant or unusual structures
TEST_F(EdgeCaseTest, VeryDistantStructures) {
    // Create two molecules with vastly different structures
    Molecule normalMolecule = parser.parsePDB(createMinimalPDB("normal_struct.pdb", 10));
    
    // Create a molecule with unusual coordinates (very distant from origin)
    std::string distantFilePath = outputDir + "distant.pdb";
    std::ofstream distantFile(distantFilePath);
    
    for (int i = 0; i < 10; i++) {
        // Write atoms with large coordinate values
        distantFile << "ATOM  " << std::setw(5) << (i + 1) << " CA  ALA A" << std::setw(4) << (i + 1) 
                    << "    " << std::fixed << std::setprecision(3) 
                    << std::setw(8) << (1000.0 + i) << std::setw(8) << (1000.0 + i) << std::setw(8) << (1000.0 + i) 
                    << "  1.00  0.00           C  " << std::endl;
    }
    
    distantFile << "END" << std::endl;
    distantFile.close();
    
    Molecule distantMolecule = parser.parsePDB(distantFilePath);
    
    // Try to align these very different structures
    std::unique_ptr<Alignment> alignment = Alignment::createSequentialAlignment(config, gplusConfig);
    alignment->setMolecules(normalMolecule, distantMolecule);
    
    AlignmentResult result = alignment->align();
    
    // Verify that alignment handles distant structures
    // Due to the nature of the problem, we can't assert specific values,
    // but we can verify that the alignment completes without errors
    EXPECT_GE(result.getAlignedCount(), 0);
    EXPECT_TRUE(std::isfinite(result.getRMSD()));
    EXPECT_TRUE(std::isfinite(result.getScore()));
}

// Test with extreme parameter values
TEST_F(EdgeCaseTest, ExtremeParameters) {
    // Create normal molecules
    Molecule molecule1 = parser.parsePDB(createMinimalPDB("normal1.pdb", 5));
    Molecule molecule2 = parser.parsePDB(createMinimalPDB("normal2.pdb", 5));
    
    // Test with very small core distance
    Config tinyConfig;
    tinyConfig.setCoreDistance(0.001);
    tinyConfig.setResidueDistance(0.001);
    
    std::unique_ptr<Alignment> tinyParamAlignment = Alignment::createSequentialAlignment(tinyConfig, gplusConfig);
    tinyParamAlignment->setMolecules(molecule1, molecule2);
    
    AlignmentResult tinyResult = tinyParamAlignment->align();
    
    // Test with very large core distance
    Config hugeConfig;
    hugeConfig.setCoreDistance(1000.0);
    hugeConfig.setResidueDistance(2000.0);
    
    std::unique_ptr<Alignment> hugeParamAlignment = Alignment::createSequentialAlignment(hugeConfig, gplusConfig);
    hugeParamAlignment->setMolecules(molecule1, molecule2);
    
    AlignmentResult hugeResult = hugeParamAlignment->align();
    
    // Test extreme GPlusConfig
    GPlusConfig extremeGPlusConfig;
    extremeGPlusConfig.setResultCount(100);
    extremeGPlusConfig.setEvaluationDepth(100);
    
    std::unique_ptr<Alignment> extremeAlignment = Alignment::createNonSequentialAlignment(config, extremeGPlusConfig);
    extremeAlignment->setMolecules(molecule1, molecule2);
    
    AlignmentResult extremeResult = extremeAlignment->align();
    
    // Verify that alignment handles extreme parameters
    // Again, we're primarily testing that it doesn't crash or produce NaN results
    EXPECT_TRUE(std::isfinite(tinyResult.getRMSD()));
    EXPECT_TRUE(std::isfinite(hugeResult.getRMSD()));
    EXPECT_TRUE(std::isfinite(extremeResult.getRMSD()));
}

// Test with molecules that have no valid alpha carbons
TEST_F(EdgeCaseTest, NoAlphaCarbons) {
    // Create a PDB file with no alpha carbons
    std::string noCAFilePath = outputDir + "no_ca.pdb";
    std::ofstream noCAFile(noCAFilePath);
    
    // Write atoms that are not alpha carbons
    noCAFile << "ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  " << std::endl;
    noCAFile << "ATOM      2  CB  ALA A   2       3.800   0.000   0.000  1.00  0.00           C  " << std::endl;
    noCAFile << "ATOM      3  CB  ALA A   3       7.600   0.000   0.000  1.00  0.00           C  " << std::endl;
    noCAFile << "END" << std::endl;
    noCAFile.close();
    
    // Parse the file without alpha carbons
    Molecule noCAMolecule = parser.parsePDB(noCAFilePath);
    
    // Verify that residues are parsed
    EXPECT_GT(noCAMolecule.getResidueCount(), 0);
    
    // Try to align with a normal molecule
    Molecule normalMolecule = parser.parsePDB(createMinimalPDB("normal_for_ca.pdb", 5));
    
    std::unique_ptr<Alignment> alignment = Alignment::createSequentialAlignment(config, gplusConfig);
    alignment->setMolecules(normalMolecule, noCAMolecule);
    
    // This should not crash, even if there are no alpha carbons in one molecule
    AlignmentResult result = alignment->align();
    
    // Expect no aligned residues since alignment requires alpha carbons
    EXPECT_EQ(0, result.getAlignedCount());
}

// Test with very large molecules
TEST_F(EdgeCaseTest, VeryLargeMolecules) {
    // Skip this test for regular runs since it creates very large files
    // Only run when specifically needed
    if (true) {
        GTEST_SKIP() << "Skipping large molecule test to avoid creating large files";
        return;
    }
    
    // Create a large molecule with 1000 residues
    std::string largeFilePath = createLargePDB("large.pdb", 1000);
    
    // Parse the large file
    Molecule largeMolecule = parser.parsePDB(largeFilePath);
    
    // Verify proper parsing
    EXPECT_EQ(1000, largeMolecule.getResidueCount());
    
    // Try to align with a smaller molecule
    Molecule smallMolecule = parser.parsePDB(createMinimalPDB("small_for_large.pdb", 10));
    
    // Sequential alignment with one large molecule
    std::unique_ptr<Alignment> seqAlignment = Alignment::createSequentialAlignment(config, gplusConfig);
    seqAlignment->setMolecules(smallMolecule, largeMolecule);
    
    // This should handle the size discrepancy gracefully
    AlignmentResult seqResult = seqAlignment->align();
    
    // Verify some alignment was found (likely matches beginning of large molecule)
    EXPECT_GT(seqResult.getAlignedCount(), 0);
    EXPECT_LE(seqResult.getAlignedCount(), smallMolecule.getResidueCount());
}

// Test with extreme vectors, rotations and transformations
TEST_F(EdgeCaseTest, ExtremeVectorsAndRotations) {
    // Test vectors with extreme values
    Vec3 extremeVec(1e10, -1e10, 1e10);
    Vec3 tinyVec(1e-10, 1e-10, 1e-10);
    
    // Normalize should handle these without producing NaN
    Vec3 normalizedExtreme = extremeVec.normalize();
    Vec3 normalizedTiny = tinyVec.normalize();
    
    EXPECT_NEAR(1.0, normalizedExtreme.length(), 1e-10);
    EXPECT_NEAR(1.0, normalizedTiny.length(), 1e-10);
    
    // Test rotation with very small angle
    Rotation tinyRotation(Vec3(0, 0, 1), 1e-10);
    Vec3 rotatedVec = tinyRotation.rotate(Vec3(1, 0, 0));
    
    // Should be very close to original vector
    EXPECT_NEAR(1.0, rotatedVec.x, 1e-9);
    EXPECT_NEAR(0.0, rotatedVec.y, 1e-9);
    EXPECT_NEAR(0.0, rotatedVec.z, 1e-9);
    
    // Test rotation with NaN values (should not crash)
    Vec3 nanAxis(std::numeric_limits<double>::quiet_NaN(), 0, 0);
    
    // This should detect the NaN and use a default axis
    Rotation nanRotation(nanAxis, 1.0);
    
    // The resulting quaternion shouldn't have NaN components
    Vec4 quat = nanRotation.getQuaternion();
    EXPECT_FALSE(std::isnan(quat.x));
    EXPECT_FALSE(std::isnan(quat.y));
    EXPECT_FALSE(std::isnan(quat.z));
    EXPECT_FALSE(std::isnan(quat.w));
}