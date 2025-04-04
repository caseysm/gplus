#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <chrono>
#include <memory>

#include "config/config.h"
#include "config/gplus_config.h"
#include "core/molecule.h"
#include "algorithm/alignment.h"
#include "utils/pdb_parser.h"
#include "utils/file.h"
#include "utils/logger.h"
#include "utils/pymol_generator.h"
#include "file_utils.h"

using namespace gangsta::config;
using namespace gangsta::core;
using namespace gangsta::algorithm;
using namespace gangsta::utils;

class IntegrationTest : public ::testing::Test {
protected:
    Config config;
    GPlusConfig gplusConfig;
    PDBParser parser;
    
    // Path to test PDB files
    std::string pdbDir;
    std::string pdb1Path;
    std::string pdb2Path;
    
    // Output paths
    std::string outputDir;
    std::string outputPdbPath;
    std::string outputPyMolPath;
    
    void SetUp() override {
        // Configure default settings
        config.setCoreDistance(4.0);
        config.setResidueDistance(8.0);
        
        gplusConfig.setResultCount(3);
        gplusConfig.setEvaluationDepth(10);
        
        // Set up file paths for test PDBs
        pdbDir = "../";  // Look for PDBs in the main project directory
        pdb1Path = pdbDir + "d1gkub1.pdb";
        pdb2Path = pdbDir + "d2uaga1.pdb";
        
        // Make sure test PDB files exist
        ASSERT_TRUE(File::exists(pdb1Path)) << "Test PDB file not found: " << pdb1Path;
        ASSERT_TRUE(File::exists(pdb2Path)) << "Test PDB file not found: " << pdb2Path;
        
        // Set up output directory
        outputDir = "./test_output/";
        if (!test::utils::directoryExists(outputDir)) {
            test::utils::createDirectory(outputDir);
        }
        
        outputPdbPath = outputDir + "test_alignment.pdb";
        outputPyMolPath = outputDir + "test_visualization.pml";
        
        // Configure logger to show info messages
        Logger::getInstance().setLevel(Logger::INFO);
    }
    
    void TearDown() override {
        // Clean up output files after tests
        if (File::exists(outputPdbPath)) {
            test::utils::removeFile(outputPdbPath);
        }
        if (File::exists(outputPyMolPath)) {
            test::utils::removeFile(outputPyMolPath);
        }
    }
};

// Full end-to-end sequential alignment test
TEST_F(IntegrationTest, SequentialAlignmentEndToEnd) {
    // Load molecules
    Molecule molecule1 = parser.parsePDB(pdb1Path);
    Molecule molecule2 = parser.parsePDB(pdb2Path);
    
    ASSERT_GT(molecule1.getResidueCount(), 0) << "Failed to load residues from " << pdb1Path;
    ASSERT_GT(molecule2.getResidueCount(), 0) << "Failed to load residues from " << pdb2Path;
    
    // Perform sequential alignment
    std::unique_ptr<Alignment> alignment = Alignment::createSequentialAlignment(config, gplusConfig);
    alignment->setMolecules(molecule1, molecule2);
    
    AlignmentResult result = alignment->align();
    
    // Verify alignment succeeded
    EXPECT_GT(result.getAlignedCount(), 0) << "Alignment failed - no aligned residues";
    EXPECT_GT(result.getScore(), 0.0) << "Alignment failed - score is zero or negative";
    
    // Verify alignment result has valid rotation and translation
    EXPECT_NEAR(result.getRotation().getQuaternion().length(), 1.0, 0.001);
    
    // Generate alignment visualization
    PyMolGenerator pymolGen;
    
    // Apply alignment to get aligned molecule1
    Molecule alignedMolecule1 = Alignment::applyAlignment(molecule1, result);
    
    // Set up PyMolGenerator
    pymolGen.setMolecules(molecule1, molecule2, alignedMolecule1);
    pymolGen.setAlignmentResult(result);
    
    // Generate script
    pymolGen.generateScript(outputPyMolPath, outputPdbPath);
    
    // Check that output files were created
    EXPECT_TRUE(File::exists(outputPyMolPath)) << "Failed to create PyMOL script";
    
    // Print alignment summary
    std::cout << "Sequential Alignment Result:" << std::endl;
    std::cout << result.getSummary() << std::endl;
}

// Full end-to-end non-sequential alignment test
TEST_F(IntegrationTest, NonSequentialAlignmentEndToEnd) {
    // Load molecules
    Molecule molecule1 = parser.parsePDB(pdb1Path);
    Molecule molecule2 = parser.parsePDB(pdb2Path);
    
    ASSERT_GT(molecule1.getResidueCount(), 0) << "Failed to load residues from " << pdb1Path;
    ASSERT_GT(molecule2.getResidueCount(), 0) << "Failed to load residues from " << pdb2Path;
    
    // Perform non-sequential alignment
    std::unique_ptr<Alignment> alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
    alignment->setMolecules(molecule1, molecule2);
    
    AlignmentResult result = alignment->align();
    
    // Verify alignment succeeded
    EXPECT_GT(result.getAlignedCount(), 0) << "Alignment failed - no aligned residues";
    EXPECT_GT(result.getScore(), 0.0) << "Alignment failed - score is zero or negative";
    
    // Verify alignment result has valid rotation and translation
    EXPECT_NEAR(result.getRotation().getQuaternion().length(), 1.0, 0.001);
    
    // Generate alignment visualization
    PyMolGenerator pymolGen;
    
    // Apply alignment to get aligned molecule1
    Molecule alignedMolecule1 = Alignment::applyAlignment(molecule1, result);
    
    // Set up PyMolGenerator
    pymolGen.setMolecules(molecule1, molecule2, alignedMolecule1);
    pymolGen.setAlignmentResult(result);
    
    // Generate script
    pymolGen.generateScript(outputPyMolPath, outputPdbPath);
    
    // Check that output files were created
    EXPECT_TRUE(File::exists(outputPyMolPath)) << "Failed to create PyMOL script";
    
    // Print alignment summary
    std::cout << "Non-Sequential Alignment Result:" << std::endl;
    std::cout << result.getSummary() << std::endl;
}

// Test to verify alignment consistency across runs
TEST_F(IntegrationTest, AlignmentConsistency) {
    // Load molecules
    Molecule molecule1 = parser.parsePDB(pdb1Path);
    Molecule molecule2 = parser.parsePDB(pdb2Path);
    
    // Run alignment multiple times and verify results are consistent
    std::vector<AlignmentResult> results;
    const int numRuns = 3;
    
    for (int i = 0; i < numRuns; i++) {
        std::unique_ptr<Alignment> alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
        alignment->setMolecules(molecule1, molecule2);
        results.push_back(alignment->align());
    }
    
    // Verify all alignments produced the same number of aligned residues
    for (int i = 1; i < numRuns; i++) {
        EXPECT_EQ(results[0].getAlignedCount(), results[i].getAlignedCount())
            << "Run " << i << " produced different number of aligned residues";
        
        // Scores should be very close, if not identical
        EXPECT_NEAR(results[0].getScore(), results[i].getScore(), 0.001)
            << "Run " << i << " produced significantly different score";
        
        // RMSDs should be very close, if not identical
        EXPECT_NEAR(results[0].getRMSD(), results[i].getRMSD(), 0.001)
            << "Run " << i << " produced significantly different RMSD";
    }
}

// File format error handling test
TEST_F(IntegrationTest, FileFormatErrorHandling) {
    // Create a corrupt/invalid PDB file
    std::string invalidPdbPath = outputDir + "invalid.pdb";
    std::ofstream invalidFile(invalidPdbPath);
    invalidFile << "This is not a valid PDB file format" << std::endl;
    invalidFile.close();
    
    // Attempt to parse invalid file
    try {
        Molecule invalidMolecule = parser.parsePDB(invalidPdbPath);
        
        // Even with an invalid file, the parser should create a molecule object
        // but it may have no residues or atoms
        EXPECT_EQ(0, invalidMolecule.getResidueCount()) 
            << "Parser should handle invalid PDB by creating empty molecule";
    }
    catch (const std::exception& e) {
        // If an exception is thrown, ensure it has a descriptive message
        EXPECT_TRUE(std::string(e.what()).size() > 0) 
            << "Exception from invalid PDB should have descriptive message";
    }
    
    // Clean up
    if (File::exists(invalidPdbPath)) {
        test::utils::removeFile(invalidPdbPath);
    }
}