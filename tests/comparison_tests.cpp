#include <gtest/gtest.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include "config/config.h"
#include "config/gplus_config.h"
#include "core/molecule.h"
#include "algorithm/alignment.h"
#include "utils/pdb_parser.h"
#include "utils/file.h"
#include "utils/logger.h"
#include "file_utils.h"

using namespace gangsta::config;
using namespace gangsta::core;
using namespace gangsta::algorithm;
using namespace gangsta::utils;

class ComparisonTest : public ::testing::Test {
protected:
    Config config;
    GPlusConfig gplusConfig;
    PDBParser parser;
    
    // Path to test PDB files
    std::string pdbDir;
    std::string pdb1Path;
    std::string pdb2Path;
    
    // Prepare test data
    void SetUp() override {
        // Configure default settings
        config.setCoreDistance(4.0);
        config.setResidueDistance(8.0);
        
        gplusConfig.setResultCount(5);
        gplusConfig.setEvaluationDepth(10);
        
        // Set up file paths for test PDBs
        pdbDir = "../";  // Look for PDBs in the main project directory
        pdb1Path = pdbDir + "d1gkub1.pdb";
        pdb2Path = pdbDir + "d2uaga1.pdb";
        
        // Make sure test PDB files exist
        ASSERT_TRUE(File::exists(pdb1Path)) << "Test PDB file not found: " << pdb1Path;
        ASSERT_TRUE(File::exists(pdb2Path)) << "Test PDB file not found: " << pdb2Path;
        
        // Configure logger to show warnings only
        Logger::getInstance().setLevel(Logger::WARNING);
    }
    
    // Helper to measure execution time
    template<typename Func>
    double measureExecutionTime(Func&& func, int repetitions = 3) {
        std::vector<double> times;
        
        for (int i = 0; i < repetitions; i++) {
            auto start = std::chrono::high_resolution_clock::now();
            func();
            auto end = std::chrono::high_resolution_clock::now();
            
            double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            times.push_back(duration);
        }
        
        // Calculate mean
        double sum = 0.0;
        for (double time : times) {
            sum += time;
        }
        double mean = sum / times.size();
        
        // Calculate standard deviation
        double sumSquaredDiff = 0.0;
        for (double time : times) {
            double diff = time - mean;
            sumSquaredDiff += diff * diff;
        }
        double stdDev = std::sqrt(sumSquaredDiff / times.size());
        
        std::cout << "Execution time over " << repetitions << " repetitions: " 
                << std::fixed << std::setprecision(3) << mean 
                << " s (±" << stdDev << " s)" << std::endl;
        
        return mean;
    }
    
    // Direct execution of both original and modular implementations
    void runDirectComparison(bool sequential = false) {
        // Load the same PDB files for both implementations
        std::cout << "Running direct comparison" << (sequential ? " (sequential)" : " (non-sequential)") << std::endl;
        
        // Run modular implementation
        std::cout << "Running modular implementation..." << std::endl;
        double moduleTime = measureExecutionTime([&]() {
            // Load molecules
            Molecule molecule1 = parser.parsePDB(pdb1Path);
            Molecule molecule2 = parser.parsePDB(pdb2Path);
            
            ASSERT_GT(molecule1.getResidueCount(), 0) << "Failed to load residues from " << pdb1Path;
            ASSERT_GT(molecule2.getResidueCount(), 0) << "Failed to load residues from " << pdb2Path;
            
            // Create alignment based on sequential flag
            std::unique_ptr<Alignment> alignment;
            if (sequential) {
                alignment = Alignment::createSequentialAlignment(config, gplusConfig);
            } else {
                alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
            }
            
            alignment->setMolecules(molecule1, molecule2);
            AlignmentResult result = alignment->align();
            
            // Check that we got a valid result
            EXPECT_GT(result.getAlignedCount(), 0) << "Modular implementation failed to align residues";
            EXPECT_GT(result.getScore(), 0.0) << "Modular implementation produced zero score";
            
            // Print results
            std::cout << "Modular implementation results:" << std::endl;
            std::cout << "- Aligned residues: " << result.getAlignedCount() << std::endl;
            std::cout << "- RMSD: " << result.getRMSD() << std::endl;
            std::cout << "- Score: " << result.getScore() << std::endl;
        });
        
        // Original implementation will be run manually for comparison
        std::cout << std::endl;
        std::cout << "To compare with original implementation, run:" << std::endl;
        
        std::string originalCommand = "../original_code/build/gplus";
        if (sequential) {
            originalCommand += " -s";
        }
        originalCommand += " -c " + std::to_string(static_cast<int>(config.getCoreDistance()));
        originalCommand += " -r " + std::to_string(config.getResidueDistance());
        originalCommand += " " + pdb1Path + " " + pdb2Path;
        
        std::cout << originalCommand << std::endl;
    }
};

// Test for comparing sequential alignment
TEST_F(ComparisonTest, SequentialAlignmentBasic) {
    runDirectComparison(true);
}

// Test for comparing non-sequential alignment
TEST_F(ComparisonTest, NonSequentialAlignmentBasic) {
    runDirectComparison(false);
}

// Test with different core distances
TEST_F(ComparisonTest, CoreDistanceComparison) {
    std::cout << "Testing with different core distances" << std::endl;
    
    std::vector<double> coreDistances = {3.0, 4.0, 5.0, 6.0};
    
    for (double coreDistance : coreDistances) {
        std::cout << std::endl << "Core Distance = " << coreDistance << std::endl;
        
        // Set the core distance
        config.setCoreDistance(coreDistance);
        
        // Load molecules
        Molecule molecule1 = parser.parsePDB(pdb1Path);
        Molecule molecule2 = parser.parsePDB(pdb2Path);
        
        // Create non-sequential alignment
        std::unique_ptr<Alignment> alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
        alignment->setMolecules(molecule1, molecule2);
        
        // Perform alignment and measure time
        double alignTime = measureExecutionTime([&]() {
            return alignment->align();
        }, 1);
        
        // Get the alignment result outside of timing
        AlignmentResult result = alignment->align();
        
        std::cout << "Results with core distance = " << coreDistance << ":" << std::endl;
        std::cout << "- Aligned residues: " << result.getAlignedCount() << std::endl;
        std::cout << "- RMSD: " << result.getRMSD() << std::endl;
        std::cout << "- Score: " << result.getScore() << std::endl;
        std::cout << "- Time: " << alignTime << " s" << std::endl;
    }
}

// Benchmark of parameter sensitivity
TEST_F(ComparisonTest, ParameterSensitivityBenchmark) {
    std::cout << "Performing parameter sensitivity analysis" << std::endl;
    
    // Create a grid of parameters to test
    std::vector<double> coreDistances = {3.0, 4.0, 6.0};
    std::vector<double> residueDistances = {6.0, 8.0, 10.0};
    
    // Store results for analysis
    struct Result {
        double coreDistance;
        double residueDistance;
        int alignedCount;
        double rmsd;
        double score;
        double time;
    };
    
    std::vector<Result> results;
    
    // Run the benchmark grid
    for (double coreDistance : coreDistances) {
        for (double residueDistance : residueDistances) {
            std::cout << std::endl 
                    << "Testing with core distance = " << coreDistance 
                    << ", residue distance = " << residueDistance << std::endl;
            
            // Set parameters
            config.setCoreDistance(coreDistance);
            config.setResidueDistance(residueDistance);
            
            // Load molecules
            Molecule molecule1 = parser.parsePDB(pdb1Path);
            Molecule molecule2 = parser.parsePDB(pdb2Path);
            
            // Create non-sequential alignment
            std::unique_ptr<Alignment> alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
            alignment->setMolecules(molecule1, molecule2);
            
            // Measure alignment time
            double alignTime = 0.0;
            AlignmentResult result;
            
            alignTime = measureExecutionTime([&]() {
                result = alignment->align();
            }, 1);
            
            // Store result
            Result benchResult = {
                coreDistance,
                residueDistance,
                result.getAlignedCount(),
                result.getRMSD(),
                result.getScore(),
                alignTime
            };
            
            results.push_back(benchResult);
            
            std::cout << "- Aligned residues: " << result.getAlignedCount() << std::endl;
            std::cout << "- RMSD: " << result.getRMSD() << std::endl;
            std::cout << "- Score: " << result.getScore() << std::endl;
            std::cout << "- Time: " << alignTime << " s" << std::endl;
        }
    }
    
    // Print summary table
    std::cout << std::endl << "Parameter Sensitivity Summary:" << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;
    std::cout << "| Core Dist | Res Dist | Aligned | RMSD    | Score   | Time   |" << std::endl;
    std::cout << "|----------|----------|---------|---------|---------|--------|" << std::endl;
    
    for (const Result& r : results) {
        std::cout << "| " << std::setw(8) << r.coreDistance 
                << " | " << std::setw(8) << r.residueDistance 
                << " | " << std::setw(7) << r.alignedCount 
                << " | " << std::fixed << std::setprecision(2) << std::setw(7) << r.rmsd 
                << " | " << std::fixed << std::setprecision(4) << std::setw(7) << r.score 
                << " | " << std::fixed << std::setprecision(3) << std::setw(6) << r.time 
                << " |" << std::endl;
    }
    
    std::cout << "--------------------------------------------------------------" << std::endl;
}