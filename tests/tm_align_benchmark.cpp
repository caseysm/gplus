#include <gtest/gtest.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include <cstdlib> // for system()

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

class TMAlignBenchmarkTest : public ::testing::Test {
protected:
    Config config;
    GPlusConfig gplusConfig;
    PDBParser parser;
    
    // Paths for data
    std::string dataDir;
    std::string benchmarkDir;
    std::string resultsDir;
    
    // Test PDB files
    std::vector<std::string> pdbPaths;
    std::string pdb1Path;
    std::string pdb2Path;
    
    void SetUp() override {
        // Configure default settings
        config.setCoreDistance(4.0);
        config.setResidueDistance(8.0);
        
        gplusConfig.setResultCount(5);
        gplusConfig.setEvaluationDepth(10);
        
        // Set up directories
        dataDir = "/home/casey/Desktop/repos/gplus/tests/data/";
        benchmarkDir = dataDir + "benchmark/";
        resultsDir = dataDir + "results/";
        
        // Create directories if they don't exist
        if (!test::utils::directoryExists(dataDir)) {
            test::utils::createDirectory(dataDir);
        }
        if (!test::utils::directoryExists(benchmarkDir)) {
            test::utils::createDirectory(benchmarkDir);
        }
        if (!test::utils::directoryExists(resultsDir)) {
            test::utils::createDirectory(resultsDir);
        }
        
        // Use existing PDB files from the project root
        pdb1Path = "/home/casey/Desktop/repos/gplus/d1gkub1.pdb";
        pdb2Path = "/home/casey/Desktop/repos/gplus/d2uaga1.pdb";
        
        // Make sure test PDB files exist
        ASSERT_TRUE(File::exists(pdb1Path)) << "Test PDB file not found: " << pdb1Path;
        ASSERT_TRUE(File::exists(pdb2Path)) << "Test PDB file not found: " << pdb2Path;
        
        // Configure logger to show warnings only
        Logger::getInstance().setLevel(Logger::WARNING);
    }
    
    // Helper to execute the original implementation and get its output
    bool runOriginalImplementation(const std::string& pdb1Path, const std::string& pdb2Path, 
                                  bool sequential, double coreDistance, double residueDistance,
                                  std::string& output) {
        // Create temporary file for output
        std::string tmpOutputFile = resultsDir + "original_output.txt";
        
        // Copy the PDB files to the original implementation build directory
        std::string pdb1FileName = File::getFileName(pdb1Path);
        std::string pdb2FileName = File::getFileName(pdb2Path);
        std::string copyCmd1 = "cp " + pdb1Path + " /home/casey/Desktop/repos/gplus/original_code/build/ && chmod +x /home/casey/Desktop/repos/gplus/original_code/build/" + pdb1FileName;
        std::string copyCmd2 = "cp " + pdb2Path + " /home/casey/Desktop/repos/gplus/original_code/build/ && chmod +x /home/casey/Desktop/repos/gplus/original_code/build/" + pdb2FileName;
        
        std::cout << "Copying PDB files with commands: " << std::endl;
        std::cout << "  " << copyCmd1 << std::endl;
        std::cout << "  " << copyCmd2 << std::endl;
        
        if (system(copyCmd1.c_str()) != 0 || system(copyCmd2.c_str()) != 0) {
            std::cerr << "Failed to copy PDB files to original implementation build directory" << std::endl;
            return false;
        }
        
        // Verify the files exist
        std::string verifyCmd = "ls -la /home/casey/Desktop/repos/gplus/original_code/build/";
        std::cout << "Verifying files with command: " << verifyCmd << std::endl;
        system(verifyCmd.c_str());
        
        // Construct command
        std::string cmd = "cd /home/casey/Desktop/repos/gplus/original_code/build && ./gplus";
        if (sequential) {
            cmd += " -s";
        }
        cmd += " -c " + std::to_string(static_cast<int>(coreDistance));
        cmd += " -r " + std::to_string(residueDistance);
        cmd += " " + pdb1FileName + " " + pdb2FileName;
        // Create a temporary output file in the original_code/build directory
        std::string localOutputFile = "/home/casey/Desktop/repos/gplus/original_code/build/test_output.txt";
        cmd += " > " + localOutputFile + " 2>&1";
        
        // Execute command
        std::cout << "Running original implementation: " << cmd << std::endl;
        int result = system(cmd.c_str());
        
        if (result != 0) {
            std::cerr << "Failed to run original implementation" << std::endl;
            return false;
        }
        
        // Read the output file from the local directory
        std::ifstream outputFile(localOutputFile);
        if (!outputFile.is_open()) {
            std::cerr << "Failed to open output file: " << localOutputFile << std::endl;
            return false;
        }
        
        output = std::string((std::istreambuf_iterator<char>(outputFile)),
                              std::istreambuf_iterator<char>());
        
        outputFile.close();
        
        // Copy the output to the results directory for reference
        std::string copyOutputCmd = "cp " + localOutputFile + " " + tmpOutputFile;
        system(copyOutputCmd.c_str());
        
        return true;
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
    
    // Parse alignment results from the original implementation output
    bool parseOriginalOutput(const std::string& output, int& alignedCount, double& rmsd, double& score) {
        // Find the relevant line with alignment results
        std::istringstream stream(output);
        std::string line;
        
        while (std::getline(stream, line)) {
            // Look for the line with "Aligned length= NN, RMSD= NN.NN, GP-Score=NN.NNNNN"
            if (line.find("Aligned length=") != std::string::npos && line.find("RMSD=") != std::string::npos) {
                // Parse aligned count
                size_t alignedPos = line.find("Aligned length=");
                size_t rmsdPos = line.find("RMSD=");
                size_t scorePos = line.find("GP-Score=");
                
                if (alignedPos != std::string::npos && rmsdPos != std::string::npos && scorePos != std::string::npos) {
                    try {
                        // Extract the aligned length
                        std::string alignedStr = line.substr(alignedPos + 15, rmsdPos - (alignedPos + 15));
                        alignedStr = alignedStr.substr(0, alignedStr.find(','));
                        alignedCount = std::stoi(alignedStr);
                        
                        // Extract the RMSD
                        std::string rmsdStr = line.substr(rmsdPos + 5, scorePos - (rmsdPos + 5));
                        rmsdStr = rmsdStr.substr(0, rmsdStr.find(','));
                        rmsd = std::stod(rmsdStr);
                        
                        // Extract the score
                        std::string scoreStr = line.substr(scorePos + 9);
                        score = std::stod(scoreStr);
                        
                        return true; // Successfully parsed all values
                    } catch (const std::exception& e) {
                        std::cerr << "Error parsing original output: " << e.what() << std::endl;
                        return false;
                    }
                }
            }
        }
        
        std::cerr << "Could not find alignment results in output:\n" << output << std::endl;
        return false;
    }
};

// Test the modular implementation with TM-align PDB files
TEST_F(TMAlignBenchmarkTest, CompareModularResults) {
    // Setup results file
    std::string resultsFile = resultsDir + "comparison_results.csv";
    std::ofstream outFile(resultsFile);
    outFile << "PDB1,PDB2,Mode,AlignedCount,RMSD,Score,Time" << std::endl;
    
    // Run comparison for sequential and non-sequential alignment
    for (bool sequential : {true, false}) {
        std::string mode = sequential ? "Sequential" : "NonSequential";
        std::cout << "\n--- " << mode << " Alignment ---" << std::endl;
        
        // Variables to store results
        int alignedCount = 0;
        double rmsd = 0.0;
        double score = 0.0;
        double time = 0.0;
        
        // Run modular implementation
        time = measureExecutionTime([&]() {
            // Load molecules
            Molecule molecule1 = parser.parsePDB(pdb1Path);
            Molecule molecule2 = parser.parsePDB(pdb2Path);
            
            ASSERT_GT(molecule1.getResidueCount(), 0);
            ASSERT_GT(molecule2.getResidueCount(), 0);
            
            // Create alignment
            std::unique_ptr<Alignment> alignment;
            if (sequential) {
                alignment = Alignment::createSequentialAlignment(config, gplusConfig);
            } else {
                alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
            }
            
            alignment->setMolecules(molecule1, molecule2);
            AlignmentResult result = alignment->align();
            
            // Store results
            alignedCount = result.getAlignedCount();
            rmsd = result.getRMSD();
            score = result.getScore();
        }, 5);
        
        std::cout << "Modular implementation results:" << std::endl;
        std::cout << "- Aligned residues: " << alignedCount << std::endl;
        std::cout << "- RMSD: " << rmsd << std::endl;
        std::cout << "- Score: " << score << std::endl;
        std::cout << "- Time: " << time << " s" << std::endl;
        
        // Write results to CSV
        outFile << "d1gkub1.pdb" << "," << "d2uaga1.pdb" << "," << mode << ","
                << alignedCount << "," << rmsd << "," << score << "," << time << std::endl;
    }
    
    outFile.close();
    std::cout << "\nComparison results saved to: " << resultsFile << std::endl;
    
    // Check for expected results
    // In sequential mode, we expect to have a high number of aligned residues
    // In non-sequential mode, we expect a better RMSD
}

// Run performance benchmark with different parameters
TEST_F(TMAlignBenchmarkTest, PerformanceBenchmark) {
    // Setup results file
    std::string resultsFile = resultsDir + "performance_benchmark.csv";
    std::ofstream outFile(resultsFile);
    outFile << "PDB1,PDB2,Mode,CoreDistance,ResidueDistance,"
            << "AlignedCount,RMSD,Score,Time" << std::endl;
    
    // Create a grid of parameters to test
    std::vector<double> coreDistances = {3.0, 4.0, 5.0};
    std::vector<double> residueDistances = {6.0, 8.0, 10.0};
    
    // Load molecules once outside the loop
    Molecule molecule1 = parser.parsePDB(pdb1Path);
    Molecule molecule2 = parser.parsePDB(pdb2Path);
    
    ASSERT_GT(molecule1.getResidueCount(), 0);
    ASSERT_GT(molecule2.getResidueCount(), 0);
    
    // Run tests for both sequential and non-sequential modes
    for (bool sequential : {true, false}) {
        std::string mode = sequential ? "Sequential" : "NonSequential";
        std::cout << "\n--- " << mode << " Alignment ---" << std::endl;
        
        // Test each parameter combination
        for (double coreDistance : coreDistances) {
            for (double residueDistance : residueDistances) {
                std::cout << "\nParameters: Core Distance = " << coreDistance
                         << ", Residue Distance = " << residueDistance << std::endl;
                
                // Set parameters
                config.setCoreDistance(coreDistance);
                config.setResidueDistance(residueDistance);
                
                // Create alignment
                std::unique_ptr<Alignment> alignment;
                if (sequential) {
                    alignment = Alignment::createSequentialAlignment(config, gplusConfig);
                } else {
                    alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
                }
                
                alignment->setMolecules(molecule1, molecule2);
                
                // Timing and results
                double alignTime = 0.0;
                int alignedCount = 0;
                double rmsd = 0.0;
                double score = 0.0;
                
                alignTime = measureExecutionTime([&]() {
                    AlignmentResult result = alignment->align();
                    alignedCount = result.getAlignedCount();
                    rmsd = result.getRMSD();
                    score = result.getScore();
                }, 5); // Increase to 5 repetitions for better timing accuracy
                
                std::cout << "Results:" << std::endl;
                std::cout << "- Aligned residues: " << alignedCount << std::endl;
                std::cout << "- RMSD: " << rmsd << std::endl;
                std::cout << "- Score: " << score << std::endl;
                std::cout << "- Time: " << alignTime << " s" << std::endl;
                
                // Write results to CSV
                outFile << "d1gkub1.pdb" << "," << "d2uaga1.pdb" << "," << mode << ","
                        << coreDistance << "," << residueDistance << ","
                        << alignedCount << "," << rmsd << "," << score << ","
                        << alignTime << std::endl;
            }
        }
    }
    
    outFile.close();
    std::cout << "\nPerformance benchmark results saved to: " << resultsFile << std::endl;
}

// Benchmark multiple executions to analyze performance stability
TEST_F(TMAlignBenchmarkTest, PerformanceStabilityBenchmark) {
    // Setup results file
    std::string resultsFile = resultsDir + "stability_benchmark.csv";
    std::ofstream outFile(resultsFile);
    outFile << "Mode,Run,Aligned,RMSD,Score,Time" << std::endl;
    
    const int numRuns = 10; // Number of times to run each implementation
    
    // Run comparison for sequential and non-sequential alignment
    for (bool sequential : {true, false}) {
        std::string mode = sequential ? "Sequential" : "NonSequential";
        std::cout << "\n--- " << mode << " Alignment Stability Test ---" << std::endl;
        
        // Test modular implementation multiple times
        std::cout << "\nTesting modular implementation stability..." << std::endl;
        
        // Load molecules once
        Molecule molecule1 = parser.parsePDB(pdb1Path);
        Molecule molecule2 = parser.parsePDB(pdb2Path);
        
        ASSERT_GT(molecule1.getResidueCount(), 0);
        ASSERT_GT(molecule2.getResidueCount(), 0);
        
        for (int run = 0; run < numRuns; run++) {
            std::cout << "Run " << (run + 1) << "/" << numRuns << std::endl;
            
            // Create alignment
            std::unique_ptr<Alignment> alignment;
            if (sequential) {
                alignment = Alignment::createSequentialAlignment(config, gplusConfig);
            } else {
                alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
            }
            
            alignment->setMolecules(molecule1, molecule2);
            
            // Timing and results
            int alignedCount = 0;
            double rmsd = 0.0, score = 0.0, time = 0.0;
            
            time = measureExecutionTime([&]() {
                AlignmentResult result = alignment->align();
                alignedCount = result.getAlignedCount();
                rmsd = result.getRMSD();
                score = result.getScore();
            }, 1);
            
            std::cout << "- Aligned residues: " << alignedCount << std::endl;
            std::cout << "- RMSD: " << rmsd << std::endl;
            std::cout << "- Score: " << score << std::endl;
            std::cout << "- Time: " << time << " s" << std::endl;
            
            // Write results
            outFile << mode << "," << (run + 1) << "," 
                    << alignedCount << "," << rmsd << "," << score << "," << time << std::endl;
        }
    }
    
    outFile.close();
    std::cout << "\nStability benchmark results saved to: " << resultsFile << std::endl;
}