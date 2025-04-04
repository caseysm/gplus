#include <gtest/gtest.h>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "config/config.h"
#include "config/gplus_config.h"
#include "core/molecule.h"
#include "algorithm/alignment.h"
#include "utils/pdb_parser.h"
#include "utils/logger.h"
#include "utils/file.h"
#include "file_utils.h"

using namespace gangsta::config;
using namespace gangsta::core;
using namespace gangsta::algorithm;
using namespace gangsta::utils;

class BenchmarkTest : public ::testing::Test {
protected:
    Config config;
    GPlusConfig gplusConfig;
    PDBParser parser;
    
    // Path to test PDB files
    std::string pdbDir;
    std::string pdb1Path;
    std::string pdb2Path;
    
    // Output directory for benchmark results
    std::string outputDir;
    std::string benchmarkResultsPath;
    
    // Helper for timing
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    
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
        outputDir = "./benchmark_results/";
        if (!test::utils::directoryExists(outputDir)) {
            test::utils::createDirectory(outputDir);
        }
        
        benchmarkResultsPath = outputDir + "benchmark_results.csv";
        
        // Turn off verbose logging during benchmarks
        Logger::getInstance().setLevel(Logger::WARNING);
    }
    
    // Helper to measure execution time
    template<typename Func>
    double measureExecutionTime(Func&& func, int repetitions = 1) {
        std::vector<double> times;
        
        for (int i = 0; i < repetitions; i++) {
            TimePoint start = Clock::now();
            func();
            TimePoint end = Clock::now();
            
            double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0;
            times.push_back(duration);
        }
        
        // Calculate mean
        double sum = 0.0;
        for (double time : times) {
            sum += time;
        }
        double mean = sum / times.size();
        
        // Calculate standard deviation if more than one repetition
        double stdDev = 0.0;
        if (repetitions > 1) {
            double sumSquaredDiff = 0.0;
            for (double time : times) {
                double diff = time - mean;
                sumSquaredDiff += diff * diff;
            }
            stdDev = std::sqrt(sumSquaredDiff / times.size());
        }
        
        // Log detailed results for multiple repetitions
        if (repetitions > 1) {
            std::cout << "Execution time over " << repetitions << " repetitions: " 
                    << std::fixed << std::setprecision(3) << mean 
                    << " s (±" << stdDev << " s)" << std::endl;
        }
        
        return mean;
    }
    
    // Write benchmark results to CSV
    void writeBenchmarkResult(const std::string& name, double time, int dataSize, 
                            const std::string& algorithm, const std::string& notes = "") {
        bool fileExists = File::exists(benchmarkResultsPath);
        
        std::ofstream file(benchmarkResultsPath, std::ios::app);
        
        // Write header if file is new
        if (!fileExists) {
            file << "Timestamp,Benchmark,Time(s),DataSize,Algorithm,Notes" << std::endl;
        }
        
        // Get current timestamp
        auto now = std::chrono::system_clock::now();
        auto now_time_t = std::chrono::system_clock::to_time_t(now);
        std::tm* now_tm = std::localtime(&now_time_t);
        
        char timeStr[20];
        std::strftime(timeStr, sizeof(timeStr), "%Y-%m-%d %H:%M:%S", now_tm);
        
        // Write data
        file << timeStr << ","
            << name << ","
            << std::fixed << std::setprecision(6) << time << ","
            << dataSize << ","
            << algorithm << ","
            << notes << std::endl;
            
        file.close();
    }
};

// Benchmark PDB parsing performance
TEST_F(BenchmarkTest, PDBParsingPerformance) {
    const int repetitions = 5;
    
    // Load PDB file and measure time
    double parseTime = measureExecutionTime([&]() {
        Molecule molecule = parser.parsePDB(pdb1Path);
        return molecule;
    }, repetitions);
    
    // Get file size for context
    std::ifstream file(pdb1Path, std::ios::binary | std::ios::ate);
    std::streamsize fileSize = file.tellg();
    file.close();
    
    // Record result
    writeBenchmarkResult("PDB_Parsing", parseTime, static_cast<int>(fileSize), 
                        "PDBParser", "Parsing " + pdb1Path);
    
    // Output result
    std::cout << "PDB parsing performance: " << std::fixed << std::setprecision(3) 
            << parseTime << " s for " << fileSize << " bytes" << std::endl;
    
    // Verify parsing worked
    Molecule molecule = parser.parsePDB(pdb1Path);
    EXPECT_GT(molecule.getResidueCount(), 0);
}

// Benchmark sequential alignment performance
TEST_F(BenchmarkTest, SequentialAlignmentPerformance) {
    const int repetitions = 3;
    
    // Load molecules
    Molecule molecule1 = parser.parsePDB(pdb1Path);
    Molecule molecule2 = parser.parsePDB(pdb2Path);
    
    int totalResidues = molecule1.getResidueCount() + molecule2.getResidueCount();
    
    // Measure sequential alignment performance
    double seqAlignTime = measureExecutionTime([&]() {
        // Create alignment instance with sequential alignment type
        std::unique_ptr<Alignment> alignment = Alignment::createSequentialAlignment(config, gplusConfig);
        alignment->setMolecules(molecule1, molecule2);
        return alignment->align();
    }, repetitions);
    
    // Record result
    writeBenchmarkResult("Sequential_Alignment", seqAlignTime, totalResidues, 
                        "SequentialAlignment", 
                        "Aligning " + molecule1.getName() + " with " + molecule2.getName());
    
    // Output result
    std::cout << "Sequential alignment performance: " << std::fixed << std::setprecision(3)
            << seqAlignTime << " s for " << totalResidues << " total residues" << std::endl;
    
    // Verify alignment worked correctly
    std::unique_ptr<Alignment> alignment = Alignment::createSequentialAlignment(config, gplusConfig);
    alignment->setMolecules(molecule1, molecule2);
    AlignmentResult result = alignment->align();
    
    EXPECT_GT(result.getAlignedCount(), 0);
}

// Benchmark non-sequential alignment performance
TEST_F(BenchmarkTest, NonSequentialAlignmentPerformance) {
    const int repetitions = 3;
    
    // Load molecules
    Molecule molecule1 = parser.parsePDB(pdb1Path);
    Molecule molecule2 = parser.parsePDB(pdb2Path);
    
    int totalResidues = molecule1.getResidueCount() + molecule2.getResidueCount();
    
    // Measure non-sequential alignment performance
    double nonSeqAlignTime = measureExecutionTime([&]() {
        // Create alignment instance with non-sequential alignment type
        std::unique_ptr<Alignment> alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
        alignment->setMolecules(molecule1, molecule2);
        return alignment->align();
    }, repetitions);
    
    // Record result
    writeBenchmarkResult("NonSequential_Alignment", nonSeqAlignTime, totalResidues, 
                        "NonSequentialAlignment", 
                        "Aligning " + molecule1.getName() + " with " + molecule2.getName());
    
    // Output result
    std::cout << "Non-sequential alignment performance: " << std::fixed << std::setprecision(3)
            << nonSeqAlignTime << " s for " << totalResidues << " total residues" << std::endl;
    
    // Verify alignment worked correctly
    std::unique_ptr<Alignment> alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
    alignment->setMolecules(molecule1, molecule2);
    AlignmentResult result = alignment->align();
    
    EXPECT_GT(result.getAlignedCount(), 0);
}

// Benchmark to compare performance with different configuration settings
TEST_F(BenchmarkTest, ConfigurationPerformanceComparison) {
    // Load molecules
    Molecule molecule1 = parser.parsePDB(pdb1Path);
    Molecule molecule2 = parser.parsePDB(pdb2Path);
    
    int totalResidues = molecule1.getResidueCount() + molecule2.getResidueCount();
    
    std::cout << "Configuration performance comparison:" << std::endl;
    
    // Test different core distance settings
    std::vector<double> coreDistances = {3.0, 4.0, 5.0, 6.0};
    
    for (double coreDistance : coreDistances) {
        config.setCoreDistance(coreDistance);
        
        double alignTime = measureExecutionTime([&]() {
            std::unique_ptr<Alignment> alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
            alignment->setMolecules(molecule1, molecule2);
            return alignment->align();
        });
        
        writeBenchmarkResult("CoreDistance_Performance", alignTime, totalResidues,
                          "NonSequentialAlignment", 
                          "CoreDistance=" + std::to_string(coreDistance));
        
        std::cout << "CoreDistance=" << coreDistance << ": " 
                << std::fixed << std::setprecision(3) << alignTime << " s" << std::endl;
    }
    
    // Reset core distance
    config.setCoreDistance(4.0);
    
    // Test different evaluation depths
    std::vector<int> evaluationDepths = {5, 10, 15, 20};
    
    for (int depth : evaluationDepths) {
        gplusConfig.setEvaluationDepth(depth);
        
        double alignTime = measureExecutionTime([&]() {
            std::unique_ptr<Alignment> alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
            alignment->setMolecules(molecule1, molecule2);
            return alignment->align();
        });
        
        writeBenchmarkResult("EvaluationDepth_Performance", alignTime, totalResidues,
                          "NonSequentialAlignment", 
                          "EvaluationDepth=" + std::to_string(depth));
        
        std::cout << "EvaluationDepth=" << depth << ": " 
                << std::fixed << std::setprecision(3) << alignTime << " s" << std::endl;
    }
}