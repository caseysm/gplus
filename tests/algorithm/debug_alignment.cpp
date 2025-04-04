#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <cstdlib>

#include "config/config.h"
#include "config/gplus_config.h"
#include "core/molecule.h"
#include "algorithm/alignment.h"
#include "utils/pdb_parser.h"
#include "utils/logger.h"

using namespace gangsta::config;
using namespace gangsta::core;
using namespace gangsta::algorithm;
using namespace gangsta::utils;

/**
 * @brief Debug program for alignment analysis
 * 
 * This program is used to debug the alignment algorithm by tracing its execution
 * and outputting debug information to a file.
 */
int main(int argc, char* argv[]) {
    // Initialize logger
    Logger::getInstance().initialize("", Logger::WARNING, false);
    
    // Check command-line arguments
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <pdb1> <pdb2> <output_file> [--sequential]" << std::endl;
        return 1;
    }
    
    // Parse arguments
    std::string pdb1Path = argv[1];
    std::string pdb2Path = argv[2];
    std::string outputPath = argv[3];
    bool sequential = false;
    
    // Check for optional sequential mode
    if (argc > 4 && std::string(argv[4]) == "--sequential") {
        sequential = true;
    }
    
    // Open output file
    std::ofstream outputFile(outputPath);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open output file " << outputPath << std::endl;
        return 1;
    }
    
    // Parse PDB files
    PDBParser parser;
    Molecule molecule1 = parser.parsePDB(pdb1Path);
    Molecule molecule2 = parser.parsePDB(pdb2Path);
    
    // Check molecules were loaded correctly
    if (molecule1.getResidueCount() == 0) {
        std::cerr << "Error: Failed to load molecule from " << pdb1Path << std::endl;
        return 1;
    }
    
    if (molecule2.getResidueCount() == 0) {
        std::cerr << "Error: Failed to load molecule from " << pdb2Path << std::endl;
        return 1;
    }
    
    // Create configuration
    Config config;
    GPlusConfig gplusConfig;
    
    // Use default values from main app
    config.setCoreDistance(4);
    config.setResidueDistance(6.0);
    gplusConfig.setResultCount(200);
    gplusConfig.setEvaluationDepth(5000);
    
    // Create alignment
    std::unique_ptr<Alignment> alignment;
    if (sequential) {
        alignment = Alignment::createSequentialAlignment(config, gplusConfig);
        std::cout << "Using sequential alignment mode" << std::endl;
    } else {
        alignment = Alignment::createNonSequentialAlignment(config, gplusConfig);
        std::cout << "Using non-sequential alignment mode" << std::endl;
    }
    
    // Set debug output
    Alignment::setDebugOutput(outputFile);
    alignment->enableDebugMode(true);
    
    // Set molecules
    alignment->setMolecules(molecule1, molecule2);
    
    // Run alignment with debug output
    std::cout << "Running alignment with debug output..." << std::endl;
    AlignmentResult result = alignment->align();
    
    // Print results to console
    std::cout << "\nAlignment Results:" << std::endl;
    std::cout << "Aligned residues: " << result.getAlignedCount() << std::endl;
    std::cout << "RMSD: " << result.getRMSD() << std::endl;
    std::cout << "Score: " << result.getScore() << std::endl;
    
    // Summary in debug file
    outputFile << "\n\n======== ALIGNMENT SUMMARY ========" << std::endl;
    outputFile << "PDB1: " << pdb1Path << std::endl;
    outputFile << "PDB2: " << pdb2Path << std::endl;
    outputFile << "Mode: " << (sequential ? "Sequential" : "Non-sequential") << std::endl;
    outputFile << "Aligned residues: " << result.getAlignedCount() << std::endl;
    outputFile << "RMSD: " << result.getRMSD() << std::endl;
    outputFile << "Score: " << result.getScore() << std::endl;
    
    std::cout << "\nDebug output written to " << outputPath << std::endl;
    outputFile.close();
    
    return 0;
}