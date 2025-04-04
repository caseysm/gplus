#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <cstring>
#include <libgen.h>
#include <fstream>

#include "config/config.h"
#include "config/gplus_config.h"
#include "core/molecule.h"
#include "utils/logger.h"
#include "utils/file.h"
#include "utils/pdb_parser.h"
#include "utils/pymol_generator.h"
#include "algorithm/kabsch.h"
#include "algorithm/alignment.h"

// Constants
#define PRODUCT_NAME "GANGSTA+ - Non-Sequential Protein Structure Alignment"
#define PRODUCT_COPY "Freie Universitaet Berlin (Biochemistry, AG Knapp)"
#define PRODUCT_AUTHORS "Aysam Guerler (aysam.guerler@gmail.com)"
#define PRODUCT_VERSION "3.08 - Modular Edition"

/**
 * @brief Extract the file name from a path
 * @param path Full path
 * @return File name without directory
 */
std::string getFileName(const std::string& path) {
    char* pathCopy = strdup(path.c_str());
    std::string result = basename(pathCopy);
    free(pathCopy);
    return result;
}

/**
 * @brief Extract the file base name (without extension)
 * @param fileName File name
 * @return Base name without extension
 */
std::string getBaseName(const std::string& fileName) {
    size_t pos = fileName.find_last_of('.');
    if (pos == std::string::npos) {
        return fileName;
    }
    return fileName.substr(0, pos);
}

/**
 * @brief Print program banner
 */
void printBanner() {
    std::cout << " **************************************************************************\n";
    std::cout << " *                              GANGSTA+                                  *\n";
    std::cout << " *                                                                        *\n";
    std::cout << " * Reference: A. Guerler and E. W. Knapp, Protein Sci. 2008 17, 1374-82   *\n";
    std::cout << " * Comments on the program, please email to: aysam.guerler@gmail.com      *\n";
    std::cout << " **************************************************************************\n\n";
}

/**
 * @brief Print usage instructions
 * @param programName Name of the program executable
 */
void printUsage(const std::string& programName) {
    std::cout << "Usage: " << programName << " <pdb1> <pdb2> [options]\n";
    std::cout << "\nOptions:\n";
    std::cout << "  --sequential               Use sequential alignment (default: non-sequential)\n";
    std::cout << "  --core-distance <value>    Set core distance (default: 4)\n";
    std::cout << "  --residue-distance <value> Set residue distance (default: 6)\n";
    std::cout << "  --allow-inversion          Allow inversion of SSEs (default: false)\n";
    std::cout << "  --no-pdb                   Don't generate PDB and PyMOL output files\n";
    std::cout << "  --debug                    Enable detailed debugging output\n";
    std::cout << "  --debug-file <filename>    Write debug output to the specified file\n";
    std::cout << "  --help                     Show this help message\n";
}

/**
 * @brief Main function
 */
int main(int argc, char* argv[]) {
    // Initialize logger
    gangsta::utils::Logger::getInstance().initialize("", gangsta::utils::Logger::INFO, true);
    
    // Print banner
    printBanner();
    
    // Parse command line arguments
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }
    
    // Get input files
    std::string pdbPath1 = argv[1];
    std::string pdbPath2 = argv[2];
    
    // Parse options
    bool useSequential = false;
    bool generateOutput = true;
    bool enableDebug = false;
    std::string debugFilePath = "debug_output.txt";
    
    // Create configuration objects with default values
    gangsta::config::Config config;
    gangsta::config::GPlusConfig gplusConfig;
    
    // Parse additional options
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--sequential") {
            useSequential = true;
        } else if (arg == "--core-distance" && i + 1 < argc) {
            config.setCoreDistance(std::atoi(argv[++i]));
        } else if (arg == "--residue-distance" && i + 1 < argc) {
            config.setResidueDistance(std::atof(argv[++i]));
        } else if (arg == "--allow-inversion") {
            config.setInversion(true);
        } else if (arg == "--no-pdb") {
            generateOutput = false;
        } else if (arg == "--debug") {
            enableDebug = true;
        } else if (arg == "--debug-file" && i + 1 < argc) {
            debugFilePath = argv[++i];
            enableDebug = true;
        } else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cout << "Unknown option: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // Display input information
    std::cout << "Chain 1:         " << pdbPath1 << "\n";
    std::cout << "Chain 2:         " << pdbPath2 << "\n";
    std::cout << "Features: " << (useSequential ? "sequential" : "non-sequential")
              << ", " << (config.getInversion() ? "" : "no-") << "inver." << "\n";
    
    // Check if input files exist
    if (!gangsta::utils::File::exists(pdbPath1)) {
        std::cout << "Error: File not found: " << pdbPath1 << "\n";
        return 1;
    }
    
    if (!gangsta::utils::File::exists(pdbPath2)) {
        std::cout << "Error: File not found: " << pdbPath2 << "\n";
        return 1;
    }
    
    // Parse PDB files
    gangsta::utils::PDBParser parser;
    gangsta::core::Molecule molecule1 = parser.parsePDB(pdbPath1);
    gangsta::core::Molecule molecule2 = parser.parsePDB(pdbPath2);
    
    // Display molecule information
    std::cout << "Chain 1: " << getFileName(pdbPath1) << "      Size= " << molecule1.getResidueCount() << "\n";
    std::cout << "Chain 2: " << getFileName(pdbPath2) << "      Size= " << molecule2.getResidueCount() << "\n";
    
    // Perform alignment
    std::unique_ptr<gangsta::algorithm::Alignment> alignment;
    
    if (useSequential) {
        alignment.reset(new gangsta::algorithm::SequentialAlignment(config, gplusConfig));
    } else {
        alignment.reset(new gangsta::algorithm::NonSequentialAlignment(config, gplusConfig));
    }
    
    // Add debug mode for algorithm tracing (for validation)
    if (enableDebug) {
        std::cout << "Debug mode enabled, writing to " << debugFilePath << "\n";
        std::ofstream* debugOutput = new std::ofstream(debugFilePath);
        if (debugOutput->is_open()) {
            gangsta::algorithm::Alignment::setDebugOutput(*debugOutput);
            alignment->enableDebugMode(true);
        } else {
            std::cout << "Warning: Failed to open debug output file: " << debugFilePath << "\n";
            delete debugOutput;
        }
    }
    
    alignment->setMolecules(molecule1, molecule2);
    
    // Run alignment
    gangsta::algorithm::AlignmentResult result = alignment->align();
    
    // Display alignment result
    std::cout << result.getSummary() << "\n\n";
    
    // Generate output files if requested
    if (generateOutput) {
        // Generate aligned molecule
        gangsta::core::Molecule alignedMolecule1 = gangsta::algorithm::Alignment::applyAlignment(
            molecule1, result);
        
        // Create output file names
        std::string baseName = getBaseName(getFileName(pdbPath1)) + "_" + getBaseName(getFileName(pdbPath2));
        std::string pdbOutputPath = baseName + ".pdb";
        std::string pyMolOutputPath = baseName + ".pml";
        
        // Generate PDB file with aligned structures
        gangsta::utils::PyMolGenerator pyMolGenerator;
        pyMolGenerator.setMolecules(molecule1, molecule2, alignedMolecule1);
        pyMolGenerator.setAlignmentResult(result);
        
        if (pyMolGenerator.generateScript(pyMolOutputPath, pdbOutputPath)) {
            std::cout << "Output files generated:\n";
            std::cout << "  PDB file: " << pdbOutputPath << "\n";
            std::cout << "  PyMOL script: " << pyMolOutputPath << "\n";
        } else {
            std::cout << "Error: Failed to generate output files\n";
        }
    }
    
    std::cout << "\nGANGSTA+ Done.\n";
    
    return 0;
}