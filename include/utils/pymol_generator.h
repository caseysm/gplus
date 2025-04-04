#ifndef GANGSTA_PYMOL_GENERATOR_H
#define GANGSTA_PYMOL_GENERATOR_H

#include "core/molecule.h"
#include "algorithm/alignment.h"
#include <string>

namespace gangsta {
namespace utils {

/**
 * @brief Generates PyMOL visualization scripts for aligned structures
 */
class PyMolGenerator 
{
public:
    /**
     * @brief Default constructor
     */
    PyMolGenerator();
    
    /**
     * @brief Set the molecules for visualization
     * @param molecule1 First molecule
     * @param molecule2 Second molecule
     * @param alignedMolecule1 Aligned version of first molecule
     */
    void setMolecules(const core::Molecule& molecule1, 
                     const core::Molecule& molecule2,
                     const core::Molecule& alignedMolecule1);
    
    /**
     * @brief Set alignment result
     * @param result Alignment result
     */
    void setAlignmentResult(const algorithm::AlignmentResult& result);
    
    /**
     * @brief Generate PyMOL script
     * @param outputPath Path to save the script
     * @param pdbOutputPath Path to save the aligned structures as PDB
     * @return True if script was generated successfully
     */
    bool generateScript(const std::string& outputPath, const std::string& pdbOutputPath);
    
    /**
     * @brief Generate PyMOL script content
     * @param pdbPath Path to the PDB file with aligned structures
     * @return PyMOL script content
     */
    std::string generateScriptContent(const std::string& pdbPath);
    
private:
    core::Molecule molecule1;           ///< First molecule
    core::Molecule molecule2;           ///< Second molecule
    core::Molecule alignedMolecule1;    ///< Aligned version of first molecule
    algorithm::AlignmentResult result;  ///< Alignment result
    
    /**
     * @brief Generate PDB file with aligned structures
     * @param outputPath Path to save the PDB file
     * @return True if PDB file was generated successfully
     */
    bool generateAlignedPDB(const std::string& outputPath);
    
    /**
     * @brief Generate color scheme settings
     * @return PyMOL color scheme commands
     */
    std::string generateColorScheme() const;
    
    /**
     * @brief Generate visualization settings
     * @return PyMOL visualization commands
     */
    std::string generateVisualizationSettings() const;
    
    /**
     * @brief Generate alignment display commands
     * @return PyMOL alignment display commands
     */
    std::string generateAlignmentCommands() const;
};

} // namespace utils
} // namespace gangsta

#endif // GANGSTA_PYMOL_GENERATOR_H