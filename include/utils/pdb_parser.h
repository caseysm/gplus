#ifndef GANGSTA_PDB_PARSER_H
#define GANGSTA_PDB_PARSER_H

#include "core/molecule.h"
#include <string>
#include <vector>

namespace gangsta {
namespace utils {

/**
 * @brief Parser for Protein Data Bank (PDB) files
 */
class PDBParser 
{
public:
    /**
     * @brief Default constructor
     */
    PDBParser();
    
    /**
     * @brief Parse a PDB file into a molecule
     * @param filePath Path to the PDB file
     * @return Molecule containing the parsed structure
     */
    core::Molecule parsePDB(const std::string& filePath);
    
    /**
     * @brief Parse PDB content directly from a string
     * @param content PDB file content as string
     * @param name Name to give the molecule
     * @return Molecule containing the parsed structure
     */
    core::Molecule parsePDBContent(const std::string& content, const std::string& name);
    
    /**
     * @brief Write a molecule to a PDB file
     * @param molecule Molecule to write
     * @param filePath Path to the output PDB file
     * @return True if the file was written successfully
     */
    bool writePDB(const core::Molecule& molecule, const std::string& filePath);
    
    /**
     * @brief Convert a molecule to PDB file content
     * @param molecule Molecule to convert
     * @return PDB file content as string
     */
    std::string moleculeToPDBContent(const core::Molecule& molecule);
    
    /**
     * @brief Get last error message
     * @return Error message from the last operation
     */
    const std::string& getErrorMessage() const;
    
private:
    std::string errorMessage;    ///< Last error message
    
    /**
     * @brief Parse an ATOM or HETATM record
     * @param line Line from PDB file
     * @param molecule Molecule to add the atom to
     * @return True if parsing was successful
     */
    bool parseAtomRecord(const std::string& line, core::Molecule& molecule);
    
    /**
     * @brief Parse a HELIX record
     * @param line Line from PDB file
     * @param molecule Molecule to add the helix to
     * @return True if parsing was successful
     */
    bool parseHelixRecord(const std::string& line, core::Molecule& molecule);
    
    /**
     * @brief Parse a SHEET record
     * @param line Line from PDB file
     * @param molecule Molecule to add the strand to
     * @return True if parsing was successful
     */
    bool parseSheetRecord(const std::string& line, core::Molecule& molecule);
    
    /**
     * @brief Trim whitespace from a string
     * @param s String to trim
     * @return Trimmed string
     */
    std::string trim(const std::string& s);
    
    /**
     * @brief Convert string to integer with range checking
     * @param s String to convert
     * @param defaultValue Default value if conversion fails
     * @return Converted integer
     */
    int safeStoi(const std::string& s, int defaultValue = 0);
    
    /**
     * @brief Convert string to double with range checking
     * @param s String to convert
     * @param defaultValue Default value if conversion fails
     * @return Converted double
     */
    double safeStod(const std::string& s, double defaultValue = 0.0);
};

} // namespace utils
} // namespace gangsta

#endif // GANGSTA_PDB_PARSER_H