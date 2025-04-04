#ifndef GANGSTA_ATOM_H
#define GANGSTA_ATOM_H

#include "math/vec3.h"
#include <string>

namespace gangsta {
namespace core {

/**
 * @brief Atom representation in protein structure
 */
class Atom 
{
public:
    /**
     * @brief Default constructor
     */
    Atom();
    
    /**
     * @brief Constructor with position and names
     * @param position 3D coordinates
     * @param atomName Atom name (e.g., "CA", "CB")
     * @param residueName Residue name (e.g., "ALA", "GLY")
     * @param residueNumber Residue sequence number
     * @param chainID Chain identifier (e.g., "A", "B")
     */
    Atom(const math::Vec3& position, 
         const std::string& atomName, 
         const std::string& residueName, 
         int residueNumber,
         char chainID = ' ');
    
    /**
     * @brief Get atom position
     * @return 3D position
     */
    const math::Vec3& getPosition() const;
    
    /**
     * @brief Set atom position
     * @param position 3D position
     */
    void setPosition(const math::Vec3& position);
    
    /**
     * @brief Get atom name
     * @return Atom name (e.g., "CA", "CB")
     */
    const std::string& getAtomName() const;
    
    /**
     * @brief Set atom name
     * @param name Atom name
     */
    void setAtomName(const std::string& name);
    
    /**
     * @brief Get residue name
     * @return Residue name (e.g., "ALA", "GLY")
     */
    const std::string& getResidueName() const;
    
    /**
     * @brief Set residue name
     * @param name Residue name
     */
    void setResidueName(const std::string& name);
    
    /**
     * @brief Get residue number
     * @return Residue sequence number
     */
    int getResidueNumber() const;
    
    /**
     * @brief Set residue number
     * @param number Residue sequence number
     */
    void setResidueNumber(int number);
    
    /**
     * @brief Get chain identifier
     * @return Chain ID character
     */
    char getChainID() const;
    
    /**
     * @brief Set chain identifier
     * @param chainID Chain ID character
     */
    void setChainID(char chainID);
    
    /**
     * @brief Get residue code (1-letter amino acid code)
     * @return Single letter amino acid code
     */
    char getResidueCode() const;
    
    /**
     * @brief Check if this is an alpha carbon
     * @return True if this is a CA atom
     */
    bool isAlphaCarbon() const;
    
    /**
     * @brief Get the secondary structure element index this atom belongs to
     * @return SSE index (-1 if not assigned)
     */
    int getSSEIndex() const;
    
    /**
     * @brief Set the secondary structure element index
     * @param index SSE index
     */
    void setSSEIndex(int index);
    
    /**
     * @brief Calculate distance to another atom
     * @param other Another atom
     * @return Euclidean distance
     */
    double distanceTo(const Atom& other) const;
    
private:
    math::Vec3 position;       ///< 3D coordinates
    std::string atomName;      ///< Atom name (e.g., "CA", "CB")
    std::string residueName;   ///< Residue name (e.g., "ALA", "GLY")
    int residueNumber;         ///< Residue sequence number
    char chainID;              ///< Chain identifier
    int sseIndex;              ///< Secondary structure element index
    
    /**
     * @brief Convert 3-letter amino acid code to 1-letter code
     * @param residueName 3-letter amino acid code
     * @return 1-letter amino acid code
     */
    static char convertResidueCode(const std::string& residueName);
};

} // namespace core
} // namespace gangsta

#endif // GANGSTA_ATOM_H