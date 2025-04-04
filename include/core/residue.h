#ifndef GANGSTA_RESIDUE_H
#define GANGSTA_RESIDUE_H

#include "core/atom.h"
#include <string>
#include <vector>
#include <memory>

namespace gangsta {
namespace core {

/**
 * @brief Represents an amino acid residue in a protein structure
 */
class Residue 
{
public:
    /**
     * @brief Default constructor
     */
    Residue();
    
    /**
     * @brief Constructor with residue information
     * @param name Three-letter amino acid code (e.g., "ALA")
     * @param number Residue sequence number
     * @param chainID Chain identifier
     */
    Residue(const std::string& name, int number, char chainID = ' ');
    
    /**
     * @brief Add an atom to the residue
     * @param atom Atom to add
     */
    void addAtom(const Atom& atom);
    
    /**
     * @brief Get atoms in the residue
     * @return Vector of atoms
     */
    const std::vector<Atom>& getAtoms() const;
    
    /**
     * @brief Get atom by name
     * @param name Atom name (e.g., "CA", "CB")
     * @return Pointer to atom if found, nullptr otherwise
     */
    Atom* getAtomByName(const std::string& name);
    
    /**
     * @brief Get atom by name (const version)
     * @param name Atom name (e.g., "CA", "CB")
     * @return Pointer to atom if found, nullptr otherwise
     */
    const Atom* getAtomByName(const std::string& name) const;
    
    /**
     * @brief Get the C-alpha atom
     * @return Pointer to C-alpha atom if found, nullptr otherwise
     */
    Atom* getAlphaCarbon();
    
    /**
     * @brief Get the C-alpha atom (const version)
     * @return Pointer to C-alpha atom if found, nullptr otherwise
     */
    const Atom* getAlphaCarbon() const;
    
    /**
     * @brief Get three-letter residue name
     * @return Three-letter amino acid code
     */
    const std::string& getName() const;
    
    /**
     * @brief Get one-letter residue code
     * @return Single letter amino acid code
     */
    char getCode() const;
    
    /**
     * @brief Get residue sequence number
     * @return Residue number
     */
    int getNumber() const;
    
    /**
     * @brief Get chain identifier
     * @return Chain ID character
     */
    char getChainID() const;
    
    /**
     * @brief Set SSE index for all atoms in this residue
     * @param index Index of the secondary structure element
     */
    void setSSEIndex(int index);
    
    /**
     * @brief Get the secondary structure element index
     * @return SSE index (-1 if not assigned)
     */
    int getSSEIndex() const;
    
    /**
     * @brief Get residue center (average position of all atoms)
     * @return Center position
     */
    math::Vec3 getCenter() const;
    
private:
    std::string name;              ///< Three-letter amino acid code
    int number;                    ///< Residue sequence number
    char chainID;                  ///< Chain identifier
    int sseIndex;                  ///< Secondary structure element index
    std::vector<Atom> atoms;       ///< Atoms in this residue
};

} // namespace core
} // namespace gangsta

#endif // GANGSTA_RESIDUE_H