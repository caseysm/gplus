#ifndef GANGSTA_MOLECULE_H
#define GANGSTA_MOLECULE_H

#include "core/atom.h"
#include "core/residue.h"
#include "core/sse.h"
#include "math/vec3.h"
#include "math/rotation.h"
#include <string>
#include <vector>
#include <memory>

namespace gangsta {
namespace core {

/**
 * @brief Complete molecular structure representation
 * 
 * Represents a protein structure with atoms, residues, and
 * secondary structure elements (SSEs).
 */
class Molecule 
{
public:
    /**
     * @brief Default constructor
     */
    Molecule();
    
    /**
     * @brief Constructor with name
     * @param name Molecule name/identifier
     */
    explicit Molecule(const std::string& name);
    
    /**
     * @brief Get molecule name
     * @return Molecule name/identifier
     */
    const std::string& getName() const;
    
    /**
     * @brief Set molecule name
     * @param name Molecule name/identifier
     */
    void setName(const std::string& name);
    
    /**
     * @brief Add an atom to the molecule
     * @param atom Atom to add
     */
    void addAtom(const Atom& atom);
    
    /**
     * @brief Get all atoms in the molecule
     * @return Vector of atoms
     */
    const std::vector<Atom>& getAtoms() const;
    
    /**
     * @brief Get atoms marked as hetero atoms (non-standard)
     * @return Vector of hetero atoms
     */
    const std::vector<Atom>& getHeteroAtoms() const;
    
    /**
     * @brief Add a hetero atom to the molecule
     * @param atom Hetero atom to add
     */
    void addHeteroAtom(const Atom& atom);
    
    /**
     * @brief Add a residue to the molecule
     * @param residue Residue to add
     */
    void addResidue(const Residue& residue);
    
    /**
     * @brief Get all residues in the molecule
     * @return Vector of residues
     */
    const std::vector<Residue>& getResidues() const;
    
    /**
     * @brief Add a secondary structure element
     * @param sse SSE to add
     */
    void addSSE(const SSE& sse);
    
    /**
     * @brief Get all secondary structure elements
     * @return Vector of SSEs
     */
    const std::vector<SSE>& getSSEs() const;
    
    /**
     * @brief Get all alpha carbons
     * @return Vector of alpha carbon atoms
     */
    const std::vector<Atom>& getAlphaCarbons() const;
    
    /**
     * @brief Detect secondary structure elements if not defined
     * @return Number of SSEs detected
     */
    int detectSSEs();
    
    /**
     * @brief Assign atoms to their corresponding secondary structure elements
     */
    void assignAtomsToSSEs();
    
    /**
     * @brief Calculate center of mass
     * @return Center of mass position
     */
    math::Vec3 getCenterOfMass() const;
    
    /**
     * @brief Apply a transformation to all atoms
     * @param rotation Rotation to apply
     * @param translation Translation to apply after rotation
     */
    void transform(const math::Rotation& rotation, const math::Vec3& translation);
    
    /**
     * @brief Calculate RMSD between this molecule and another
     * @param other Another molecule
     * @param alignedIndices Vector of paired atom indices
     * @return Root mean square deviation
     */
    double calculateRMSD(const Molecule& other, const std::vector<std::pair<int, int>>& alignedIndices) const;
    
    /**
     * @brief Finalize the molecule structure
     * 
     * This should be called after adding all atoms to ensure
     * internal data structures are consistent.
     */
    void finalize();
    
    /**
     * @brief Get number of atoms
     * @return Atom count
     */
    size_t getAtomCount() const;
    
    /**
     * @brief Get number of residues
     * @return Residue count
     */
    size_t getResidueCount() const;
    
    /**
     * @brief Get number of SSEs
     * @return SSE count
     */
    size_t getSSECount() const;
    
    /**
     * @brief Get information about the molecule
     * @return Formatted string with molecule information
     */
    std::string getInfo() const;
    
private:
    std::string name;                ///< Molecule name/identifier
    std::vector<Atom> atoms;         ///< All standard atoms
    std::vector<Atom> heteroAtoms;   ///< Non-standard (HETATM) atoms
    std::vector<Residue> residues;   ///< All residues
    std::vector<SSE> sses;           ///< Secondary structure elements
    std::vector<Atom> alphaCarbons;  ///< Alpha carbon atoms
    
    /**
     * @brief Extract residues from atom list
     */
    void extractResidues();
    
    /**
     * @brief Extract alpha carbons from atoms
     */
    void extractAlphaCarbons();
};

} // namespace core
} // namespace gangsta

#endif // GANGSTA_MOLECULE_H