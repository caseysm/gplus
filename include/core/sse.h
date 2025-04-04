#ifndef GANGSTA_SSE_H
#define GANGSTA_SSE_H

#include <string>

namespace gangsta {
namespace core {

/**
 * @brief Secondary Structure Element representation
 * 
 * Represents a protein secondary structure element (SSE) such as
 * alpha-helix or beta-strand.
 */
class SSE 
{
public:
    // Secondary structure types
    enum Type {
        UNKNOWN = 0,
        HELIX = 1,
        STRAND = 2,
        LOOP = 3
    };
    
    /**
     * @brief Default constructor
     */
    SSE();
    
    /**
     * @brief Constructor with type and residue range
     * @param type SSE type (helix, strand, etc.)
     * @param startResidueIndex Start residue index
     * @param endResidueIndex End residue index
     */
    SSE(Type type, int startResidueIndex, int endResidueIndex);
    
    /**
     * @brief Get the type of secondary structure
     * @return SSE type
     */
    Type getType() const;
    
    /**
     * @brief Set the type of secondary structure
     * @param type SSE type
     */
    void setType(Type type);
    
    /**
     * @brief Get the starting residue index
     * @return Start residue index
     */
    int getStart() const;
    
    /**
     * @brief Set the starting residue index
     * @param index Start residue index
     */
    void setStart(int index);
    
    /**
     * @brief Get the ending residue index
     * @return End residue index
     */
    int getEnd() const;
    
    /**
     * @brief Set the ending residue index
     * @param index End residue index
     */
    void setEnd(int index);
    
    /**
     * @brief Get the length of the SSE in residues
     * @return Length of the SSE (end - start + 1)
     */
    int getLength() const;
    
    /**
     * @brief Get string representation of the SSE type
     * @return String name ("Helix", "Strand", "Loop", or "Unknown")
     */
    std::string getTypeName() const;
    
    /**
     * @brief Get one letter code for SSE type
     * @return Single character representation (H, E, L, or -)
     */
    char getTypeCode() const;
    
    /**
     * @brief Check if a residue index is within this SSE
     * @param residueIndex The residue index to check
     * @return True if the residue is part of this SSE
     */
    bool containsResidue(int residueIndex) const;
    
private:
    Type type;            ///< Type of secondary structure
    int startResidueIndex; ///< Starting residue index
    int endResidueIndex;   ///< Ending residue index
};

} // namespace core
} // namespace gangsta

#endif // GANGSTA_SSE_H