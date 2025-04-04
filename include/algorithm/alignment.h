#ifndef GANGSTA_ALIGNMENT_H
#define GANGSTA_ALIGNMENT_H

#include "core/molecule.h"
#include "math/rotation.h"
#include "algorithm/kabsch.h"
#include "config/config.h"
#include "config/gplus_config.h"
#include <vector>
#include <memory>

namespace gangsta {
namespace algorithm {

/**
 * @brief Result of a protein structure alignment
 */
class AlignmentResult 
{
public:
    /**
     * @brief Default constructor
     */
    AlignmentResult();
    
    /**
     * @brief Constructor with all parameters
     * @param rmsd The RMSD value
     * @param score The score value
     * @param alignedCount Number of aligned residues
     * @param alignedPairs Vector of aligned residue pairs
     * @param rotation Rotation to apply
     * @param translation Translation to apply
     */
    AlignmentResult(double rmsd, double score, int alignedCount,
                   const std::vector<std::pair<int, int>>& alignedPairs,
                   const math::Rotation& rotation,
                   const math::Vec3& translation);
    
    /**
     * @brief Get RMSD (root mean square deviation)
     * @return RMSD value
     */
    double getRMSD() const;
    
    /**
     * @brief Set RMSD
     * @param rmsd RMSD value
     */
    void setRMSD(double rmsd);
    
    /**
     * @brief Get alignment score
     * @return Alignment score
     */
    double getScore() const;
    
    /**
     * @brief Set alignment score
     * @param score Alignment score
     */
    void setScore(double score);
    
    /**
     * @brief Get number of aligned residues
     * @return Number of aligned residues
     */
    int getAlignedCount() const;
    
    /**
     * @brief Set number of aligned residues
     * @param count Number of aligned residues
     */
    void setAlignedCount(int count);
    
    /**
     * @brief Get rotation matrix for alignment
     * @return Rotation matrix
     */
    const math::Rotation& getRotation() const;
    
    /**
     * @brief Set rotation matrix
     * @param rotation Rotation matrix
     */
    void setRotation(const math::Rotation& rotation);
    
    /**
     * @brief Get translation vector for alignment
     * @return Translation vector
     */
    const math::Vec3& getTranslation() const;
    
    /**
     * @brief Set translation vector
     * @param translation Translation vector
     */
    void setTranslation(const math::Vec3& translation);
    
    /**
     * @brief Get aligned residue pairs
     * @return Vector of aligned residue indices (molecule1_index, molecule2_index)
     */
    const std::vector<std::pair<int, int>>& getAlignedPairs() const;
    
    /**
     * @brief Set aligned residue pairs
     * @param pairs Vector of aligned residue indices
     */
    void setAlignedPairs(const std::vector<std::pair<int, int>>& pairs);
    
    /**
     * @brief Add an aligned residue pair
     * @param molecule1Index Index in first molecule
     * @param molecule2Index Index in second molecule
     */
    void addAlignedPair(int molecule1Index, int molecule2Index);
    
    /**
     * @brief Clear all aligned pairs
     */
    void clearAlignedPairs();
    
    /**
     * @brief Compare with another result
     * @param other Another alignment result
     * @return True if this result has a better score
     */
    bool isBetterThan(const AlignmentResult& other) const;
    
    /**
     * @brief Get a summary of the alignment result
     * @return Formatted string with alignment information
     */
    std::string getSummary() const;
    
private:
    double rmsd;                               ///< RMSD of alignment
    double score;                              ///< Alignment score
    int alignedCount;                          ///< Number of aligned residues
    math::Rotation rotation;                   ///< Rotation matrix
    math::Vec3 translation;                    ///< Translation vector
    std::vector<std::pair<int, int>> alignedPairs; ///< Aligned residue pairs
};

/**
 * @brief Base class for protein structure alignment algorithms
 */
class Alignment 
{
public:
    /**
     * @brief Default constructor
     */
    Alignment();
    
    /**
     * @brief Constructor with configuration
     * @param config Alignment configuration
     * @param gplusConfig GPlusConfig configuration
     */
    Alignment(const config::Config& config, const config::GPlusConfig& gplusConfig);
    
    /**
     * @brief Enable or disable debug output for algorithm tracing
     * @param enable Whether to enable debug output
     */
    void enableDebugMode(bool enable);
    
    /**
     * @brief Set the output stream for debug messages
     * @param output Reference to output stream
     */
    static void setDebugOutput(std::ostream& output);
    
    /**
     * @brief Virtual destructor
     */
    virtual ~Alignment();
    
    /**
     * @brief Set the molecules to align
     * @param molecule1 First molecule
     * @param molecule2 Second molecule
     */
    void setMolecules(const core::Molecule& molecule1, const core::Molecule& molecule2);
    
    /**
     * @brief Get the first molecule
     * @return First molecule
     */
    const core::Molecule& getMolecule1() const;
    
    /**
     * @brief Get the second molecule
     * @return Second molecule
     */
    const core::Molecule& getMolecule2() const;
    
    /**
     * @brief Set alignment configuration
     * @param config Alignment configuration
     */
    void setConfig(const config::Config& config);
    
    /**
     * @brief Get alignment configuration
     * @return Alignment configuration
     */
    const config::Config& getConfig() const;
    
    /**
     * @brief Set GPlusConfig configuration
     * @param config GPlusConfig configuration
     */
    void setGPlusConfig(const config::GPlusConfig& config);
    
    /**
     * @brief Get GPlusConfig configuration
     * @return GPlusConfig configuration
     */
    const config::GPlusConfig& getGPlusConfig() const;
    
    /**
     * @brief Perform the alignment
     * @return Alignment result
     */
    virtual AlignmentResult align() = 0;
    
    /**
     * @brief Apply the alignment result to a molecule
     * @param molecule Molecule to transform
     * @param result Alignment result to apply
     * @return Transformed molecule
     */
    static core::Molecule applyAlignment(const core::Molecule& molecule, const AlignmentResult& result);
    
    /**
     * @brief Create a sequential alignment instance
     * @param config Alignment configuration
     * @param gplusConfig GPlusConfig configuration
     * @return Unique pointer to sequential alignment
     */
    static std::unique_ptr<Alignment> createSequentialAlignment(
        const config::Config& config, 
        const config::GPlusConfig& gplusConfig);
        
    /**
     * @brief Create a non-sequential alignment instance
     * @param config Alignment configuration
     * @param gplusConfig GPlusConfig configuration
     * @return Unique pointer to non-sequential alignment
     */
    static std::unique_ptr<Alignment> createNonSequentialAlignment(
        const config::Config& config, 
        const config::GPlusConfig& gplusConfig);
    
protected:
    core::Molecule molecule1;          ///< First molecule
    core::Molecule molecule2;          ///< Second molecule
    config::Config config;             ///< Alignment configuration
    config::GPlusConfig gplusConfig;   ///< GPlusConfig configuration
    Kabsch kabsch;                     ///< Kabsch algorithm for RMSD calculation
    bool debugEnabled;                 ///< Flag indicating if debug mode is enabled
    static std::ostream* debugStream;  ///< Stream for debug output
    
    /**
     * @brief Log a debug message
     * @param message Message to log
     */
    void debugLog(const std::string& message) const;
    
    /**
     * @brief Log molecule information for debugging
     * @param molecule Molecule to log
     * @param name Name to identify the molecule in the log
     */
    void debugLogMolecule(const core::Molecule& molecule, const std::string& name) const;
    
    /**
     * @brief Log SSE information for debugging
     * @param sseList List of SSEs to log
     * @param name Name to identify the SSE list in the log
     */
    void debugLogSSEs(const std::vector<core::SSE>& sseList, const std::string& name) const;
    
    /**
     * @brief Log alignment pairs for debugging
     * @param pairs Vector of alignment pairs
     * @param description Description of alignment pairs
     */
    void debugLogAlignmentPairs(const std::vector<std::pair<int, int>>& pairs, 
                              const std::string& description) const;
    
    /**
     * @brief Calculate alignment score
     * @param rmsd RMSD of alignment
     * @param alignedCount Number of aligned residues
     * @param totalCount Total number of residues
     * @return Alignment score
     */
    double calculateScore(double rmsd, int alignedCount, int totalCount) const;
    
    /**
     * @brief Calculate RMSD of a set of aligned residue pairs
     * @param alignedPairs Vector of aligned residue indices
     * @return RMSD value and transformation
     */
    AlignmentResult calculateRMSDAndTransform(const std::vector<std::pair<int, int>>& alignedPairs);
    
    /**
     * @brief Check if two residues can be aligned
     * @param residue1 First residue
     * @param residue2 Second residue
     * @return True if residues can be aligned
     */
    bool canAlignResidues(const core::Residue& residue1, const core::Residue& residue2) const;
};

/**
 * @brief Sequential alignment algorithm
 */
class SequentialAlignment : public Alignment 
{
public:
    /**
     * @brief Default constructor
     */
    SequentialAlignment();
    
    /**
     * @brief Constructor with configuration
     * @param config Alignment configuration
     * @param gplusConfig GPlusConfig configuration
     */
    SequentialAlignment(const config::Config& config, const config::GPlusConfig& gplusConfig);
    
    /**
     * @brief Perform sequential alignment
     * @return Alignment result
     */
    AlignmentResult align() override;
    
private:
    /**
     * @brief Align molecules using a sequential approach
     * @return Alignment result
     */
    AlignmentResult alignSequential();
};

/**
 * @brief Non-sequential alignment algorithm (GANGSTA+)
 */
class NonSequentialAlignment : public Alignment 
{
public:
    /**
     * @brief Default constructor
     */
    NonSequentialAlignment();
    
    /**
     * @brief Constructor with configuration
     * @param config Alignment configuration
     * @param gplusConfig GPlusConfig configuration
     */
    NonSequentialAlignment(const config::Config& config, const config::GPlusConfig& gplusConfig);
    
    /**
     * @brief Perform non-sequential alignment
     * @return Alignment result
     */
    AlignmentResult align() override;
    
private:
    /**
     * @brief Align secondary structure elements
     * @return Vector of aligned SSE pairs
     */
    std::vector<std::pair<int, int>> alignSSEs();
    
    /**
     * @brief Extend SSE alignment to include residues
     * @param sseAlignments Vector of aligned SSE pairs
     * @return Alignment result
     */
    AlignmentResult extendSSEAlignment(const std::vector<std::pair<int, int>>& sseAlignments);
    
    /**
     * @brief Optimize alignment by iterative refinement
     * @param initialResult Initial alignment result
     * @return Optimized alignment result
     */
    AlignmentResult optimizeAlignment(const AlignmentResult& initialResult);
    
    /**
     * @brief Helper function to extract points for a subset of atoms
     * @param atoms Vector of atoms
     * @param pairs Vector of atom index pairs
     * @param useFirst Whether to use the first or second element of the pairs
     * @return Vector of atom positions
     */
    std::vector<math::Vec3> transformPoints(const std::vector<core::Atom>& atoms, 
                                         const std::vector<std::pair<int, int>>& pairs,
                                         bool useFirst);
};

} // namespace algorithm
} // namespace gangsta

#endif // GANGSTA_ALIGNMENT_H