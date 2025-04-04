#ifndef GANGSTA_CONFIG_H
#define GANGSTA_CONFIG_H

namespace gangsta {
namespace config {

/**
 * @brief Configuration for alignment detection parameters
 */
class Config 
{
public:
    /**
     * @brief Constructor with default values
     */
    Config();
    
    /**
     * @brief Set core distance parameter
     * @param distance maximal core distance
     */
    void setCoreDistance(int distance);
    
    /**
     * @brief Set residue distance parameter
     * @param distance maximal residue assignment cut off
     */
    void setResidueDistance(double distance);
    
    /**
     * @brief Set inversion parameter
     * @param allow whether to allow inversion of SSEs
     */
    void setInversion(bool allow);
    
    /**
     * @brief Get core distance parameter
     * @return maximal core distance value
     */
    int getCoreDistance() const;
    
    /**
     * @brief Get residue distance parameter
     * @return maximal residue assignment cut off value
     */
    double getResidueDistance() const;
    
    /**
     * @brief Get inversion parameter
     * @return whether inversion of SSEs is allowed
     */
    bool getInversion() const;
    
private:
    int alCoreDistance;      ///< maximal core distance
    double alResidueDistance; ///< maximal residue assignment cut off
    bool alInversion;        ///< allow inversion of SSEs
};

} // namespace config
} // namespace gangsta

#endif // GANGSTA_CONFIG_H