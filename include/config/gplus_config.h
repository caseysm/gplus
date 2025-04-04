#ifndef GANGSTA_GPLUS_CONFIG_H
#define GANGSTA_GPLUS_CONFIG_H

namespace gangsta {
namespace config {

/**
 * @brief Advanced configuration parameters for the GANGSTA+ algorithm
 */
class GPlusConfig 
{
public:
    /**
     * @brief Constructor with default values
     */
    GPlusConfig();
    
    // Getters
    int getEvaluationDepth() const;
    int getResultCount() const;
    int getCoreDelta() const;
    double getDistMax() const;
    double getRescale() const;
    
    // Setters
    void setEvaluationDepth(int depth);
    void setResultCount(int count);
    void setCoreDelta(int delta);
    void setDistMax(double max);
    void setRescale(double scale);
    
private:
    // General configuration
    int coResults;      ///< number of results to be refined
    
    // Info
    double inDistMax;   ///< maximum distance for info object
    double inRescale;   ///< rescale factor for info object
    
    // Prepare
    int ppCoreDelta;    ///< maximal SSE length difference
    
    // Evaluation
    int evDepth;        ///< evaluation depth
};

} // namespace config
} // namespace gangsta

#endif // GANGSTA_GPLUS_CONFIG_H