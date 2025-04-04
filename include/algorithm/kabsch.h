#ifndef GANGSTA_KABSCH_H
#define GANGSTA_KABSCH_H

#include "math/vec3.h"
#include "math/rotation.h"
#include <vector>

namespace gangsta {
namespace algorithm {

/**
 * @brief Implementation of the Kabsch algorithm for optimal alignment of two point sets
 * 
 * The Kabsch algorithm finds the optimal rotation matrix that minimizes the
 * root mean squared deviation (RMSD) between two paired sets of points.
 */
class Kabsch 
{
public:
    /**
     * @brief Default constructor
     */
    Kabsch();
    
    /**
     * @brief Calculate the optimal rotation and translation to align two point sets
     * @param points1 First set of points
     * @param points2 Second set of points (must be the same size as points1)
     * @param rotation Output rotation
     * @param translation Output translation
     * @return RMSD after alignment
     */
    double align(const std::vector<math::Vec3>& points1, 
                const std::vector<math::Vec3>& points2,
                math::Rotation& rotation,
                math::Vec3& translation);
    
    /**
     * @brief Calculate RMSD between two point sets without aligning
     * @param points1 First set of points
     * @param points2 Second set of points (must be the same size as points1)
     * @return RMSD between the point sets
     */
    double calculateRMSD(const std::vector<math::Vec3>& points1, 
                        const std::vector<math::Vec3>& points2);
    
    /**
     * @brief Calculate the centroid of a set of points
     * @param points Set of points
     * @return Centroid position
     */
    math::Vec3 calculateCentroid(const std::vector<math::Vec3>& points);
    
private:
    /**
     * @brief Compute the covariance matrix between two point sets
     * @param points1 First set of points (centered at origin)
     * @param points2 Second set of points (centered at origin)
     * @param covariance Output 3x3 covariance matrix (as 9-element array)
     */
    void computeCovarianceMatrix(const std::vector<math::Vec3>& points1,
                                const std::vector<math::Vec3>& points2,
                                double covariance[9]);
    
    /**
     * @brief Compute the optimal rotation matrix using singular value decomposition
     * @param covariance 3x3 covariance matrix (as 9-element array)
     * @param rotation Output rotation matrix
     */
    void computeOptimalRotation(const double covariance[9], math::Rotation& rotation);
};

} // namespace algorithm
} // namespace gangsta

#endif // GANGSTA_KABSCH_H