#ifndef GANGSTA_RANDOM_H
#define GANGSTA_RANDOM_H

#include <random>
#include <algorithm>
#include <vector>
#include <limits>

namespace gangsta {
namespace utils {

/**
 * @brief Random number generator utility
 */
class Random 
{
public:
    /**
     * @brief Get the singleton instance
     * @return Reference to the random generator
     */
    static Random& getInstance();
    
    /**
     * @brief Initialize with a specific seed
     * @param seed Random seed
     */
    void initialize(unsigned int seed);
    
    /**
     * @brief Initialize with a random seed
     */
    void initializeRandom();
    
    /**
     * @brief Get a random integer in a range
     * @param min Minimum value (inclusive)
     * @param max Maximum value (inclusive)
     * @return Random integer in the range [min, max]
     */
    int getInt(int min, int max);
    
    /**
     * @brief Get a random float in a range
     * @param min Minimum value (inclusive)
     * @param max Maximum value (inclusive)
     * @return Random float in the range [min, max]
     */
    float getFloat(float min, float max);
    
    /**
     * @brief Get a random double in a range
     * @param min Minimum value (inclusive)
     * @param max Maximum value (inclusive)
     * @return Random double in the range [min, max]
     */
    double getDouble(double min, double max);
    
    /**
     * @brief Get a random boolean
     * @param probability Probability of true (0.0-1.0)
     * @return Random boolean
     */
    bool getBool(double probability = 0.5);
    
    /**
     * @brief Get a random sign (-1 or 1)
     * @return Random sign
     */
    int getSign();
    
    /**
     * @brief Randomly shuffle a vector
     * @param vec Vector to shuffle
     */
    template<typename T>
    void shuffle(std::vector<T>& vec) {
        std::shuffle(vec.begin(), vec.end(), engine);
    }
    
private:
    /**
     * @brief Private constructor (singleton)
     */
    Random();
    
    /**
     * @brief Private destructor
     */
    ~Random();
    
    std::mt19937 engine;        ///< Mersenne Twister engine
    bool initialized;           ///< Whether the generator is initialized
    
    // Make the generator a singleton
    Random(const Random&) = delete;
    Random& operator=(const Random&) = delete;
};

} // namespace utils
} // namespace gangsta

#endif // GANGSTA_RANDOM_H