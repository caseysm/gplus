#include "utils/random.h"
#include <ctime>
#include <algorithm>
#include <vector>

namespace gangsta {
namespace utils {

Random::Random() 
    : initialized(false) 
{
    initializeRandom();
}

Random::~Random() 
{
}

Random& Random::getInstance() 
{
    static Random instance;
    return instance;
}

void Random::initialize(unsigned int seed) 
{
    engine.seed(seed);
    initialized = true;
}

void Random::initializeRandom() 
{
    // Use a combination of time and address of object for better randomness
    unsigned int seed = static_cast<unsigned int>(
        std::time(nullptr) ^ 
        reinterpret_cast<uintptr_t>(this)
    );
    
    initialize(seed);
}

int Random::getInt(int min, int max) 
{
    std::uniform_int_distribution<int> dist(min, max);
    return dist(engine);
}

float Random::getFloat(float min, float max) 
{
    std::uniform_real_distribution<float> dist(min, max);
    return dist(engine);
}

double Random::getDouble(double min, double max) 
{
    std::uniform_real_distribution<double> dist(min, max);
    return dist(engine);
}

bool Random::getBool(double probability) 
{
    std::bernoulli_distribution dist(probability);
    return dist(engine);
}

int Random::getSign() 
{
    return getBool() ? 1 : -1;
}

} // namespace utils
} // namespace gangsta