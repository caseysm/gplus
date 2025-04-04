#include "config/gplus_config.h"

namespace gangsta {
namespace config {

GPlusConfig::GPlusConfig() 
{
    evDepth = 5000;     // evaluation depth
    coResults = 200;    // number of results to be refined
    
    // prepare
    ppCoreDelta = 7;    // maximal SSE length difference
    
    // info object parameters
    inDistMax = 11;
    inRescale = 5;
}

// Getters
int GPlusConfig::getEvaluationDepth() const
{
    return evDepth;
}

int GPlusConfig::getResultCount() const
{
    return coResults;
}

int GPlusConfig::getCoreDelta() const
{
    return ppCoreDelta;
}

double GPlusConfig::getDistMax() const
{
    return inDistMax;
}

double GPlusConfig::getRescale() const
{
    return inRescale;
}

// Setters
void GPlusConfig::setEvaluationDepth(int depth)
{
    evDepth = depth;
}

void GPlusConfig::setResultCount(int count)
{
    coResults = count;
}

void GPlusConfig::setCoreDelta(int delta)
{
    ppCoreDelta = delta;
}

void GPlusConfig::setDistMax(double max)
{
    inDistMax = max;
}

void GPlusConfig::setRescale(double scale)
{
    inRescale = scale;
}

} // namespace config
} // namespace gangsta