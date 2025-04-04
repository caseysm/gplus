#include "config/config.h"

namespace gangsta {
namespace config {

Config::Config() 
{
    // Default alignment detection parameters
    alCoreDistance = 4;       // maximal core distance
    alResidueDistance = 6;    // maximal residue assignment cut off
    alInversion = false;      // allow inversion of SSEs
}

void Config::setCoreDistance(int distance)
{
    alCoreDistance = distance;
}

void Config::setResidueDistance(double distance)
{
    alResidueDistance = distance;
}

void Config::setInversion(bool allow)
{
    alInversion = allow;
}

int Config::getCoreDistance() const
{
    return alCoreDistance;
}

double Config::getResidueDistance() const
{
    return alResidueDistance;
}

bool Config::getInversion() const
{
    return alInversion;
}

} // namespace config
} // namespace gangsta