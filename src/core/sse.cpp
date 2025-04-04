#include "core/sse.h"

namespace gangsta {
namespace core {

SSE::SSE() 
    : type(UNKNOWN), startResidueIndex(0), endResidueIndex(0) 
{
}

SSE::SSE(Type type, int startResidueIndex, int endResidueIndex) 
    : type(type), startResidueIndex(startResidueIndex), endResidueIndex(endResidueIndex) 
{
}

SSE::Type SSE::getType() const 
{
    return type;
}

void SSE::setType(Type type) 
{
    this->type = type;
}

int SSE::getStart() const 
{
    return startResidueIndex;
}

void SSE::setStart(int index) 
{
    startResidueIndex = index;
}

int SSE::getEnd() const 
{
    return endResidueIndex;
}

void SSE::setEnd(int index) 
{
    endResidueIndex = index;
}

int SSE::getLength() const 
{
    return endResidueIndex - startResidueIndex + 1;
}

std::string SSE::getTypeName() const 
{
    switch (type) {
        case HELIX:
            return "Helix";
        case STRAND:
            return "Strand";
        case LOOP:
            return "Loop";
        case UNKNOWN:
        default:
            return "Unknown";
    }
}

char SSE::getTypeCode() const 
{
    switch (type) {
        case HELIX:
            return 'H';
        case STRAND:
            return 'E';
        case LOOP:
            return 'L';
        case UNKNOWN:
        default:
            return '-';
    }
}

bool SSE::containsResidue(int residueIndex) const 
{
    return (residueIndex >= startResidueIndex && residueIndex <= endResidueIndex);
}

} // namespace core
} // namespace gangsta