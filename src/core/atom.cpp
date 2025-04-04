#include "core/atom.h"
#include <unordered_map>
#include <algorithm>

namespace gangsta {
namespace core {

// Static mapping of 3-letter to 1-letter amino acid codes
static const std::unordered_map<std::string, char> RESIDUE_CODES = {
    {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'},
    {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'},
    {"HIS", 'H'}, {"ILE", 'I'}, {"LEU", 'L'}, {"LYS", 'K'},
    {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'}, {"SER", 'S'},
    {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'},
    // Non-standard amino acids
    {"MSE", 'M'}, {"HSD", 'H'}, {"HSE", 'H'}, {"HSP", 'H'},
    {"HID", 'H'}, {"HIE", 'H'}, {"HIP", 'H'}, {"UNK", 'X'}
};

Atom::Atom() 
    : position(0.0, 0.0, 0.0), 
      atomName(""), 
      residueName(""), 
      residueNumber(0), 
      chainID(' '), 
      sseIndex(-1) 
{
}

Atom::Atom(const math::Vec3& position, 
           const std::string& atomName, 
           const std::string& residueName, 
           int residueNumber,
           char chainID) 
    : position(position), 
      atomName(atomName), 
      residueName(residueName), 
      residueNumber(residueNumber), 
      chainID(chainID), 
      sseIndex(-1) 
{
}

const math::Vec3& Atom::getPosition() const 
{
    return position;
}

void Atom::setPosition(const math::Vec3& position) 
{
    this->position = position;
}

const std::string& Atom::getAtomName() const 
{
    return atomName;
}

void Atom::setAtomName(const std::string& name) 
{
    atomName = name;
}

const std::string& Atom::getResidueName() const 
{
    return residueName;
}

void Atom::setResidueName(const std::string& name) 
{
    residueName = name;
}

int Atom::getResidueNumber() const 
{
    return residueNumber;
}

void Atom::setResidueNumber(int number) 
{
    residueNumber = number;
}

char Atom::getChainID() const 
{
    return chainID;
}

void Atom::setChainID(char chainID) 
{
    this->chainID = chainID;
}

char Atom::getResidueCode() const 
{
    return convertResidueCode(residueName);
}

bool Atom::isAlphaCarbon() const 
{
    return (atomName == "CA");
}

int Atom::getSSEIndex() const 
{
    return sseIndex;
}

void Atom::setSSEIndex(int index) 
{
    sseIndex = index;
}

double Atom::distanceTo(const Atom& other) const 
{
    return position.dist(other.position);
}

char Atom::convertResidueCode(const std::string& residueName) 
{
    // Convert to uppercase for lookup
    std::string upperRes = residueName;
    std::transform(upperRes.begin(), upperRes.end(), upperRes.begin(), ::toupper);
    
    // Look up the code
    auto it = RESIDUE_CODES.find(upperRes);
    if (it != RESIDUE_CODES.end()) {
        return it->second;
    }
    
    // Return 'X' for unknown residues
    return 'X';
}

} // namespace core
} // namespace gangsta