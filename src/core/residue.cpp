#include "core/residue.h"
#include <algorithm>
#include <unordered_map>

namespace gangsta {
namespace core {

Residue::Residue() 
    : name("UNK"), number(0), chainID(' '), sseIndex(-1) 
{
}

Residue::Residue(const std::string& name, int number, char chainID) 
    : name(name), number(number), chainID(chainID), sseIndex(-1) 
{
}

void Residue::addAtom(const Atom& atom) 
{
    atoms.push_back(atom);
}

const std::vector<Atom>& Residue::getAtoms() const 
{
    return atoms;
}

Atom* Residue::getAtomByName(const std::string& name) 
{
    for (auto& atom : atoms) {
        if (atom.getAtomName() == name) {
            return &atom;
        }
    }
    return nullptr;
}

const Atom* Residue::getAtomByName(const std::string& name) const 
{
    for (const auto& atom : atoms) {
        if (atom.getAtomName() == name) {
            return &atom;
        }
    }
    return nullptr;
}

Atom* Residue::getAlphaCarbon() 
{
    return getAtomByName("CA");
}

const Atom* Residue::getAlphaCarbon() const 
{
    return getAtomByName("CA");
}

const std::string& Residue::getName() const 
{
    return name;
}

char Residue::getCode() const 
{
    // Use the residue code function from any atom (they should all be the same)
    if (!atoms.empty()) {
        return atoms[0].getResidueCode();
    }
    
    // Fallback case if there are no atoms
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
    
    auto it = RESIDUE_CODES.find(name);
    if (it != RESIDUE_CODES.end()) {
        return it->second;
    }
    
    return 'X'; // Unknown residue
}

int Residue::getNumber() const 
{
    return number;
}

char Residue::getChainID() const 
{
    return chainID;
}

void Residue::setSSEIndex(int index) 
{
    sseIndex = index;
    
    // Update all atoms in this residue
    for (auto& atom : atoms) {
        atom.setSSEIndex(index);
    }
}

int Residue::getSSEIndex() const 
{
    return sseIndex;
}

math::Vec3 Residue::getCenter() const 
{
    if (atoms.empty()) {
        return math::Vec3(0, 0, 0);
    }
    
    math::Vec3 center(0, 0, 0);
    for (const auto& atom : atoms) {
        center = center + atom.getPosition();
    }
    
    return center / atoms.size();
}

} // namespace core
} // namespace gangsta