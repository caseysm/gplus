#!/bin/bash

# This script adds instrumentation to the original GANGSTA+ implementation
# to trace algorithm execution and to make it easier to compare with
# the modular implementation.

# Ensure script is run from the correct directory
if [ ! -f "gplus.cpp" ]; then
    echo "Error: This script must be run from the original_code directory"
    exit 1
fi

# Make a backup of the original file
cp gplus.cpp gplus.cpp.bak

# Add debug output macros at the beginning of the file
sed -i '1s/^/#define DEBUG_MODE\n\n#include <fstream>\n/' gplus.cpp

# Add debug output function
sed -i '/^using namespace std;/a \
\
// Debug output function\
void debug_log(const std::string& message) {\
    static std::ofstream debug_file("debug_original.txt", std::ios::out | std::ios::app);\
    debug_file << "[DEBUG] " << message << std::endl;\
    debug_file.flush();\
}' gplus.cpp

# Add instrumentation to SSE detection
sed -i 's/void construct (SpecMolecule\* target, string idb)/void construct (SpecMolecule\* target, string idb)\n{\n#ifdef DEBUG_MODE\n    debug_log("Starting alignment");\n    debug_log("Molecule 1 size: " + to_string(target->lcalpha.size()));\n#endif/;s/return construct(target, idb);/return construct(target, idb);\n}/;s/double construct (SpecMolecule\* target, SpecMolecule\* mol)/double construct (SpecMolecule\* target, SpecMolecule\* mol)\n{/' gplus.cpp

# Add instrumentation to the SSE alignment stage
sed -i '/topo.build/a\\n#ifdef DEBUG_MODE\n    debug_log("Building topology");\n    debug_log("Target SSE count: " + to_string(target->lsse.size()));\n    debug_log("Mol SSE count: " + to_string(mol->lsse.size()));\n#endif' gplus.cpp

# Add instrumentation to key alignment steps
sed -i '/Aligned length=  /a\\n#ifdef DEBUG_MODE\n    debug_log("Alignment result:");\n    debug_log("  Aligned length: " + to_string(pairs));\n    debug_log("  RMSD: " + to_string(rmsd));\n    debug_log("  GP-Score: " + to_string(score));\n#endif' gplus.cpp

# Add instrumentation to the alignment result
sed -i '/cout << "GANGSTA+ Done.";/a\\n#ifdef DEBUG_MODE\n    debug_log("Alignment completed");\n#endif' gplus.cpp

echo "Instrumentation added to gplus.cpp"
echo "Original file backed up as gplus.cpp.bak"
echo ""
echo "Build the instrumented version with:"
echo "  mkdir -p build && cd build && cmake .. && make"
echo ""
echo "Run with:"
echo "  ./gplus d1gkub1.pdb d2uaga1.pdb"
echo ""
echo "Debug output will be written to debug_original.txt"