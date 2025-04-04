# GANGSTA+ Modularization Plan

## Overview
This document outlines a plan to restructure the monolithic `gplus.cpp` file into a modular architecture with separate header and implementation files.

## Directory Structure
```
gplus/
├── include/                  # Header files
│   ├── config/               # Configuration settings
│   ├── math/                 # Mathematical utilities 
│   ├── core/                 # Core protein structures
│   ├── algorithm/            # Alignment algorithms
│   └── utils/                # General utilities
├── src/                      # Implementation files
│   ├── config/
│   ├── math/
│   ├── core/
│   ├── algorithm/
│   └── utils/
├── examples/                 # Example PDB files
├── tests/                    # Unit tests
└── apps/                     # Applications (main executable)
```

## Module Breakdown

### 1. Configuration (config/)
- **config.h/cpp**: Basic alignment parameters
- **gplus_config.h/cpp**: Algorithm-specific parameters

### 2. Math Utilities (math/)
- **vec3.h/cpp**: 3D vector operations
- **vec4.h/cpp**: 4D vector operations (quaternions)
- **mat4.h/cpp**: 4x4 matrix operations
- **rotation.h/cpp**: Rotation utilities
- **transformation.h/cpp**: Coordinate transformations
- **kabsch.h/cpp**: Kabsch algorithm for RMSD calculation

### 3. Core Structures (core/)
- **atom.h/cpp**: Atom representation
- **residue.h/cpp**: Amino acid residue representation
- **sse.h/cpp**: Secondary structure elements
- **molecule.h/cpp**: Complete molecular structure
- **forceset.h/cpp**: Forces and energy calculation

### 4. Algorithms (algorithm/)
- **alignment.h/cpp**: Core alignment algorithm
- **scoring.h/cpp**: Alignment scoring functions
- **search.h/cpp**: Structure search algorithms

### 5. Utilities (utils/)
- **file.h/cpp**: File I/O operations
- **logging.h/cpp**: Logging and messaging
- **random.h/cpp**: Random number generation
- **string_utils.h/cpp**: String manipulation
- **pymol.h/cpp**: PyMol script generation

### 6. Application (apps/)
- **gplus.cpp**: Main application

## Implementation Strategy

### Phase 1: Initial Modularization
1. Create directory structure
2. Extract configuration classes to separate files
3. Move mathematical utilities to math module
4. Extract core structures to core module

### Phase 2: Algorithm Separation
1. Separate alignment algorithms
2. Move utility functions to appropriate modules
3. Establish proper header dependencies

### Phase 3: Main Application
1. Modify main application to use modular components
2. Update CMake build system to build from modules

### Phase 4: Testing
1. Create unit tests for individual modules
2. Ensure alignment results match original implementation
3. Validate on example PDB files

### Phase 5: Algorithm Validation and Alignment
1. Implement debug mode for algorithm tracing in modular code
2. Instrument original code at key decision points
3. Create side-by-side comparison framework for detailed analysis
4. Identify and document all algorithm discrepancies
5. Update modular implementation to match original behavior
6. Develop comprehensive regression test suite

### Validation Strategy
1. **Debug Instrumentation**
   - Add detailed logging to both implementations
   - Create output tracing for SSE detection, pairing, and alignment steps
   - Log all parameter values and intermediate calculations

2. **Systematic Comparison**
   - Create test fixtures with diverse protein structures
   - Compare algorithm outputs step-by-step for each test case
   - Identify exact points of divergence between implementations

3. **Parameter Verification**
   - Create tests to verify parameter equivalence
   - Ensure default values match between implementations
   - Document any parameter differences that affect output

4. **Modular Code Adaptation**
   - Update SSE detection and matching logic if needed
   - Fix residue pairing algorithms to match original
   - Align transformation and RMSD calculation methods
   - Update scoring function to match original algorithm

5. **Regression Testing**
   - Develop automated tests comparing outputs against original
   - Create benchmarks for performance comparison
   - Test edge cases to ensure consistent behavior

## Code Guidelines
- Use include guards or pragma once in all headers
- Minimize dependencies between modules
- Use forward declarations where appropriate
- Define clear interfaces for each module
- Maintain existing performance optimizations
- Add appropriate namespaces

## Benefits
- Improved code maintainability
- Easier understanding of the codebase
- Better separation of concerns
- Potential for reusing components in other applications
- Easier to extend with new features