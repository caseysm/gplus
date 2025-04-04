# GANGSTA+ Protein Structure Alignment

GANGSTA+ (GANGSTA Non-Sequential Search and Alignment) is a non-sequential protein structure alignment algorithm that considers the arrangement of secondary structure elements (SSEs) for structure comparison.

## Features

- Sequential and non-sequential protein structure alignment
- Secondary structure element (SSE) based alignment
- Integrated PDB file parser and writer
- PyMOL output for visualization
- Customizable alignment parameters

## Building

```bash
# Configure and build the project
cmake -B build
cmake --build build

# Run the executable
./build/gplus examples/d2uaga1.pdb examples/d1gkub1.pdb
```

## Usage

```
./build/gplus [options] <reference_pdb> <target_pdb>

Options:
  -s, --sequential      Use sequential alignment instead of non-sequential
  -o, --output <prefix> Output prefix for result files
  -c, --core-dist <d>   Maximum core distance (default: 4.0)
  -r, --res-dist <d>    Maximum residue distance (default: 8.0)
  -i, --iterations <n>  Maximum optimization iterations (default: 10)
  -h, --help            Display this help message
```

## Example

```bash
# Run a non-sequential alignment
./build/gplus -o result examples/d2uaga1.pdb examples/d1gkub1.pdb

# Run a sequential alignment
./build/gplus -s -o result examples/d2uaga1.pdb examples/d1gkub1.pdb
```

## Testing

The project includes a comprehensive testing framework with unit tests, integration tests, edge case tests, and benchmarks using Google Test:

```bash
# Run all tests (excluding benchmarks)
./run_tests.sh

# Run all tests including benchmarks
./run_tests.sh --with-benchmarks

# Run specific test categories
./build/run_tests --gtest_filter=Vec3Test.*        # Individual unit tests
./build/run_tests --gtest_filter=IntegrationTest.* # Integration tests
./build/run_tests --gtest_filter=EdgeCaseTest.*    # Edge case tests
./build/run_tests --gtest_filter=BenchmarkTest.*   # Performance benchmarks

# Run with detailed output
ctest -VV
```

### Testing Framework

The testing framework is divided into several components:

1. **Unit Tests**
   - 71 basic unit tests covering all components
   - Verifies correctness of individual functions and classes
   - Ensures each component works correctly in isolation

2. **Integration Tests**
   - End-to-end workflow tests for complete algorithms
   - Tests interaction between multiple components
   - Validates sequential and non-sequential alignment algorithms
   - Verifies output generation (PDB files, PyMOL scripts)
   - Tests alignment consistency across multiple runs

3. **Edge Case Tests**
   - Tests for boundary conditions and unusual inputs
   - Validates behavior with empty or minimal PDB files
   - Tests extreme parameter values
   - Verifies numerical stability with large coordinate values
   - Handles missing or invalid data gracefully

4. **Performance Benchmarks**
   - Measures execution time for critical operations
   - Compares performance with different configuration parameters
   - Records benchmark results in CSV format for tracking
   - Establishes baseline metrics for future optimization
   - Statistical analysis of multiple runs for reliable measurements

### Test Coverage

The project includes a comprehensive test suite with over 90 test cases:

- **Math Module (26 tests)**
  - Vec3Test: 15 tests for vector operations
  - Mat4Test: 11 tests for matrix operations

- **Core Module (23 tests)**
  - AtomTest: 11 tests for atom functionality
  - ResidueTest: 12 tests for residue management

- **Algorithm Module (16 tests)**
  - KabschTest: 6 tests for RMSD calculations
  - AlignmentTest: 10 tests for alignment algorithms

- **Utils Module (6 tests)**
  - PDBParserTest: 6 tests for PDB file parsing and generation

- **Advanced Tests**
  - IntegrationTest: 4 end-to-end workflow tests
  - EdgeCaseTest: 7 boundary condition tests
  - BenchmarkTest: 4 performance measurement tests

### Test Structure

The test suite is organized into modules that mirror the project structure:

- `tests/math/`: Tests for vector, matrix, and rotation classes
- `tests/core/`: Tests for atom, residue, molecule, and SSE classes
- `tests/algorithm/`: Tests for alignment and Kabsch algorithms
- `tests/utils/`: Tests for file I/O, PDB parsing, and PyMOL generation
- `tests/integration_tests.cpp`: End-to-end workflow tests
- `tests/edge_cases.cpp`: Boundary condition tests
- `tests/benchmarks.cpp`: Performance measurement tests

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

The original GANGSTA+ algorithm was developed by:
- M. Guerler, A. and Knapp, E.W. GEN: a Conformational Database Filter for Ligand-Bound Protein Structures. In: J Comput Chem 29.12 (2008), pp. 2012-2022.
- A. Guerler, and E.W. Knapp. Novel Protein Folds and Their Non-sequential Structural Analogs. In: Protein Science 17.8 (2008), pp. 1374-1382.