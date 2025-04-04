#!/bin/bash

# Build the project with tests
cmake -B build
cmake --build build

# Create test output directories if they don't exist
mkdir -p build/test_output
mkdir -p build/benchmark_results
mkdir -p build/edge_case_tests
mkdir -p build/comparison_tests

# Define color codes for better visibility
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Function to run a specific test and report results
run_test_suite() {
    local filter=$1
    local description=$2
    
    echo -e "\n${CYAN}Running $description tests${NC}"
    
    # Run the tests with the specified filter
    pushd build >/dev/null
    ./run_tests --gtest_filter="$filter.*" --gtest_brief=1
    local result=$?
    popd >/dev/null
    
    # Check if tests passed
    if [ $result -eq 0 ]; then
        echo -e "${GREEN}✓ All $description tests passed${NC}"
    else
        echo -e "${YELLOW}⚠ Some $description tests failed${NC}"
    fi
    
    return $result
}

# Run basic unit tests
echo -e "${CYAN}Running unit tests...${NC}"
pushd build >/dev/null
./run_tests --gtest_filter="-IntegrationTest.*:BenchmarkTest.*:EdgeCaseTest.*" --gtest_brief=1
basic_test_result=$?
popd >/dev/null

# Run integration tests
echo -e "\n${CYAN}Running integration tests...${NC}"
pushd build >/dev/null
run_test_suite "IntegrationTest" "integration"
integration_test_result=$?
popd >/dev/null

# Run edge case tests
echo -e "\n${CYAN}Running edge case tests...${NC}"
pushd build >/dev/null
run_test_suite "EdgeCaseTest" "edge case"
edge_test_result=$?
popd >/dev/null

# Run benchmarks (if --with-benchmarks is provided)
benchmark_test_result=0
if [ "$1" == "--with-benchmarks" ] || [ "$2" == "--with-benchmarks" ]; then
    echo -e "\n${CYAN}Running benchmarks...${NC}"
    pushd build >/dev/null
    run_test_suite "BenchmarkTest" "benchmark"
    benchmark_test_result=$?
    
    echo -e "\n${CYAN}Benchmark results saved to:${NC} build/benchmark_results/benchmark_results.csv"
    popd >/dev/null
else
    echo -e "\n${YELLOW}Benchmarks skipped. Use --with-benchmarks to run them.${NC}"
fi

# Run comparison tests (if --with-comparison is provided)
comparison_test_result=0
if [ "$1" == "--with-comparison" ] || [ "$2" == "--with-comparison" ]; then
    echo -e "\n${CYAN}Running original code comparison tests...${NC}"
    pushd build >/dev/null
    run_test_suite "ComparisonTest" "comparison"
    comparison_test_result=$?
    popd >/dev/null
else
    echo -e "\n${YELLOW}Comparison tests skipped. Use --with-comparison to run them.${NC}"
fi

# Calculate total tests
total_unit_tests=71
total_integration_tests=$(grep -c "TEST_F(IntegrationTest," tests/integration_tests.cpp)
total_edge_tests=$(grep -c "TEST_F(EdgeCaseTest," tests/edge_cases.cpp)
total_benchmark_tests=$(grep -c "TEST_F(BenchmarkTest," tests/benchmarks.cpp)
total_comparison_tests=$(grep -c "TEST_F(ComparisonTest," tests/comparison_tests.cpp)
total_tests=$((total_unit_tests + total_integration_tests + total_edge_tests + total_benchmark_tests + total_comparison_tests))

# Print summary
echo ""
echo -e "${CYAN}Test summary:${NC}"
echo "============="
echo -e "Total tests: ${total_tests} tests in 13 test suites"
echo ""
echo "Basic unit tests: ($total_unit_tests tests)"
echo "- Math module tests: (26 tests)"
echo "  - Vec3Test: 15 tests for vector operations"
echo "  - Mat4Test: 11 tests for matrix operations"
echo ""
echo "- Core module tests: (23 tests)"
echo "  - AtomTest: 11 tests for atom functionality"
echo "  - ResidueTest: 12 tests for residue management"
echo ""
echo "- Algorithm module tests: (16 tests)"
echo "  - KabschTest: 6 tests for RMSD calculations"
echo "  - AlignmentTest: 10 tests for alignment algorithms"
echo ""
echo "- Utils module tests: (6 tests)"
echo "  - PDBParserTest: 6 tests for PDB file parsing and generation"
echo ""
echo "Advanced tests: ($((total_integration_tests + total_edge_tests + total_benchmark_tests + total_comparison_tests)) tests)"
echo "- Integration tests: $total_integration_tests end-to-end tests"
echo "- Edge case tests: $total_edge_tests tests for boundary conditions"
echo "- Benchmark tests: $total_benchmark_tests performance tests"
echo "- Comparison tests: $total_comparison_tests original vs. modular tests"
echo ""
echo "Test commands:"
echo "- Run all tests: ./run_tests.sh"
echo "- Run with benchmarks: ./run_tests.sh --with-benchmarks"
echo "- Run with comparison: ./run_tests.sh --with-comparison"
echo "- Run with both: ./run_tests.sh --with-benchmarks --with-comparison"
echo "- Run specific test: ./build/run_tests --gtest_filter=<TestName>.*"
echo "  Example: ./build/run_tests --gtest_filter=Vec3Test.*"
echo "- Run integration tests only: ./build/run_tests --gtest_filter=IntegrationTest.*"
echo "- Run edge case tests only: ./build/run_tests --gtest_filter=EdgeCaseTest.*"
echo "- Run benchmark tests only: ./build/run_tests --gtest_filter=BenchmarkTest.*"
echo "- Run comparison tests only: ./build/run_tests --gtest_filter=ComparisonTest.*"
echo ""
echo "For detailed test output: ctest -VV"