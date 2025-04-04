#!/bin/bash

# Run benchmark with local PDB files
# This script tests both implementations with PDB files already in the repository

# Colors for better visibility
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Create results directory
mkdir -p benchmark_results

# List of PDB files to test
PDB_FILES=(
    "./examples/d1gkub1.pdb"
    "./examples/d2uaga1.pdb"
    "./build/edge_case_tests/normal1.pdb"
    "./build/edge_case_tests/normal2.pdb"
    "./build/edge_case_tests/normal.pdb"
    "./build/edge_case_tests/normal_struct.pdb"
)

# Check if edge case directory exists, if not, build tests
if [ ! -d "./build/edge_case_tests" ]; then
    echo -e "${CYAN}Building edge case tests...${NC}"
    cmake -B build
    cmake --build build
    ./build/run_tests --gtest_filter=EdgeCaseTest.*
fi

# Create result file header
echo "protein1,protein2,original_aligned,original_rmsd,original_score,modular_aligned,modular_rmsd,modular_score,match" > benchmark_results/comparison_results.csv

# Build modular implementation if needed
echo -e "\n${CYAN}Building modular implementation...${NC}"
cmake -B build
cmake --build build

# Check if original implementation exists, and if not, build it if possible
if [ -d "./original_code" ] && [ ! -f "./original_code/build/gplus" ]; then
    echo -e "\n${CYAN}Building original implementation...${NC}"
    mkdir -p ./original_code/build
    cd ./original_code/build
    cmake ..
    cmake --build .
    cd ../../
fi

# Run pairwise comparisons
echo -e "\n${CYAN}Running benchmarks...${NC}"
total_pairs=0
matching_results=0

# Test with all combinations of PDB files
for i in "${!PDB_FILES[@]}"; do
    protein1=${PDB_FILES[$i]}
    
    # Test against all other proteins
    for j in "${!PDB_FILES[@]}"; do
        # Skip comparing a protein with itself
        if [ $i -eq $j ]; then
            continue
        fi
        
        protein2=${PDB_FILES[$j]}
        protein1_name=$(basename $protein1 .pdb)
        protein2_name=$(basename $protein2 .pdb)
        
        total_pairs=$((total_pairs + 1))
        echo -e "\n${YELLOW}Testing ${protein1_name} vs ${protein2_name}${NC}"
        
        # Run modular implementation
        echo "Running modular implementation..."
        ./build/gplus ${protein1} ${protein2} --residue-distance 8.0 > modular_output.txt
        
        # Extract results
        if grep -q "Aligned length=" modular_output.txt; then
            modular_aligned=$(grep "Aligned length=" modular_output.txt | awk '{print $3}')
            modular_rmsd=$(grep "RMSD=" modular_output.txt | awk '{print $5}')
            modular_score=$(grep "GP-Score=" modular_output.txt | awk '{print $7}')
            echo "Modular results: Aligned=${modular_aligned}, RMSD=${modular_rmsd}, Score=${modular_score}"
        else
            echo "Error in modular implementation"
            cat modular_output.txt
            modular_aligned="ERROR"
            modular_rmsd="ERROR"
            modular_score="ERROR"
        fi
        
        # Check if original implementation is available
        if [ -f "./original_code/build/gplus" ]; then
            echo "Running original implementation..."
            
            # Copy PDB files to original_code/build directory for the original implementation
            protein1_basename=$(basename ${protein1})
            protein2_basename=$(basename ${protein2})
            cp ${protein1} ./original_code/build/${protein1_basename}
            cp ${protein2} ./original_code/build/${protein2_basename}
            
            # Change to original code directory and run
            cd ./original_code/build/
            ./gplus -r 8.0 ${protein1_basename} ${protein2_basename} > ../../original_output.txt
            cd ../../
            
            # Extract results
            if grep -q "Aligned length=" original_output.txt; then
                original_aligned=$(grep "Aligned length=" original_output.txt | awk '{print $3}')
                original_rmsd=$(grep "RMSD=" original_output.txt | awk '{print $5}')
                original_score=$(grep "GP-Score=" original_output.txt | awk '{print $7}')
                echo "Original results: Aligned=${original_aligned}, RMSD=${original_rmsd}, Score=${original_score}"
            else
                echo "Error in original implementation"
                cat original_output.txt
                original_aligned="ERROR"
                original_rmsd="ERROR"
                original_score="ERROR"
            fi
            
            # Check if results match
            if [ "$original_aligned" == "$modular_aligned" ] && \
               [ "$original_rmsd" == "$modular_rmsd" ] && \
               [ "$original_score" == "$modular_score" ]; then
                match="TRUE"
                matching_results=$((matching_results + 1))
            else
                match="FALSE"
            fi
            
            # Save results
            echo "${protein1_name},${protein2_name},${original_aligned},${original_rmsd},${original_score},${modular_aligned},${modular_rmsd},${modular_score},${match}" >> benchmark_results/comparison_results.csv
        else
            echo "Original implementation not found, skipping comparison"
            # Save only modular results
            echo "${protein1_name},${protein2_name},N/A,N/A,N/A,${modular_aligned},${modular_rmsd},${modular_score},N/A" >> benchmark_results/comparison_results.csv
        fi
    done
done

# Generate a summary report
echo -e "\n${CYAN}Benchmark Summary:${NC}"
echo "======================"
echo "Total protein pairs tested: ${total_pairs}"
if [ -f "./original_code/build/gplus" ]; then
    echo "Results matching original implementation: ${matching_results}/${total_pairs}"
    match_percentage=$((matching_results * 100 / total_pairs))
    echo "Match percentage: ${match_percentage}%"
else
    echo "Original implementation not found, no comparison was made"
    echo "Generated results with modular implementation:"
    
    # Format the benchmark results as a table
    echo -e "\n${CYAN}Results:${NC}"
    echo "---------------------------------------------------"
    echo "| Protein 1      | Protein 2      | Aligned | RMSD   | Score    |"
    echo "|----------------|----------------|---------|--------|----------|"
    
    # Skip header and process each line
    tail -n +2 benchmark_results/comparison_results.csv | while IFS=, read -r p1 p2 o_aligned o_rmsd o_score m_aligned m_rmsd m_score match; do
        printf "| %-14s | %-14s | %7s | %6s | %8s |\n" "$p1" "$p2" "$m_aligned" "$m_rmsd" "$m_score"
    done
    
    echo "---------------------------------------------------"
fi
echo "Results saved to: benchmark_results/comparison_results.csv"

# Clean up
rm -f modular_output.txt original_output.txt