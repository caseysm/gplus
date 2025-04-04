#!/bin/bash

# GANGSTA+ Benchmark Script
# This script runs benchmarks for the GANGSTA+ structural alignment tool
# It compares the modularized implementation with various parameters

# Configuration
OUTPUT_DIR="/home/casey/Desktop/repos/gplus/tests/data/results"
GPLUS_PATH="/home/casey/Desktop/repos/gplus/build/gplus"
PDB_DIR="/home/casey/Desktop/repos/gplus/tests/data/benchmark/pdb"

# Create directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$PDB_DIR"

# Check if the TM-align dataset has been downloaded
if [ ! -f "$PDB_DIR/list.txt" ]; then
    echo "Downloading TM-align benchmark dataset..."
    mkdir -p "$PDB_DIR/tmp"
    cd "$PDB_DIR/tmp"
    
    # Download a subset of the TM-align benchmark dataset
    wget -q https://zhanggroup.org/TM-align/benchmark/raft1.tar.gz
    
    # Extract the dataset
    tar -xzf raft1.tar.gz
    
    # Move PDB files to the benchmark directory
    find . -name "*.pdb" -exec cp {} "$PDB_DIR" \;
    
    # Create a list of all PDB files
    find "$PDB_DIR" -name "*.pdb" | sort > "$PDB_DIR/list.txt"
    
    # Clean up
    cd -
    rm -rf "$PDB_DIR/tmp"
    
    echo "Download complete. PDB files are stored in $PDB_DIR"
fi

# Create CSV file for results
RESULTS_FILE="$OUTPUT_DIR/tm_benchmark_results.csv"
echo "PDB1,PDB2,Mode,CoreDistance,ResidueDistance,AlignedCount,RMSD,Score,Time" > "$RESULTS_FILE"

# Function to run benchmark for a pair of PDB files
run_benchmark() {
    local pdb1="$1"
    local pdb2="$2"
    local mode="$3"
    local core_dist="$4"
    local res_dist="$5"
    
    local pdb1_name=$(basename "$pdb1")
    local pdb2_name=$(basename "$pdb2")
    
    # Build the command
    local cmd="$GPLUS_PATH"
    if [ "$mode" == "Sequential" ]; then
        cmd="$cmd --sequential"
    fi
    cmd="$cmd --core-distance $core_dist --residue-distance $res_dist $pdb1 $pdb2"
    
    # Run the command and time it
    echo "Running benchmark: $cmd"
    local start_time=$(date +%s.%N)
    local output=$($cmd)
    local end_time=$(date +%s.%N)
    local time_diff=$(echo "$end_time - $start_time" | bc)
    
    # Parse the results
    local aligned_count=$(echo "$output" | grep "Aligned length" | awk '{print $3}')
    local rmsd=$(echo "$output" | grep "RMSD" | awk '{print $2}')
    local score=$(echo "$output" | grep "GP-Score" | awk '{print $2}')
    
    # Write results to CSV
    echo "$pdb1_name,$pdb2_name,$mode,$core_dist,$res_dist,$aligned_count,$rmsd,$score,$time_diff" >> "$RESULTS_FILE"
    
    echo "  - Results: Aligned=$aligned_count, RMSD=$rmsd, Score=$score, Time=$time_diff"
}

# Get the list of PDB files
pdb_files=($(cat "$PDB_DIR/list.txt"))

# If there are not enough PDB files, copy the existing ones to use as test cases
if [ ${#pdb_files[@]} -lt 5 ]; then
    cp "/home/casey/Desktop/repos/gplus/d1gkub1.pdb" "$PDB_DIR/"
    cp "/home/casey/Desktop/repos/gplus/d2uaga1.pdb" "$PDB_DIR/"
    pdb_files=($(find "$PDB_DIR" -name "*.pdb" | sort))
    echo "Using local test PDB files: ${pdb_files[@]}"
fi

# Run benchmarks for a selection of PDB pairs
echo "Running benchmarks for ${#pdb_files[@]} PDB files..."

# Use the first 5 PDB files as reference structures
max_files=$((${#pdb_files[@]} > 5 ? 5 : ${#pdb_files[@]}))
for (( i=0; i<$max_files; i++ )); do
    pdb1="${pdb_files[$i]}"
    
    # Compare with the next 5 PDB files
    for (( j=i+1; j<$max_files; j++ )); do
        pdb2="${pdb_files[$j]}"
        
        # Test different parameters
        for mode in "Sequential" "NonSequential"; do
            for core_dist in 3 4 5; do
                for res_dist in 6 8 10; do
                    run_benchmark "$pdb1" "$pdb2" "$mode" "$core_dist" "$res_dist"
                done
            done
        done
    done
done

echo "Benchmark complete. Results saved to $RESULTS_FILE"