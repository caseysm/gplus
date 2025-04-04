# GANGSTA+ Benchmark Tools

This directory contains tools for benchmarking the GANGSTA+ protein structure alignment algorithm. The benchmarks use protein structures from the TM-align benchmark dataset to evaluate the performance and results of the modularized GANGSTA+ implementation.

## Available Tools

1. **run_benchmark.sh** - Script to run benchmarks on protein structures
2. **analyze_results.py** - Script to analyze benchmark results and generate visualizations

## Running Benchmarks

### Basic Benchmark

To run a basic benchmark using the test PDB files:

```bash
cd /home/casey/Desktop/repos/gplus
./tests/data/benchmark/run_benchmark.sh
```

This will:
1. Download a subset of the TM-align benchmark dataset (if not already downloaded)
2. Run the modularized GANGSTA+ implementation on multiple PDB pairs
3. Test different parameter combinations (core distance, residue distance)
4. Compare sequential and non-sequential alignment modes
5. Generate CSV files with benchmark results

### Analyzing Results

To analyze the benchmark results and generate visualization plots:

```bash
cd /home/casey/Desktop/repos/gplus
./tests/data/benchmark/analyze_results.py
```

This will:
1. Load the benchmark data from CSV files
2. Generate comparison plots for sequential vs. non-sequential alignment
3. Create heatmaps showing the effect of different parameters
4. Generate stability plots to verify consistent results

The generated plots will be saved in the `/home/casey/Desktop/repos/gplus/tests/data/results` directory.

## Understanding the Results

The benchmark results provide insights into several aspects of the GANGSTA+ algorithm:

1. **Alignment Quality**:
   - Sequential vs. non-sequential alignment differences
   - Number of aligned residues
   - RMSD (Root Mean Square Deviation) values
   - Alignment scores

2. **Parameter Effects**:
   - Impact of core distance parameter
   - Impact of residue distance parameter
   - Optimal parameter values for different protein structures

3. **Performance**:
   - Execution time for different protein pairs
   - Performance scaling with protein size
   - Impact of parameters on performance

4. **Stability**:
   - Consistency of results across multiple runs
   - Deterministic behavior verification

## Adding Custom PDB Files

To benchmark with your own PDB files:

1. Copy your PDB files to the `/home/casey/Desktop/repos/gplus/tests/data/benchmark/pdb` directory
2. Create a `list.txt` file in that directory with the full paths to your PDB files
3. Run the benchmark script as described above

## Requirements

- Bash shell
- Python 3.6+ with the following packages:
  - pandas
  - matplotlib
  - numpy
  - seaborn