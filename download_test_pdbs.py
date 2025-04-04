#!/usr/bin/env python3
"""
Script to download PDB files from the RCSB PDB for GANGSTA+ testing.

This script downloads PDB files directly from the RCSB PDB,
and prepares them for testing with both the original and modular GANGSTA+ implementations.
"""

import os
import sys
import urllib.request
import shutil
import argparse
import subprocess
import time
import re

# Base URL for RCSB PDB
BASE_URL = "https://files.rcsb.org/download/"

# Default PDB IDs to download (representative proteins)
DEFAULT_PROTEINS = [
    "1TIM", "4HHB", "1UBQ", "3SSI", "1LZM",  # Common test proteins
]

def download_pdb(pdb_id, output_dir, overwrite=False):
    """Download a PDB file from RCSB PDB."""
    pdb_id = pdb_id.upper()
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    
    # Skip if file exists and overwrite is False
    if os.path.exists(output_path) and not overwrite:
        print(f"Skipping {pdb_id}.pdb (already exists)")
        return True
    
    url = f"{BASE_URL}{pdb_id}.pdb"
    try:
        print(f"Downloading {pdb_id}.pdb...")
        urllib.request.urlretrieve(url, output_path)
        
        # Verify the file was downloaded properly
        if os.path.getsize(output_path) == 0:
            print(f"Error: Downloaded file {pdb_id}.pdb is empty")
            return False
            
        return True
    except Exception as e:
        print(f"Error downloading {pdb_id}.pdb: {e}")
        return False

def prepare_directories():
    """Prepare directories for testing."""
    # Create directories if they don't exist
    pdb_dir = "pdb_benchmark"
    os.makedirs(pdb_dir, exist_ok=True)
    os.makedirs("comparison_results", exist_ok=True)
    
    # Create original code build directory if it doesn't exist
    original_build_dir = os.path.join("original_code", "build")
    if not os.path.exists(original_build_dir):
        os.makedirs(original_build_dir, exist_ok=True)
    
    # Create modular code build directory if it doesn't exist
    modular_build_dir = "build"
    if not os.path.exists(modular_build_dir):
        os.makedirs(modular_build_dir, exist_ok=True)
    
    return pdb_dir, original_build_dir, modular_build_dir

def build_implementations():
    """Build both implementations."""
    print("Building modular implementation...")
    subprocess.run(["cmake", "-B", "build"], check=True)
    subprocess.run(["cmake", "--build", "build"], check=True)
    
    # Check if original_code exists and build it
    if os.path.exists("original_code") and os.path.isfile(os.path.join("original_code", "CMakeLists.txt")):
        print("Building original implementation...")
        subprocess.run(["cmake", "-B", "original_code/build", "-S", "original_code"], check=True)
        subprocess.run(["cmake", "--build", "original_code/build"], check=True)
        return True
    else:
        print("Original code directory not found, skipping build")
        return False

def copy_pdbs_to_test_dirs(pdb_dir, original_build_dir, modular_build_dir):
    """Copy PDB files to test directories."""
    # Get all PDB files
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(".pdb")]
    
    print(f"Copying {len(pdb_files)} PDB files to test directories...")
    
    # Copy to original build directory
    for pdb_file in pdb_files:
        src = os.path.join(pdb_dir, pdb_file)
        dst_original = os.path.join(original_build_dir, pdb_file)
        dst_modular = os.path.join(modular_build_dir, pdb_file)
        
        # Copy to original directory
        shutil.copy2(src, dst_original)
        
        # Copy to modular directory if different
        if original_build_dir != modular_build_dir:
            shutil.copy2(src, dst_modular)

def run_comparison_tests(pdb_dir, original_build_dir, has_original_impl):
    """Run comparison tests between original and modular implementations."""
    # Get all PDB files
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(".pdb")]
    
    if len(pdb_files) < 2:
        print("Error: Need at least 2 PDB files for comparison tests")
        return
    
    # Create output CSV file
    results_file = os.path.join("comparison_results", "comparison_results.csv")
    with open(results_file, 'w') as f:
        f.write("protein1,protein2,original_aligned,original_rmsd,original_score,modular_aligned,modular_rmsd,modular_score,match\n")
    
    # Total pairs and matches
    total_pairs = 0
    matching_results = 0
    
    print(f"\nRunning comparison tests with {len(pdb_files)} PDB files...")
    
    # Run with all possible pairs for a comprehensive test
    test_pairs = []
    for i, pdb1 in enumerate(pdb_files):
        for pdb2 in pdb_files[i+1:]:
            if pdb1 != pdb2:
                test_pairs.append((pdb1, pdb2))
    
    print(f"Testing {len(test_pairs)} protein pairs...")
    
    for pdb1, pdb2 in test_pairs:
        total_pairs += 1
        pdb1_name = os.path.splitext(pdb1)[0]
        pdb2_name = os.path.splitext(pdb2)[0]
        
        print(f"\nTesting {pdb1_name} vs {pdb2_name}")
        
        # Run modular implementation
        print("Running modular implementation...")
        pdb1_path = os.path.join(pdb_dir, pdb1)
        pdb2_path = os.path.join(pdb_dir, pdb2)
        modular_cmd = [
            "./build/gplus",
            pdb1_path,
            pdb2_path,
            "--residue-distance", "8.0"
        ]
        
        modular_output = os.path.join("comparison_results", f"modular_{pdb1_name}_{pdb2_name}.txt")
        try:
            with open(modular_output, 'w') as outfile:
                subprocess.run(modular_cmd, stdout=outfile, check=True)
            
            # Extract results
            modular_aligned, modular_rmsd, modular_score = extract_results(modular_output)
            print(f"Modular results: Aligned={modular_aligned}, RMSD={modular_rmsd}, Score={modular_score}")
        except Exception as e:
            print(f"Error running modular implementation: {e}")
            modular_aligned, modular_rmsd, modular_score = "ERROR", "ERROR", "ERROR"
        
        # Run original implementation if available
        original_aligned = "N/A"
        original_rmsd = "N/A"
        original_score = "N/A"
        match = "N/A"
        
        if has_original_impl:
            print("Running original implementation...")
            original_output = os.path.join("comparison_results", f"original_{pdb1_name}_{pdb2_name}.txt")
            
            try:
                # We need to copy the PDB files to the original build directory and run from there
                # The original implementation requires files to be in the same directory
                original_cmd = [
                    "./gplus",
                    "-r", "8.0",
                    pdb1,
                    pdb2
                ]
                
                with open(original_output, 'w') as outfile:
                    # Change to original build directory
                    subprocess.run(original_cmd, stdout=outfile, check=True, cwd=original_build_dir)
                
                # Extract results
                original_aligned, original_rmsd, original_score = extract_results(original_output)
                print(f"Original results: Aligned={original_aligned}, RMSD={original_rmsd}, Score={original_score}")
            
                # Check if results match
                if (original_aligned == modular_aligned and 
                    original_rmsd == modular_rmsd and 
                    original_score == modular_score):
                    match = "TRUE"
                    matching_results += 1
                else:
                    match = "FALSE"
            except Exception as e:
                print(f"Error running original implementation: {e}")
                original_aligned, original_rmsd, original_score = "ERROR", "ERROR", "ERROR"
                match = "ERROR"
        
        # Write results to CSV
        with open(results_file, 'a') as f:
            f.write(f"{pdb1_name},{pdb2_name},{original_aligned},{original_rmsd},{original_score},{modular_aligned},{modular_rmsd},{modular_score},{match}\n")
    
    # Print summary
    print("\nComparison Summary:")
    print(f"Total protein pairs tested: {total_pairs}")
    if has_original_impl:
        print(f"Results matching original implementation: {matching_results}/{total_pairs}")
        match_percentage = int(matching_results * 100 / total_pairs) if total_pairs > 0 else 0
        print(f"Match percentage: {match_percentage}%")
    else:
        print("Original implementation not found, no comparison was made")
    
    print(f"Results saved to: {results_file}")
    
    # Generate a summary table for the README
    generate_summary_table(results_file)

def extract_results(output_file):
    """Extract alignment results from output file."""
    aligned = "ERROR"
    rmsd = "ERROR"
    score = "ERROR"
    
    try:
        with open(output_file, 'r') as f:
            content = f.read()
            
            # Extract aligned residues
            aligned_match = re.search(r'Aligned length=\s*(\d+)', content)
            if aligned_match:
                aligned = aligned_match.group(1)
            
            # Extract RMSD
            rmsd_match = re.search(r'RMSD=\s*([\d\.]+)', content)
            if rmsd_match:
                rmsd = rmsd_match.group(1)
            
            # Extract score
            score_match = re.search(r'GP-Score=([\d\.]+)', content)
            if score_match:
                score = score_match.group(1)
    except Exception as e:
        print(f"Error extracting results from {output_file}: {e}")
    
    return aligned, rmsd, score

def generate_summary_table(results_file):
    """Generate a summary table of comparison results."""
    try:
        # Read the results CSV
        with open(results_file, 'r') as f:
            lines = f.readlines()[1:]  # Skip header
        
        # Create markdown table
        summary_file = os.path.join("comparison_results", "summary.md")
        with open(summary_file, 'w') as f:
            f.write("# GANGSTA+ Implementation Comparison\n\n")
            f.write("## Alignment Results\n\n")
            f.write("| Proteins | Original | Modular | Match |\n")
            f.write("|----------|----------|---------|-------|\n")
            
            for line in lines:
                parts = line.strip().split(',')
                if len(parts) >= 9:
                    prot1, prot2, orig_aligned, orig_rmsd, orig_score, mod_aligned, mod_rmsd, mod_score, match = parts
                    
                    # Format the original and modular results
                    if orig_aligned != "N/A" and orig_aligned != "ERROR":
                        orig_result = f"{orig_aligned} residues, RMSD {orig_rmsd}, Score {orig_score}"
                    else:
                        orig_result = orig_aligned
                    
                    if mod_aligned != "N/A" and mod_aligned != "ERROR":
                        mod_result = f"{mod_aligned} residues, RMSD {mod_rmsd}, Score {mod_score}"
                    else:
                        mod_result = mod_aligned
                    
                    # Write the row
                    f.write(f"| {prot1} vs {prot2} | {orig_result} | {mod_result} | {match} |\n")
        
        print(f"Summary table generated: {summary_file}")
    except Exception as e:
        print(f"Error generating summary table: {e}")

def main():
    parser = argparse.ArgumentParser(description='Download and prepare PDB files for GANGSTA+ testing')
    parser.add_argument('--proteins', nargs='+', default=DEFAULT_PROTEINS, 
                        help='PDB IDs to download (default: 5 representative proteins)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing PDB files')
    
    args = parser.parse_args()
    
    # Prepare directories
    pdb_dir, original_build_dir, modular_build_dir = prepare_directories()
    
    # Download PDB files
    print(f"Downloading {len(args.proteins)} PDB files...")
    
    for pdb_id in args.proteins:
        download_pdb(pdb_id, pdb_dir, args.overwrite)
    
    # Build implementations
    has_original_impl = build_implementations()
    
    # Copy PDB files to test directories
    copy_pdbs_to_test_dirs(pdb_dir, original_build_dir, modular_build_dir)
    
    # Run comparison tests
    run_comparison_tests(pdb_dir, original_build_dir, has_original_impl)

if __name__ == "__main__":
    main()