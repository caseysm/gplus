#!/usr/bin/env python3

"""
Implementation Comparison Tool for GANGSTA+

This script compares the debug outputs from the original and modular
implementations of GANGSTA+ to identify algorithmic differences.
"""

import os
import re
import argparse
from collections import defaultdict
import difflib

def parse_debug_file(filepath):
    """Parse a debug output file into structured sections."""
    if not os.path.exists(filepath):
        print(f"Error: File {filepath} not found")
        return None
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Extract all debug lines
    debug_lines = [line.strip() for line in content.split('\n') if line.strip().startswith('[DEBUG]')]
    
    # Organize debug lines into sections
    sections = defaultdict(list)
    
    current_section = "initialization"
    for line in debug_lines:
        # Clean up the line
        line = line.replace('[DEBUG] ', '')
        
        # Check for section headers
        if line.startswith('Starting '):
            current_section = "start"
        elif line.startswith('Step 1:') or line.startswith('Building topology'):
            current_section = "sse_alignment"
        elif line.startswith('Step 2:') or line.startswith('Extending'):
            current_section = "residue_alignment"
        elif line.startswith('Step 3:') or line.startswith('Optimizing'):
            current_section = "optimization"
        elif line.startswith('RMSD calculation') or line.startswith('Calculating RMSD'):
            current_section = "rmsd_calculation"
        elif line.startswith('Alignment result') or line.startswith('Final alignment'):
            current_section = "results"
        
        # Store line in appropriate section
        sections[current_section].append(line)
    
    return sections

def extract_metrics(sections):
    """Extract key metrics from the debug sections."""
    metrics = {
        'sse_count1': None,
        'sse_count2': None,
        'aligned_sses': None,
        'aligned_residues': None,
        'rmsd': None,
        'score': None
    }
    
    # Try to extract metrics from various sections
    for section_name, lines in sections.items():
        for line in lines:
            # SSE counts
            if 'SSE count:' in line:
                count = int(re.search(r'SSE count: (\d+)', line).group(1))
                if metrics['sse_count1'] is None:
                    metrics['sse_count1'] = count
                else:
                    metrics['sse_count2'] = count
            
            # Aligned SSEs
            if 'SSE alignment complete, found' in line:
                metrics['aligned_sses'] = int(re.search(r'found (\d+) SSE pairs', line).group(1))
            
            # Final metrics
            if 'Aligned' in line and 'residues' in line:
                try:
                    metrics['aligned_residues'] = int(re.search(r'Aligned .*?: (\d+)', line).group(1))
                except:
                    pass
            
            if 'RMSD' in line and ':' in line:
                try:
                    metrics['rmsd'] = float(re.search(r'RMSD.*?: ([\d\.]+)', line).group(1))
                except:
                    pass
                
            if 'Score' in line or 'GP-Score' in line:
                try:
                    metrics['score'] = float(re.search(r'Score.*?: ([\d\.]+)', line).group(1))
                except:
                    pass
    
    return metrics

def compare_implementations(original_file, modular_file, output_file):
    """Compare the original and modular implementations."""
    # Parse debug files
    original_sections = parse_debug_file(original_file)
    modular_sections = parse_debug_file(modular_file)
    
    if not original_sections or not modular_sections:
        print("Error: Failed to parse one or both debug files")
        return
    
    # Extract metrics
    original_metrics = extract_metrics(original_sections)
    modular_metrics = extract_metrics(modular_sections)
    
    # Write comparison to output file
    with open(output_file, 'w') as f:
        f.write("# GANGSTA+ Implementation Comparison\n\n")
        
        # Write metrics comparison
        f.write("## Alignment Metrics\n\n")
        f.write("| Metric | Original | Modular | Difference |\n")
        f.write("|--------|----------|---------|------------|\n")
        
        for metric in ['sse_count1', 'sse_count2', 'aligned_sses', 'aligned_residues', 'rmsd', 'score']:
            orig_val = original_metrics.get(metric, "N/A")
            mod_val = modular_metrics.get(metric, "N/A")
            
            if orig_val != "N/A" and mod_val != "N/A":
                if isinstance(orig_val, (int, float)) and isinstance(mod_val, (int, float)):
                    diff = mod_val - orig_val
                    diff_str = f"{diff:+g}" if metric != 'score' else f"{diff:+.4f}"
                else:
                    diff_str = "N/A"
            else:
                diff_str = "N/A"
            
            f.write(f"| {metric.replace('_', ' ').title()} | {orig_val} | {mod_val} | {diff_str} |\n")
        
        # Write section comparisons
        for section in ['initialization', 'sse_alignment', 'residue_alignment', 'rmsd_calculation', 'results']:
            f.write(f"\n## {section.replace('_', ' ').title()} Section\n\n")
            
            if section in original_sections and section in modular_sections:
                orig_text = "\n".join(original_sections[section])
                mod_text = "\n".join(modular_sections[section])
                
                # Generate diff
                diff = difflib.ndiff(orig_text.splitlines(), mod_text.splitlines())
                f.write("```diff\n" + "\n".join(diff) + "\n```\n")
            else:
                if section not in original_sections:
                    f.write("Section not found in original implementation\n")
                if section not in modular_sections:
                    f.write("Section not found in modular implementation\n")
    
    print(f"Comparison written to {output_file}")
    
    # Print summary to console
    print("\nAlignment Metrics Comparison:")
    print(f"Original: {original_metrics['aligned_residues']} residues, RMSD: {original_metrics['rmsd']}, Score: {original_metrics['score']}")
    print(f"Modular:  {modular_metrics['aligned_residues']} residues, RMSD: {modular_metrics['rmsd']}, Score: {modular_metrics['score']}")

def main():
    parser = argparse.ArgumentParser(description="Compare GANGSTA+ implementations")
    parser.add_argument("original", help="Debug output from original implementation")
    parser.add_argument("modular", help="Debug output from modular implementation")
    parser.add_argument("--output", "-o", default="comparison_report.md", help="Output file for comparison report")
    
    args = parser.parse_args()
    compare_implementations(args.original, args.modular, args.output)

if __name__ == "__main__":
    main()