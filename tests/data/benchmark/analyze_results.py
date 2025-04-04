#!/usr/bin/env python3

"""
GANGSTA+ Benchmark Analysis Script

This script analyzes benchmark results from the GANGSTA+ structural alignment tool
and generates visualizations to help understand the performance characteristics.
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Configuration
RESULTS_DIR = "/home/casey/Desktop/repos/gplus/tests/data/results"
OUTPUT_DIR = RESULTS_DIR

def load_benchmark_data(filename):
    """Load benchmark data from CSV file."""
    file_path = os.path.join(RESULTS_DIR, filename)
    if not os.path.exists(file_path):
        print(f"Error: File {file_path} does not exist.")
        return None
    
    try:
        data = pd.read_csv(file_path)
        print(f"Loaded {len(data)} benchmark results from {filename}")
        return data
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def generate_comparison_plot(data, filename="comparison_plot.png"):
    """Generate a comparison plot for sequential vs. non-sequential alignment."""
    if data is None or len(data) == 0:
        return
        
    # Group data by mode and compute averages
    grouped = data.groupby(['Mode']).agg({
        'AlignedCount': 'mean',
        'RMSD': 'mean',
        'Score': 'mean',
        'Time': 'mean'
    }).reset_index()
    
    # Plot comparison
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Sequential vs. Non-Sequential Alignment Comparison', fontsize=16)
    
    # Plot aligned count
    axs[0, 0].bar(grouped['Mode'], grouped['AlignedCount'])
    axs[0, 0].set_title('Average Aligned Residues')
    axs[0, 0].set_ylabel('Count')
    
    # Plot RMSD
    axs[0, 1].bar(grouped['Mode'], grouped['RMSD'])
    axs[0, 1].set_title('Average RMSD')
    axs[0, 1].set_ylabel('RMSD (Å)')
    
    # Plot score
    axs[1, 0].bar(grouped['Mode'], grouped['Score'])
    axs[1, 0].set_title('Average Score')
    axs[1, 0].set_ylabel('Score')
    
    # Plot time
    axs[1, 1].bar(grouped['Mode'], grouped['Time'])
    axs[1, 1].set_title('Average Runtime')
    axs[1, 1].set_ylabel('Time (s)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, filename))
    print(f"Saved comparison plot to {filename}")

def generate_parameter_heatmap(data, filename="parameter_heatmap.png"):
    """Generate heatmaps showing the effect of parameters on alignment results."""
    if data is None or len(data) == 0:
        return
    
    # Check if we have parameter data
    if 'CoreDistance' not in data.columns or 'ResidueDistance' not in data.columns:
        print("Parameter data not available for heatmap visualization")
        return
    
    # Create separate plots for sequential and non-sequential modes
    for mode in data['Mode'].unique():
        mode_data = data[data['Mode'] == mode]
        
        try:
            # Pivot data for heatmaps
            aligned_pivot = mode_data.pivot_table(
                index='CoreDistance', columns='ResidueDistance', values='AlignedCount', aggfunc='mean')
            rmsd_pivot = mode_data.pivot_table(
                index='CoreDistance', columns='ResidueDistance', values='RMSD', aggfunc='mean')
            score_pivot = mode_data.pivot_table(
                index='CoreDistance', columns='ResidueDistance', values='Score', aggfunc='mean')
            
            # Create figure
            fig, axs = plt.subplots(1, 3, figsize=(18, 5))
            fig.suptitle(f'Parameter Effects on {mode} Alignment', fontsize=16)
            
            # Plot aligned count
            sns.heatmap(aligned_pivot, annot=True, fmt=".1f", ax=axs[0], cmap="YlGnBu")
            axs[0].set_title('Aligned Residues')
            
            # Plot RMSD
            sns.heatmap(rmsd_pivot, annot=True, fmt=".2f", ax=axs[1], cmap="YlOrRd")
            axs[1].set_title('RMSD')
            
            # Plot score
            sns.heatmap(score_pivot, annot=True, fmt=".3f", ax=axs[2], cmap="YlGnBu")
            axs[2].set_title('Score')
            
            plt.tight_layout()
            mode_filename = f"{mode.lower()}_{filename}"
            plt.savefig(os.path.join(OUTPUT_DIR, mode_filename))
            print(f"Saved parameter heatmap for {mode} to {mode_filename}")
        except Exception as e:
            print(f"Error generating parameter heatmap for {mode}: {e}")

def generate_stability_plot(data, filename="stability_plot.png"):
    """Generate plots showing stability of results across multiple runs."""
    if data is None or len(data) == 0:
        return
    
    # Check if we have run number data for stability analysis
    if 'Run' not in data.columns:
        stability_data = None
        # Try to load stability benchmark data directly
        try:
            stability_data = load_benchmark_data("stability_benchmark.csv")
        except:
            print("Stability data not available")
            return
            
        if stability_data is None:
            print("Stability data not available")
            return
        data = stability_data
    
    try:
        # Create plots for sequential and non-sequential modes
        for mode in data['Mode'].unique():
            mode_data = data[data['Mode'] == mode]
            
            # Check if we have PDB1, PDB2 columns for pairing
            if 'PDB1' in data.columns and 'PDB2' in data.columns:
                # Filter by parameters if available
                if 'CoreDistance' in data.columns and 'ResidueDistance' in data.columns:
                    std_params = mode_data[(mode_data['CoreDistance'] == 4) & (mode_data['ResidueDistance'] == 8)]
                else:
                    std_params = mode_data
                
                # Group by PDB pair
                pdb_pairs = std_params.groupby(['PDB1', 'PDB2'])
                
                # Select up to 5 PDB pairs for plotting
                sample_pairs = list(pdb_pairs.groups.keys())[:5]
                
                # Create figure
                fig, axs = plt.subplots(len(sample_pairs), 3, figsize=(15, 3*len(sample_pairs)))
                fig.suptitle(f'Stability Analysis for {mode} Alignment', fontsize=16)
                
                # Single row case
                if len(sample_pairs) == 1:
                    axs = [axs]
                
                for i, (pdb1, pdb2) in enumerate(sample_pairs):
                    pair_data = std_params[(std_params['PDB1'] == pdb1) & (std_params['PDB2'] == pdb2)]
                    
                    # Plot aligned count stability
                    axs[i][0].plot(pair_data['AlignedCount'], 'o-')
                    axs[i][0].set_title(f'{pdb1} vs {pdb2}\nAligned Count')
                    axs[i][0].set_ylim(0, pair_data['AlignedCount'].max() * 1.2)
                    
                    # Plot RMSD stability
                    axs[i][1].plot(pair_data['RMSD'], 'o-')
                    axs[i][1].set_title('RMSD')
                    axs[i][1].set_ylim(0, pair_data['RMSD'].max() * 1.2)
                    
                    # Plot score stability
                    axs[i][2].plot(pair_data['Score'], 'o-')
                    axs[i][2].set_title('Score')
                    axs[i][2].set_ylim(0, pair_data['Score'].max() * 1.2)
                
                plt.tight_layout()
                mode_filename = f"{mode.lower()}_pairs_{filename}"
                plt.savefig(os.path.join(OUTPUT_DIR, mode_filename))
                print(f"Saved stability plot for {mode} pairs to {mode_filename}")
            else:
                # Simple stability plot based on run number
                fig, axs = plt.subplots(1, 3, figsize=(15, 5))
                fig.suptitle(f'Stability Analysis for {mode} Alignment', fontsize=16)
                
                # Plot aligned count stability
                axs[0].plot(mode_data['Run'], mode_data['Aligned'], 'o-')
                axs[0].set_title('Aligned Count')
                axs[0].set_xlabel('Run')
                axs[0].set_ylabel('Count')
                
                # Plot RMSD stability
                axs[1].plot(mode_data['Run'], mode_data['RMSD'], 'o-')
                axs[1].set_title('RMSD')
                axs[1].set_xlabel('Run')
                axs[1].set_ylabel('RMSD')
                
                # Plot score stability
                axs[2].plot(mode_data['Run'], mode_data['Score'], 'o-')
                axs[2].set_title('Score')
                axs[2].set_xlabel('Run')
                axs[2].set_ylabel('Score')
                
                plt.tight_layout()
                mode_filename = f"{mode.lower()}_runs_{filename}"
                plt.savefig(os.path.join(OUTPUT_DIR, mode_filename))
                print(f"Saved stability plot for {mode} runs to {mode_filename}")
    except Exception as e:
        print(f"Error generating stability plot: {e}")

def main():
    """Main function to analyze benchmark results."""
    print("GANGSTA+ Benchmark Analysis")
    print("===========================")
    
    # Check if benchmark data files exist
    benchmark_files = [
        "tm_benchmark_results.csv",
        "comparison_results.csv",
        "performance_benchmark.csv",
        "stability_benchmark.csv"
    ]
    
    # Try to load the benchmark results
    data = None
    for file in benchmark_files:
        data = load_benchmark_data(file)
        if data is not None and len(data) > 0:
            break
    
    if data is None or len(data) == 0:
        print("No benchmark data found. Please run benchmarks first.")
        sys.exit(1)
    
    # Generate plots
    print("\nGenerating visualization plots...")
    generate_comparison_plot(data)
    generate_parameter_heatmap(data)
    generate_stability_plot(data)
    
    print("\nAnalysis complete. Visualizations saved to:", OUTPUT_DIR)

if __name__ == "__main__":
    main()