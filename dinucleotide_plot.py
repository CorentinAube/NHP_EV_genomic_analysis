#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to analyze dinucleotide frequencies by annotation
and create visualizations with error bars.
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.lines import Line2D

"""
Please note to add a final column Annotation to the input csv file.
"""



# Marker definition for each annotation
MARKER_MAP = {
    'EV-A clade 2': '^',      # triangle
    'EV-A clade 1': 's',      # square
    'EV-A clade 3': 'v',      # inverted triangle
    'EV-A clade 4': 'D',      # diamond
    'EV-B Hs': 's',           # square
    'EV-B NHP': 'D',          # diamond
    'EV-J': 'o',              # circle
    'EV-N': 'o',              # circle
    'EV-H': 'o',              # circle
    'EV-L': 'o',              # circle
    'EV-C NC': '^',           # triangle
    'EV-C Cult': 's',         # square
    'EV-C resp': 'v',         # inverted triangle
    'EV-D': 'o',              # circle
    'EV-M': 'o',              # circle
}

# Color definition for each annotation
COLOR_MAP = {
    'EV-A clade 2': '#3dff32',
    'EV-A clade 1': '#3dff32',
    'EV-A clade 3': '#3dff32',
    'EV-A clade 4': '#3dff32',
    'EV-B Hs': '#0000ff',
    'EV-B NHP': '#0000ff',
    'EV-J': '#fe6fcf',
    'EV-N': '#ff8001',
    'EV-H': '#68cb6e',
    'EV-L': '#830cff',
    'EV-C NC': '#ff211b',
    'EV-C Cult': '#ff211b',
    'EV-C resp': '#ff211b',
    'EV-D': '#66ccff',
    'EV-M': '#820a07',
}

# List of dinucleotides
DINUCLEOTIDES = ['AA', 'AC', 'AG', 'AU', 'CA', 'CC', 'CG', 'CU', 
                 'GA', 'GC', 'GG', 'GU', 'UA', 'UC', 'UG', 'UU']


def read_data(input_file):
    """
    Reads the data file.
    
    Args:
        input_file (str): Path to the input file
        
    Returns:
        pd.DataFrame: DataFrame with the data
    """
    try:
        df = pd.read_csv(input_file, sep='\t')
        
        # Check that the Annotation column exists
        if 'Annotation' not in df.columns:
            print("Error: Column 'Annotation' does not exist", file=sys.stderr)
            sys.exit(1)
        
        # Check that dinucleotide columns exist
        missing_dinuc = [d for d in DINUCLEOTIDES if d not in df.columns]
        if missing_dinuc:
            print(f"Warning: Missing dinucleotides: {missing_dinuc}")
        
        print(f"Data loaded: {len(df)} sequences")
        print(f"Annotations found: {df['Annotation'].nunique()}")
        
        return df
        
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)


def calculate_statistics(df, dinucleotides):
    """
    Calculates mean and standard deviation for each dinucleotide by annotation.
    
    Args:
        df (pd.DataFrame): DataFrame with the data
        dinucleotides (list): List of dinucleotides to analyze
        
    Returns:
        pd.DataFrame: DataFrame with columns [Annotation, Dinucleotide, Mean, Std, N]
    """
    results = []
    
    for annotation in df['Annotation'].unique():
        # Filter data for this annotation
        subset = df[df['Annotation'] == annotation]
        
        for dinuc in dinucleotides:
            if dinuc in df.columns:
                # Calculate statistics
                values = subset[dinuc].dropna()
                mean_val = values.mean()
                std_val = values.std()
                n = len(values)
                
                results.append({
                    'Annotation': annotation,
                    'Dinucleotide': dinuc,
                    'Mean': mean_val,
                    'Std': std_val,
                    'N': n
                })
    
    stats_df = pd.DataFrame(results)
    print(f"\nStatistics calculated for {len(stats_df)} annotation-dinucleotide combinations")
    
    return stats_df


def save_statistics(stats_df, output_file):
    """
    Saves statistics to a file.
    
    Args:
        stats_df (pd.DataFrame): DataFrame with statistics
        output_file (str): Output file path
    """
    try:
        stats_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
        print(f"Statistics saved: {output_file}")
    except Exception as e:
        print(f"Error saving file: {e}", file=sys.stderr)


def create_full_plot(stats_df, output_file, title="Relative frequencies of dinucleotides"):
    """
    Creates the complete plot with all dinucleotides.
    
    Args:
        stats_df (pd.DataFrame): DataFrame with statistics
        output_file (str): Output file path
        title (str): Plot title
    """
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Get unique list of annotations
    annotations = stats_df['Annotation'].unique()
    dinucleotides = stats_df['Dinucleotide'].unique()
    
    # X positions for dinucleotides
    x_positions = np.arange(len(dinucleotides))
    
    # Offset width to separate points
    n_annotations = len(annotations)
    offset_width = 0.6 / n_annotations if n_annotations > 1 else 0
    
    # For each annotation
    for i, annotation in enumerate(annotations):
        subset = stats_df[stats_df['Annotation'] == annotation]
        
        # Calculate offset for this annotation
        offset = (i - n_annotations / 2) * offset_width + offset_width / 2
        
        # X positions with offset
        x_vals = x_positions + offset
        
        # Extract values for each dinucleotide
        means = []
        stds = []
        for dinuc in dinucleotides:
            dinuc_data = subset[subset['Dinucleotide'] == dinuc]
            if not dinuc_data.empty:
                means.append(dinuc_data['Mean'].values[0])
                stds.append(dinuc_data['Std'].values[0])
            else:
                means.append(np.nan)
                stds.append(np.nan)
        
        # Get marker shape and color
        marker = MARKER_MAP.get(annotation, 'o')
        color = COLOR_MAP.get(annotation, 'gray')
        
        # Plot points with error bars
        ax.errorbar(x_vals, means, yerr=stds, 
                   fmt=marker, 
                   color=color,
                   markersize=8,
                   markeredgecolor='black',
                   markeredgewidth=0.5,
                   capsize=3,
                   capthick=1,
                   elinewidth=1,
                   alpha=0.8,
                   label=annotation)
    
    # Customization
    ax.set_xlabel('Dinucleotides', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean relative frequency', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # X axis
    ax.set_xticks(x_positions)
    ax.set_xticklabels(dinucleotides, rotation=0, ha='center')
    
    # Horizontal line at y=1
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5, zorder=0)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='y')
    
    # Create custom handles for legend (without error bars)
    custom_handles = []
    custom_labels = []
    for annotation in sorted(annotations):
        marker = MARKER_MAP.get(annotation, 'o')
        color = COLOR_MAP.get(annotation, 'gray')
        handle = Line2D([0], [0], marker=marker, color='w', 
                       markerfacecolor=color, markeredgecolor='black',
                       markeredgewidth=0.5, markersize=8, label=annotation)
        custom_handles.append(handle)
        custom_labels.append(annotation)

    # Legend on the right - alphabetically sorted
    ax.legend(custom_handles, custom_labels,
        bbox_to_anchor=(1.05, 1), loc='upper left', 
        fontsize=9, framealpha=0.9, edgecolor='black')
    
    # Improve appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save
    try:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Full plot saved: {output_file}")
    except Exception as e:
        print(f"Error saving file: {e}", file=sys.stderr)
    
    plt.close()


def create_cg_plot(stats_df, output_file, title="Relative frequency of CG dinucleotide"):
    """
    Creates the plot with only the CG dinucleotide.
    
    Args:
        stats_df (pd.DataFrame): DataFrame with statistics
        output_file (str): Output file path
        title (str): Plot title
    """
    # Filter for CG only
    cg_data = stats_df[stats_df['Dinucleotide'] == 'CG'].copy()
    
    if cg_data.empty:
        print("Warning: No data for CG dinucleotide", file=sys.stderr)
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get list of annotations
    annotations = sorted(cg_data['Annotation'].unique())
    
    # X positions for annotations (with spacing)
    x_positions = np.arange(len(annotations))
    
    # For each annotation
    means = []
    stds = []
    colors = []
    markers = []
    
    for annotation in annotations:
        subset = cg_data[cg_data['Annotation'] == annotation]
        means.append(subset['Mean'].values[0])
        stds.append(subset['Std'].values[0])
        colors.append(COLOR_MAP.get(annotation, 'gray'))
        markers.append(MARKER_MAP.get(annotation, 'o'))
    
    # Plot each point individually to have different markers
    for i, annotation in enumerate(annotations):
        ax.errorbar(x_positions[i], means[i], yerr=stds[i],
                   fmt=markers[i],
                   color=colors[i],
                   markersize=12,
                   markeredgecolor='black',
                   markeredgewidth=0.7,
                   capsize=5,
                   capthick=1.5,
                   elinewidth=1.5,
                   alpha=0.8,
                   label=annotation)
    
    # Customization
    ax.set_xlabel('Annotations', fontsize=12, fontweight='bold')
    ax.set_ylabel('CG relative frequency (mean ± std)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # X axis
    ax.set_xticks(x_positions)
    ax.set_xticklabels(annotations, rotation=45, ha='right')
    
    # Horizontal line at y=1
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1, alpha=0.5, zorder=0)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='y')
    
    # Create custom handles for legend (without error bars)
    custom_handles = []
    custom_labels = []
    for annotation in annotations:  # Already alphabetically sorted
        marker = MARKER_MAP.get(annotation, 'o')
        color = COLOR_MAP.get(annotation, 'gray')
        handle = Line2D([0], [0], marker=marker, color='w', 
                       markerfacecolor=color, markeredgecolor='black',
                       markeredgewidth=0.7, markersize=12, label=annotation)
        custom_handles.append(handle)
        custom_labels.append(annotation)

    # Legend on the right - alphabetically sorted
    ax.legend(custom_handles, custom_labels,
             bbox_to_anchor=(1.05, 1), loc='upper left',
             fontsize=9, framealpha=0.9, edgecolor='black')
    # Improve appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save
    try:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"CG plot saved: {output_file}")
    except Exception as e:
        print(f"Error saving file: {e}", file=sys.stderr)
    
    plt.close()


def parse_arguments():
    """
    Parses command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Analysis of dinucleotide frequencies by annotation."
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help="Input TSV file"
    )
    
    parser.add_argument(
        '-o', '--output-prefix',
        default="dinucleotide",
        help="Prefix for output files (default: dinucleotide)"
    )
    
    parser.add_argument(
        '--stats',
        action='store_true',
        help="Save statistics to a TSV file"
    )
    
    parser.add_argument(
        '--title-full',
        default="Relative frequencies of dinucleotides",
        help="Title for the full plot"
    )
    
    parser.add_argument(
        '--title-cg',
        default="Relative frequency of CG dinucleotide",
        help="Title for the CG plot"
    )
    
    return parser.parse_args()


def main():
    """
    Main function of the script.
    """
    # Parse arguments
    args = parse_arguments()
    
    print("="*70)
    print("Dinucleotide analysis by annotation")
    print("="*70)
    
    # Read data
    print(f"\nReading file: {args.input}")
    df = read_data(args.input)
    
    # Display present annotations
    print(f"\nAnnotations present in file:")
    for annotation in sorted(df['Annotation'].unique()):
        count = len(df[df['Annotation'] == annotation])
        print(f"  - {annotation}: {count} sequence(s)")
    
    # Calculate statistics
    print(f"\nCalculating statistics...")
    stats_df = calculate_statistics(df, DINUCLEOTIDES)
    
    # Save statistics if requested
    if args.stats:
        stats_file = f"{args.output_prefix}_statistics.tsv"
        save_statistics(stats_df, stats_file)
    
    # Create full plot
    print(f"\nCreating full plot...")
    full_plot_file = f"{args.output_prefix}_all_dinucleotides.png"
    create_full_plot(stats_df, full_plot_file, args.title_full)
    
    # Create CG plot
    print(f"Creating CG plot...")
    cg_plot_file = f"{args.output_prefix}_CG_only.png"
    create_cg_plot(stats_df, cg_plot_file, args.title_cg)
    
    print("\n" + "="*70)
    print("Processing completed successfully!")
    print("="*70)


if __name__ == "__main__":
    main()