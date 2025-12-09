#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ENC vs GC3 Plot Generator - User Guide

Description
This script visualizes the relationship between ENC (Effective Number of Codons) and GC3 
(GC content at the third codon position). It displays an expected curve and overlays 
observed data points, allowing you to identify sequences with codon usage bias.

You need to modify the output file of the GC3_gen.py code to add 3 columns: Sequence_Name   GC3_Percentage  ENC

python3 enc_gc3_plotter.py -e expected_curve.csv -d observed_data.tsv \
    -o my_plot.png -c "#FF5733" -t "ENC vs GC3 Analysis"
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def read_expected_curve(csv_file):
    """
    Reads the CSV file containing the expected curve.
    
    Args:
        csv_file (str): Path to the CSV file with columns GC3,ENC_expected
        
    Returns:
        pd.DataFrame: DataFrame with columns GC3 and ENC_expected
    """
    try:
        df = pd.read_csv(csv_file)
        
        # Check that required columns are present
        required_cols = ['GC3', 'ENC_expected']
        for col in required_cols:
            if col not in df.columns:
                print(f"Error: Column '{col}' is not present in {csv_file}", 
                      file=sys.stderr)
                sys.exit(1)
        
        # Sort by GC3 for a smooth curve
        df = df.sort_values('GC3')
        
        print(f"Expected curve loaded: {len(df)} points")
        return df
        
    except FileNotFoundError:
        print(f"Error: The file '{csv_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CSV file: {e}", file=sys.stderr)
        sys.exit(1)


def read_observed_data(tsv_file):
    """
    Reads the TSV file containing observed data.
    
    Args:
        tsv_file (str): Path to the TSV file with columns Sequence_Name, GC3_Percentage, ENC
        
    Returns:
        pd.DataFrame: DataFrame with observed data
    """
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        
        # Check that required columns are present
        required_cols = ['Sequence_Name', 'GC3_Percentage', 'ENC']
        for col in required_cols:
            if col not in df.columns:
                print(f"Error: Column '{col}' is not present in {tsv_file}", 
                      file=sys.stderr)
                sys.exit(1)
        
        # Remove missing values
        df_clean = df.dropna(subset=['GC3_Percentage', 'ENC'])
        
        if len(df_clean) < len(df):
            print(f"Warning: {len(df) - len(df_clean)} rows with missing values ignored")
        
        print(f"Observed data loaded: {len(df_clean)} sequences")
        return df_clean
        
    except FileNotFoundError:
        print(f"Error: The file '{tsv_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading TSV file: {e}", file=sys.stderr)
        sys.exit(1)


def create_plot(expected_df, observed_df, point_color, output_file, title):
    """
    Creates the ENC vs GC3 plot.
    
    Args:
        expected_df (pd.DataFrame): DataFrame with expected curve
        observed_df (pd.DataFrame): DataFrame with observed data
        point_color (str): Color of observed points
        output_file (str): Output file path
        title (str): Plot title
    """
    # Create the figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot the expected curve in blue
    ax.plot(expected_df['GC3'], 
            expected_df['ENC_expected'], 
            color='blue', 
            linewidth=2, 
            label=None,#'Expected ENC',
            zorder=1)
    
    # Plot the observed points
    ax.scatter(observed_df['GC3_Percentage'], 
              observed_df['ENC'],
              c=point_color,
              s=100,  # Point size
              alpha=0.6,  # Transparency
              edgecolors='black',  # Black border
              linewidth=0.5,
              label='EV-A clade 1',
              zorder=2)
    
    # Customize the plot
    ax.set_xlabel('GC3 (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('ENC (Effective Number of Codons)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Axis limits
    ax.set_xlim(0, 1)
    ax.set_ylim(20, 61)  # ENC generally ranges between 20 and 61
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # Legend
    ax.legend(loc='best', fontsize=10, framealpha=0.9)
    
    # Improve appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save
    try:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved: {output_file}")
    except Exception as e:
        print(f"Error saving plot: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Display (optional)
    # plt.show()
    
    plt.close()


def parse_arguments():
    """
    Parses command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Visualizes the relationship between ENC and GC3 with expected curve and observed data."
    )
    
    parser.add_argument(
        '-e', '--expected',
        required=True,
        help="CSV file with expected curve (columns: GC3,ENC_expected)"
    )
    
    parser.add_argument(
        '-d', '--data',
        required=True,
        help="TSV file with observed data (columns: Sequence_Name, GC3_Percentage, ENC)"
    )
    
    parser.add_argument(
        '-o', '--output',
        default="enc_vs_gc3_plot.png",
        help="Output file for the plot (default: enc_vs_gc3_plot.png)"
    )
    
    parser.add_argument(
        '-c', '--color',
        default="red",
        help="Color of observed points (default: red). E.g.: red, green, #FF5733, etc."
    )
    
    parser.add_argument(
        '-t', '--title',
        default=None,
        help="Plot title (default: 'ENC vs GC3')"
    )
    
    return parser.parse_args()


def main():
    """
    Main function of the script.
    """
    # Parse arguments
    args = parse_arguments()
    
    print("="*50)
    print("Generating ENC vs GC3 plot")
    print("="*50)
    
    # Read expected curve
    print(f"\nReading expected curve: {args.expected}")
    expected_df = read_expected_curve(args.expected)
    
    # Read observed data
    print(f"Reading observed data: {args.data}")
    observed_df = read_observed_data(args.data)
    
    # Create the plot
    print(f"\nCreating plot...")
    print(f"  - Point color: {args.color}")
    print(f"  - Output file: {args.output}")
    
    create_plot(expected_df, observed_df, args.color, args.output, args.title)
    
    # Statistics
    print("\n" + "="*50)
    print("Statistics:")
    print(f"  - Number of observed points: {len(observed_df)}")
    print(f"  - Mean GC3: {observed_df['GC3_Percentage'].mean():.2f}%")
    print(f"  - Mean ENC: {observed_df['ENC'].mean():.2f}")
    print(f"  - Min ENC: {observed_df['ENC'].min():.2f}")
    print(f"  - Max ENC: {observed_df['ENC'].max():.2f}")
    print("="*50)
    print("\nProcessing completed successfully!")


if __name__ == "__main__":
    main()