#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

"""
usage python3 enc_plot.py input_data.tsv -o statistics.txt -g plot.png

Please note that the input data need to be: Sp annotation  acc ENC (Sp for species, acc for accession number)
"""

def calculate_statistics(input_file, output_file):
    """
    Calculates the mean and standard deviation of ENC for each Sp/annotation combination
    """
    # Read the file
    df = pd.read_csv(input_file, sep='\t')
    
    # Calculate mean and standard deviation for each group
    stats = df.groupby(['Sp', 'annotation'])['ENC'].agg(['mean', 'std']).reset_index()
    stats.columns = ['Sp', 'annotation', 'mean', 'std_dev']
    
    # Replace NaN with 0 for standard deviations (case of a single value)
    stats['std_dev'] = stats['std_dev'].fillna(0)
    
    # Save to file
    stats.to_csv(output_file, sep='\t', index=False)
    print(f"Statistics saved to: {output_file}")
    
    return stats

def create_plot(stats_file, plot_file='ENC_plot.png'):
    """
    Creates a plot with means and standard deviations
    """
    # Read the statistics
    df = pd.read_csv(stats_file, sep='\t')
    
    # Shape and color configuration
    shape_config = {
        'EV-A clade 2': '^',
        'EV-A clade 1': 's',
        'EV-A clade 3': 'v',
        'EV-A clade 4': 'D',
        'EV-B Hs': 's',
        'EV-B NHP': 'D',
        'EV-J': 'o',
        'EV-N': 'o',
        'EV-H': 'o',
        'EV-L': 'o',
        'EV-C NC': '^',
        'EV-C Cult': 's',
        'EV-C resp': 'v',
        'EV-D': 'o',
        'EV-M': 'o'
    }
    
    color_config = {
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
        'EV-M': '#820a07'
    }
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Get unique species and create a numerical position for each
    species = df['Sp'].unique()
    species_positions = {sp: i for i, sp in enumerate(species)}
    
    # To handle point offset
    max_offset = 0.15
    
    # Dictionary to track already plotted annotations (for legend)
    plotted_annotations = {}
    
    # Group by species
    for sp in species:
        df_species = df[df['Sp'] == sp]
        nb_annotations = len(df_species)
        
        # Calculate offsets if multiple annotations
        if nb_annotations > 1:
            offsets = np.linspace(-max_offset, max_offset, nb_annotations)
        else:
            offsets = [0]
        
        # Plot each annotation
        for idx, (_, row) in enumerate(df_species.iterrows()):
            annotation = row['annotation']
            x_pos = species_positions[sp] + offsets[idx]
            y_val = row['mean']
            y_err = row['std_dev']
            
            # Get shape and color
            shape = shape_config.get(annotation, 'o')
            color = color_config.get(annotation, '#000000')
            
            # Plot the point with error bar
            # Add label only if not yet plotted (to avoid duplicates in legend)
            label = annotation if annotation not in plotted_annotations else None
            
            ax.errorbar(x_pos, y_val, yerr=y_err, 
                       fmt=shape, color=color, markersize=10,
                       markeredgecolor='black', markeredgewidth=0.5,
                       capsize=5, capthick=2, elinewidth=2,
                       label=label, alpha=0.8)
            
            if annotation not in plotted_annotations:
                plotted_annotations[annotation] = True
    
    # Axis configuration
    ax.set_xticks(range(len(species)))
    ax.set_xticklabels(species, rotation=45, ha='right')
    ax.set_xlabel('Species', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean ENC', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of ENC values by species and annotation', 
                 fontsize=14, fontweight='bold')
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
             frameon=True, fontsize=9, title='Annotations')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {plot_file}")
    
    # Display
    plt.show()

def main():
    parser = argparse.ArgumentParser(
        description='Analysis of ENC values by species and annotation'
    )
    parser.add_argument('input_file', 
                       help='Input file (TSV format with columns: Sp, annotation, acc, ENC)')
    parser.add_argument('-o', '--output', 
                       default='ENC_statistics.txt',
                       help='Output file for statistics (default: ENC_statistics.txt)')
    parser.add_argument('-g', '--plot', 
                       default='ENC_plot.png',
                       help='Output file for the plot (default: ENC_plot.png)')
    parser.add_argument('--no-plot', 
                       action='store_true',
                       help='Do not create the plot')
    
    args = parser.parse_args()
    
    try:
        # Calculate statistics
        print("Calculating statistics...")
        stats = calculate_statistics(args.input_file, args.output)
        
        # Create the plot
        if not args.no_plot:
            print("\nCreating plot...")
            create_plot(args.output, args.plot)
        
        print("\nProcessing completed successfully!")
        
    except FileNotFoundError:
        print(f"Error: The file '{args.input_file}' does not exist.", 
              file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()