#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GC3 Calculator - User Guide

Description
This script calculates the GC content at the third codon position (GC3) for DNA sequences 
in FASTA format. GC3 is a useful metric in molecular evolution and codon usage bias studies.

Usage
python3 gc3_gen.py -i input.fasta -o my_results.txt
"""

import argparse
import sys
from pathlib import Path


def read_fasta(fasta_file):
    """
    Reads a FASTA file and returns a dictionary {name: sequence}.
    
    Args:
        fasta_file (str): Path to the FASTA file
        
    Returns:
        dict: Dictionary with sequence names as keys and sequences as values
    """
    sequences = {}
    current_name = None
    current_seq = []
    
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # New sequence
                    if current_name is not None:
                        sequences[current_name] = ''.join(current_seq)
                    current_name = line[1:]  # Remove the '>'
                    current_seq = []
                elif line:
                    current_seq.append(line.upper())
            
            # Add the last sequence
            if current_name is not None:
                sequences[current_name] = ''.join(current_seq)
                
    except FileNotFoundError:
        print(f"Error: The file '{fasta_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return sequences


def calculate_gc3(sequence):
    """
    Calculates the GC percentage at the 3rd codon position.
    
    Args:
        sequence (str): Nucleotide sequence
        
    Returns:
        float: GC3 percentage, or None if the sequence is invalid
    """
    # Check that the sequence is not empty
    if not sequence:
        return None
    
    # Extract nucleotides at position 3 of each codon
    third_positions = sequence[2::3]  # Starts at index 2, then every 3 nucleotides
    
    # Check that there are position 3s
    if not third_positions:
        return None
    
    # Count G and C
    gc_count = third_positions.count('G') + third_positions.count('C')
    total_count = len(third_positions)
    
    # Calculate the percentage
    gc3_percentage = (gc_count / total_count) if total_count > 0 else 0
    #gc3_percentage = (gc_count / total_count) * 100 if total_count > 0 else 0

    return gc3_percentage


def calculate_gc3_for_all(sequences):
    """
    Calculates GC3 for all sequences.
    
    Args:
        sequences (dict): Dictionary {name: sequence}
        
    Returns:
        dict: Dictionary {name: gc3_percentage}
    """
    results = {}
    
    for name, sequence in sequences.items():
        gc3 = calculate_gc3(sequence)
        if gc3 is not None:
            results[name] = gc3
        else:
            print(f"Warning: Unable to calculate GC3 for '{name}'", file=sys.stderr)
    
    return results


def write_results(results, output_file):
    """
    Writes the results to a text file.
    
    Args:
        results (dict): Dictionary {name: gc3_percentage}
        output_file (str): Path to the output file
    """
    try:
        with open(output_file, 'w') as f:
            f.write("Sequence_Name\tGC3_Percentage\n")
            for name, gc3 in results.items():
                f.write(f"{name}\t{gc3:.2f}\n")
        print(f"Results saved to '{output_file}'")
    except Exception as e:
        print(f"Error writing file: {e}", file=sys.stderr)
        sys.exit(1)


def parse_arguments():
    """
    Parses command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Calculates the GC percentage at the 3rd codon position (GC3) for FASTA sequences."
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        '-o', '--output',
        default="gc3_results.txt",
        help="Output file (default: gc3_results.txt)"
    )
    
    return parser.parse_args()


def main():
    """
    Main function of the script.
    """
    # Parse arguments
    args = parse_arguments()
    
    print(f"Reading FASTA file: {args.input}")
    
    # Read the FASTA file
    sequences = read_fasta(args.input)
    
    if not sequences:
        print("No sequences found in the file.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Number of sequences found: {len(sequences)}")
    
    # Calculate GC3 for all sequences
    print("Calculating GC3...")
    results = calculate_gc3_for_all(sequences)
    
    if not results:
        print("No results to save.", file=sys.stderr)
        sys.exit(1)
    
    # Write the results
    write_results(results, args.output)
    
    print(f"Processing completed successfully! ({len(results)} sequences)")


if __name__ == "__main__":
    main()