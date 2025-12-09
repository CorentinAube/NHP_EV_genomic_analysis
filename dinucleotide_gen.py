from Bio import SeqIO
import pandas as pd
from collections import Counter
import itertools
import matplotlib.pyplot as plt
import numpy as np

def calculate_dinucleotide_ratios(fasta_file, output_file):
    """
    Calculates dinucleotide ratios for sequences in a FASTA file.
    
    The ratio is computed as: observed frequency / (freq_nuc1 * freq_nuc2)
    where freq_nuc1 and freq_nuc2 are the frequencies of individual nucleotides.
    
    Args:
        fasta_file (str): Path to input FASTA file
        output_file (str): Path to output TSV file
    """
    nucleotides = "ACGT"
    dinucleotides = [a + b for a, b in itertools.product(nucleotides, repeat=2)]
    results = []
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        seq_length = len(seq)
        
        if seq_length < 2:
            continue
        
        # Calculate single nucleotide frequencies
        nucleotide_counts = Counter(seq)
        total_nucleotides = sum(nucleotide_counts.values())
        nucleotide_freqs = {nuc: nucleotide_counts[nuc] / total_nucleotides for nuc in nucleotides}

        
        # Calculate dinucleotide frequencies
        dinucleotide_counts = Counter(seq[i:i+2] for i in range(seq_length - 1) if seq[i:i+2] in dinucleotides)
        total_dinucleotides = sum(dinucleotide_counts.values())
        dinucleotide_freqs = {dinuc: (dinucleotide_counts[dinuc] / total_dinucleotides) if total_dinucleotides > 0 else 0 for dinuc in dinucleotides}
        
        # Calculate ratio (observed/expected)
        dinucleotide_ratios = {dinuc: (dinucleotide_freqs[dinuc] / (nucleotide_freqs[dinuc[0]] * nucleotide_freqs[dinuc[1]])) if nucleotide_freqs[dinuc[0]] * nucleotide_freqs[dinuc[1]] > 0 else 0 for dinuc in dinucleotides}
        
        results.append({"Sequence": record.id, **dinucleotide_ratios})
    
    # Create DataFrame and save to file
    df = pd.DataFrame(results)
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Dinucleotide ratios saved to: {output_file}")

def main():
    """
    Main function to run dinucleotide ratio analysis.
    """
    fasta_file = "All_seq_oneline.fa"
    output_file = "All_seq_oneline_formean.csv"
    calculate_dinucleotide_ratios(fasta_file, output_file)


if __name__ == '__main__':
    main()