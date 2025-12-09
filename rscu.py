import pandas as pd
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_rscu(sequence):
    """Calculates RSCU values for a DNA sequence."""
    # Check that the sequence has correct length (multiple of 3)
    if len(sequence) % 3 != 0:
        sequence = sequence[:-(len(sequence) % 3)]
    
    # Dictionary to store the number of occurrences of each codon
    codon_count = {}
    
    # Count occurrences of each codon
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3].upper()
        if 'N' not in codon:  # Ignore codons with unknown nucleotides
            if codon in codon_count:
                codon_count[codon] += 1
            else:
                codon_count[codon] = 1
    # Dictionary mapping codons to their amino acids
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    # Calculate RSCU values
    rscu_values = {}
    
    # Group codons by amino acid
    aa_to_codons = {}
    for codon, aa in genetic_code.items():
        if aa not in aa_to_codons:
            aa_to_codons[aa] = []
        aa_to_codons[aa].append(codon)

    # Calculate RSCU for each codon
    for aa, codons in aa_to_codons.items():
        if aa in ['M', 'W', '*']:  # Exclude Met, Trp and stop codons
            continue
        
        # Total number of occurrences of this amino acid
        total_count = sum(codon_count.get(codon, 0) for codon in codons)
        
        # If the amino acid is not present in the sequence, assign default value
        if total_count == 0:
            for codon in codons:
                rscu_values[codon] = 1.0  # Default value
        else:
            # Calculate RSCU for each synonymous codon
            for codon in codons:
                observed = codon_count.get(codon, 0)
                expected = total_count / len(codons)
                rscu_values[codon] = observed / expected if expected > 0 else 0

    return rscu_values

def process_fasta(fasta_file):
    """Loads sequences from a fasta file and calculates RSCU for each sequence."""
    sequences = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        
        # Calculate RSCU for this sequence
        rscu_values = calculate_rscu(sequence)
        sequences[seq_id] = rscu_values
    
    return sequences

def prepare_rscu_matrix(sequences):
    """Prepares a matrix of RSCU values for correspondence analysis."""
    # List of codons to include (exclude Met, Trp and stop codons)
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    codons_to_include = [codon for codon, aa in genetic_code.items() if aa not in ['M', 'W', '*']]
    
    # Create a DataFrame to store RSCU values
    rscu_df = pd.DataFrame(index=sequences.keys(), columns=codons_to_include)
    
    # Fill the DataFrame with RSCU values
    for seq_id, rscu_values in sequences.items():
        for codon in codons_to_include:
            rscu_df.loc[seq_id, codon] = rscu_values.get(codon, 0)
    return rscu_df


def plot_hierarchical_clustering(rscu_matrix, output_file="hierarchical_clustering.png", 
                                 method='average', metric='euclidean', color_threshold=None,
                                 figsize=(12, 8)):
    """
    Performs hierarchical clustering of sequences based on their RSCU values
    and generates a dendrogram with sequence names placed below the figure.
    """
    import scipy.cluster.hierarchy as sch
    from scipy.spatial.distance import pdist
    import numpy as np
    
    # Convert matrix to numpy array of type float64
    numeric_matrix = rscu_matrix.copy()
    
    # Convert each column to numeric
    for col in numeric_matrix.columns:
        numeric_matrix[col] = pd.to_numeric(numeric_matrix[col], errors='coerce').fillna(0)
    
    # Convert DataFrame to numpy array of type float64
    numeric_array = numeric_matrix.astype(np.float64).values
    
    # Replace NaN or inf values with zero
    numeric_array = np.nan_to_num(numeric_array, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Calculate distance matrix between sequences
    distances = pdist(numeric_array, metric=metric)
    
    # Hierarchical clustering
    Z = sch.linkage(distances, method=method)
    
    # Adjust figure size to leave space for names at the bottom
    adjusted_figsize = (figsize[0], figsize[1] + 2)  # Increase height for name space
    
    # Create figure with two subfigures (one for dendrogram, one for labels)
    fig = plt.figure(figsize=adjusted_figsize)
    
    # Create two areas: one for dendrogram (90% top) and one for names (10% bottom)
    # First area goes from [left, bottom, width, height] = [0.1, 0.3, 0.8, 0.6]
    dendrogram_ax = fig.add_axes([0.1, 0.25, 0.8, 0.65])
    
    # Define colors for different clusters
    if color_threshold is None:
        # Automatic threshold estimation
        color_threshold = 0.7 * max(Z[:, 2])
    
    # Create dendrogram without labels (they will be added separately)
    dendrogram_result = sch.dendrogram(
        Z,
        labels=None,  # No labels here, we'll add them manually below
        orientation='top',  # Upward orientation (branches upward)
        color_threshold=0,#color_threshold,
        above_threshold_color='#2ca02c',
        ax=dendrogram_ax,
        no_labels=True  # Completely remove labels
    )
    
    # Customize cluster colors (blue, orange, green)
    #for i, c in zip(range(len(dendrogram_result['color_list'])), ['#1f77b4', '#ff7f0e', '#2ca02c']):
    #    if i < len(dendrogram_result['color_list']):
    #        dendrogram_result['color_list'][i] = c
    
    # Adjust figure parameters
    #dendrogram_ax.set_title('Hierarchical Clustering of sequences', fontsize=16)
    dendrogram_ax.set_xlabel('')
    dendrogram_ax.set_ylabel('Distance', fontsize=14)
    
    # Adjust Y axis
    dendrogram_ax.set_ylim(0, max(Z[:, 2]) * 1.1)
    
    # Remove X axis ticks for cleaner appearance
    dendrogram_ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    
    # Add light grid to facilitate distance reading
    dendrogram_ax.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Find leaf positions (where sequences are in the dendrogram)
    leaf_positions = np.array(dendrogram_result['ivl'], dtype=int)
    
    # Second area for sequence names (bottom)
    # Create a new area below the dendrogram
    names_ax = fig.add_axes([0.1, 0.05, 0.8, 0.1])
    names_ax.set_axis_off()  # No visible axes
    
    # Calculate horizontal positions for names
    # These positions must correspond to leaf positions in the dendrogram
    x_positions = np.arange(len(numeric_matrix.index))
    
    # Get names in dendrogram order
    ordered_names = [numeric_matrix.index[i] for i in leaf_positions]
    
    # Add names at the bottom
    for i, name in enumerate(ordered_names):
        names_ax.text(i, 0.8, name, ha='center', va='center', rotation=90,
                      fontsize=10, fontweight='normal')
    
    # Set name area limits
    names_ax.set_xlim(-0.5, len(numeric_matrix.index) - 0.5)
    names_ax.set_ylim(0, 1)
    
    # Save image
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Hierarchical clustering dendrogram has been saved to {output_file}")
    
    return Z, dendrogram_result

def run_codon_coa(fasta_file):
    """
    Runs complete correspondence analysis on a fasta file,
    with optional coloring by species.
    
    Args:
        fasta_file: Path to FASTA file containing sequences
    """
    print("Loading and processing sequences...")
    sequences = process_fasta(fasta_file)
    
    print("Preparing RSCU matrix...")
    rscu_matrix = prepare_rscu_matrix(sequences)
    
    print("Performing hierarchical clustering of sequences...")
    Z, dendrogram = plot_hierarchical_clustering(
        rscu_matrix, 
        output_file="hierarchical_clustering.png",
        method='average',
        color_threshold=3.5
    )
    
if __name__ == "__main__":
    # Replace "your_file.fasta" with your own FASTA file
    fasta_file = "allSeq_rscu.fasta"
    run_codon_coa(fasta_file)