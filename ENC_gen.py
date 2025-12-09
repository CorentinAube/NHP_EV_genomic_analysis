from Bio import SeqIO
from collections import Counter

def codon_usage(seq):
    """Counts the occurrence of each valid codon in the sequence."""
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3) if len(seq[i:i+3]) == 3]
    return Counter(codons)

def calc_f_hat(codon_counts, total_codons):
    """Calculates F_hat values for synonymous codon classes."""
    synonym_groups = {
        "F": ["TTT", "TTC"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "I": ["ATT", "ATC", "ATA"], "V": ["GTT", "GTC", "GTA", "GTG"], "M": ["ATG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "P": ["CCT", "CCC", "CCA", "CCG"],
        "T": ["ACT", "ACC", "ACA", "ACG"], "A": ["GCT", "GCC", "GCA", "GCG"], "Y": ["TAT", "TAC"],
        "H": ["CAT", "CAC"], "Q": ["CAA", "CAG"], "N": ["AAT", "AAC"], "K": ["AAA", "AAG"],
        "D": ["GAT", "GAC"], "E": ["GAA", "GAG"], "C": ["TGT", "TGC"], "W": ["TGG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], "G": ["GGT", "GGC", "GGA", "GGG"]
    }
    f_hat_values = {}
    for aa, codons in synonym_groups.items():
        if len(codons) > 1:
            total = sum(codon_counts.get(c, 0) for c in codons)
            if total > 1:
                p_values = [(codon_counts.get(c, 0) / total) for c in codons]
                f_hat = ((total * sum(p**2 for p in p_values) - 1) / (total - 1))
                f_hat_values[len(codons)] = f_hat_values.get(len(codons), []) + [f_hat]
    return {k: sum(v) / len(v) if v else 1 for k, v in f_hat_values.items()}

def calculate_enc(seq):
    """Calculates ENC."""
    codon_counts = codon_usage(seq)
    total_codons = sum(codon_counts.values())
    f_hat = calc_f_hat(codon_counts, total_codons)
    enc = 2 + 9 / f_hat.get(2, 1) + 1 / f_hat.get(3, 1) + 5 / f_hat.get(4, 1) + 3 / f_hat.get(6, 1)
    return min(max(enc, 20), 61)  # Constrain within theoretical values

def process_fasta(file_path, output_path):
    """Reads a FASTA file and calculates ENC for each sequence, saves results to a file."""
    with open(output_path, "w") as output_file:
        with open(file_path, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                enc_value = calculate_enc(str(record.seq))
                output_file.write(f"{record.id}\t{enc_value:.2f}\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python3 script.py sequences.fasta output.txt")
    else:
        process_fasta(sys.argv[1], sys.argv[2])