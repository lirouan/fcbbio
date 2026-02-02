"""
Translate a single gene FASTA file (SHH.fa) and print in FASTA format
– Use codon_table.txt
"""

import os
import sys

def load_codon_table(filename):
    table = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                table[parts[0]] = parts[1]
    return table

def translate(sequence, codon_table):
    protein = []

    # Process in triplets (codons)
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3].upper()
        amino_acid = codon_table.get(codon, "?")
        protein.append(amino_acid)
    return "".join(protein)

def main():
    codon_table = load_codon_table("codon_table.txt")
    current_id = None
    sequence = []

    for line in sys.stdin:
        line = line.strip()
        if line.startswith(">"):
            current_id = line[1:]
        elif line:
            sequence.append(line)

    full_seq = "".join(sequence)
    protein_seq = translate(full_seq, codon_table)

    print(f">{current_id}")
    print(protein_seq)

if __name__ == "__main__": main()