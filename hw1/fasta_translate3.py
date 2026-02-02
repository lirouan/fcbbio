"""
Translate a multiple-gene FASTA file (genes.fa) and print in FASTA format.
Uses codon_table.txt
"""

import sys

def load_codon_table(filename):
    """
    Load in codon table and return as a dictionary object of type {codon:amino_acid_symbol}

    Argument:
    filename -- name of codon table in txt format 
    """
    codon_table = {}
    with open(filename, 'r') as f:

        for line in f:
            parts = line.strip().split()
            codon = parts[0] # codon
            aa = parts[1] # amino acid symbol
            codon_table[codon] = aa

    return codon_table


def translate(sequence, codon_table):
    """
    Translate nucleotide sequence into an amino acid sequence
    Assumes no missing/unknown nucleotides.

    Argument:
    sequence -- nucleotide sequence of type string
    codon_table -- dictionary object of type {codon:amino_acid_symbol}
    """
    sequence = sequence.upper()
    protein = []

    # Process in triplets (codons)
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3]
        aa = codon_table[codon] # aa symbol
        protein.append(aa)

    return "".join(protein)


def main():
    codon_table = load_codon_table("codon_table.txt")
    current_id = None
    sequence_data = []

    for line in sys.stdin:
        line = line.strip()

        if line.startswith(">"):
            if current_id:
                print(f">{current_id}")
                print(translate("".join(sequence_data), codon_table))
            
            current_id = line[1:]
            sequence_data = []
        elif line:
            sequence_data.append(line)

    if current_id:
        print(f">{current_id}")
        print(translate("".join(sequence_data), codon_table))


if __name__ == "__main__": main()