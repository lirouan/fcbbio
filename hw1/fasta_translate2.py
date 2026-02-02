"""
Translate a single gene FASTA file (SHH.fa) and print in FASTA format
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
            
            if len(parts) == 2:
                codon = parts[0] # amino acid name
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

    current_id = sys.stdin.readline().strip()[1:]
    sequence = []

    # Assumes a single gene FASTA file
    for line in sys.stdin:
        line = line.strip()
        sequence.append(line) # Assumes all lines valid

    print(f">{current_id}")
    print(translate("".join(sequence), codon_table))

if __name__ == "__main__": main()