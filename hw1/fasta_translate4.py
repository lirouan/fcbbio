"""
Translate a multiple-gene FASTA file containing missing and unknown nucleotides (messy.fa) and print in FASTA format
Uses: codon_table_hard.txt
"""

import sys

def load_codon_table(filename):
    """
    Load in codon table (hard) and return as a dictionary object of type {codon:amino_acid_symbol}

    Argument:
    filename -- name of codon table in txt format 
    """
    codon_table = {}

    print 

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Split by whitespace (tabs or spaces)
            parts = line.split()

            aa = parts[0] # amino acid name
            symbol = parts[1] # amino acid symbol
            codons = parts[2].split(",") # possible codons
            
            for codon in codons:
                codon_table[codon.strip().upper()] = symbol

    return codon_table
    

def translate(sequence, codon_table):
    """
    Translate nucleotide sequence into an amino acid sequence
    Handle missing and unknown nucleotides by using an X for the translated amino acid

    Argument:
    sequence -- nucleotide sequence of type string
    codon_table -- dictionary object of type {codon:amino_acid_symbol}
    """
    sequence = sequence.upper()
    protein = []

    # Process sequence in chunks of 3 nucleotides
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3]

        # Translate codon to amino acid
        if codon in codon_table:
            protein.append(codon_table[codon])
        else:
            # sys.stderr.write(f"Missing/unknown nucleotide {codon}.\n")
            protein.append("X")

    return "".join(protein)


def main():
    codon_table = load_codon_table("codon_table_hard.txt")
    current_id = None
    sequence_data = []

    # Read in lines of the FASTA nucleotide sequence
    for line in sys.stdin:
        line = line.strip()
        
        if line.startswith(">"):
            if current_id:
                print(f">{current_id}")
                print(translate("".join(sequence_data), codon_table))
            
            current_id = line[1:]
            sequence_data = []
        elif line:
            sequence_data.append(line) # Append nucleotides for the current sequence ID

    # Print the last sequence
    if current_id:
        print(f">{current_id}")
        print(translate("".join(sequence_data), codon_table))


if __name__ == "__main__": main()
