"""
Translate a multiple-gene FASTA file containing missing and unknown
nucleotides (messy.fa) and print in FASTA format
- Read-in a more complex codon table: codon_table_hard.txt
"""

import sys

def load_codon_table_hard(filename):
    table = {}

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Split by whitespace (tabs or spaces)
                parts = line.split()
                if len(parts) >= 3:
                    """
                    parts[0]: Amino acid name (e.g., Isoleucine)
                    parts[1]: Symbol (e.g., I)
                    parts[2]: Codon list (e.g., ATT,ATC,ATA)
                    """
                    amino_acid = parts[1]
                    codon = parts[2].split(',')
                    for codon in codons:
                        table[codon.strip().upper()] = amino_acid
    except Exception as e:
        sys.stderr.write(f"Error reading codon_table_hard.txt: {e}\n")
    return table

def translate_sequence(sequence, codon_table):
    protein = []

    # Process sequence in chunks of 3
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3].upper()

        # Handle unknown missing or unknown nucleotides
        if codon in codon_table:
            protein.append(codon_table[codon])
        else:
            sys.stderr.write(f"Missing or unknown codon '{codon}.\n")
            protein.append("X")
    return "".join(protein)

def main():
    codon_table = load_codon_table_hard("codon_table_hard.txt")
    current_id = None
    sequence_data = []

    try:
        for line in sys.stdin:
            line = line.strip()
            
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    print(f">{current_id}")
                    print(translate_sequence("".join(sequence_data), codon_table))
                
                current_id = line[1:]
                sequence_data = []
            else:
                sequence_data.append(line)

        # Sequence the last line in file
        if current_id:
            print(f">{current_id}")
            print(translate_sequence("".join(sequence_data), codon_table))
    
    except Exception as e:
        sys.stderr.write(f"Runtime error: {e}\n")

if __name__ == "__main__": main()
