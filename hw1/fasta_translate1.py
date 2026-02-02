"""
Reading in a FASTA file (SHH.fa) into a dictionary object.

Notes: SHH.fa is a single gene FASTA file, but this script considers if we are passing in a multiple gene FASTA file.
"""

import sys

def parse_fasta():
    """ 
    Parse FASTA file into a dictionary object of type {descriptor_str:sequence_str}
    """
    fasta_dict = {}
    current_id = None
    sequence = []

    for line in sys.stdin:
        line = line.strip()

        if line.startswith(">"):
            if current_id: # Stores the previous sequence in dictionary
                fasta_dict[current_id] = "".join(sequence)

            current_id = line[1:]
            sequence = []
        elif line:
            sequence.append(line)
    
    # Save the final sequence
    if current_id:
        fasta_dict[current_id] = "".join(sequence)
    
    return fasta_dict


def main():
    dict = parse_fasta()
    print(dict)


if __name__ == "__main__": main()