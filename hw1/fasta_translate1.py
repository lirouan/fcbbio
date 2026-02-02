"""
Reading in a FASTA file (SHH.fa) into a dictionary object.
"""
import os
import sys

def parse_fasta():
    """ Parse FASTA file into a dictionary object of type {descriptor_str:sequence_str} """
    fasta_dict = {}
    current_id = None
    sequence = []

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id:
                fast_dict[current_id] = "".join(sequence)
            current_id = line[1:]
            sequence = []
        else:
            sequence.append(line)
    
    # Save the final sequence
    if current_id:
        fasta_dict[current_id] = "".join(sequence)
    
    return fasta_dict

def main():
    data = parse_fasta()
    print(data)

if __name__ = "__main__": main()