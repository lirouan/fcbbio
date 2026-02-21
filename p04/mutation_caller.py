"""
mutation_caller.py
A simplified somatic mutation caller inspired by GATK.

Usage:
    python3 mutation_caller.py -n normal.bam -c cancer.bam
"""

# import libraries
import pysam
# print(pysam.__version__)
import math
import sys
import argparse

# constants
MIN_COVERAGE = 20
GENOTYPES = [('A','A'),('C','C'),('G','G'),('T','T'),('A','C'),('A','G'),('A','T'),('C','G'), ('C','T'),('G','T')]

ERROR_RATE = 0.1 # fixed sequencing error rate (e)

LOG_LIKELIHOOD_SKIP = -50 # threshold for ambiguous genotype
LOG_LIKELIHOOD_MUT = -75 # threshold for labeling as somatic mutation candidate


# likelihood helper functions
def base_likelihood(base, genotype = tuple, e: float = ERROR_RATE) -> float:
    
    A1, A2 = genotype

    # probability of base given allele
    # P(b | A) = 1-e if base matches allele, e/3 otherwise
    p_A1 = (1.0 - e) if base == A1 else (e / 3.0)
    p_A2 = (1.0 - e) if base == A2 else (e / 3.0)

    # probability of base given genotype
    # P(b | G) = 1/2 * P(b | A1) + 1/2 * P(b | A2)
    return 0.5 * p_A1 + 0.5 * p_A2

# log P(D | genotype) = sum of log P(base_i | genotype) over all reads
def log_likelihood_genotype(pileup: list, genotype: tuple, e: float = ERROR_RATE) -> float:

    log_ll = 0.0
    for base in pileup:
        p = base_likelihood(base, genotype, e)
        log_ll += math.log(p)

    return log_ll

# return (best_genotype_tuple, best_log_likelihood).
def best_genotype(pileup: list, e: float = ERROR_RATE) -> tuple:

    best_gt = None
    best_ll = float('-inf')

    for gt in GENOTYPES:
        ll = log_likelihood_genotype(pileup, gt, e)
        if ll > best_ll:
            best_ll = ll
            best_gt = gt

    return best_gt, best_ll


# pileup extraction
def get_pileup_bases(bam: pysam.AlignmentFile, chrom: str, pos: int) -> list:
    """
    Return a list of base calls at the given position.
    """
    bases = []
    for pileupcolumn in bam.pileup(chrom, pos, pos+1, truncate=True):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
                base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                bases.append(base)

    return bases


# main run logic
def run(normal_bam_path: str, cancer_bam_path: str):

    # open BAM files
    normal_bam = pysam.AlignmentFile(normal_bam_path, 'rb')
    cancer_bam = pysam.AlignmentFile(cancer_bam_path, 'rb')

    # build a unified set of (chrom, pos) from both BAM pileups
    positions = set()
    for pileupcolumn in normal_bam.pileup():
        positions.add((pileupcolumn.reference_name, pileupcolumn.pos))
    for pileupcolumn in cancer_bam.pileup():
        positions.add((pileupcolumn.reference_name, pileupcolumn.pos))
    positions = sorted(positions)

    normal_bam.reset()
    cancer_bam.reset()

    # get for each position
    for chrom, pos in positions:

        # get pileup bases
        normal_bases = get_pileup_bases(normal_bam, chrom, pos)
        cancer_bases = get_pileup_bases(cancer_bam, chrom, pos)

        if len(normal_bases) < MIN_COVERAGE or len(cancer_bases) < MIN_COVERAGE:
            print(f"Insufficient coverage at position {pos}")
            continue

        # find the best normal genotype
        best_gt, best_ll_normal = best_genotype(normal_bases)

        if best_ll_normal < LOG_LIKELIHOOD_SKIP:
            print(f"Position {pos} has ambiguous genotype")
            continue

        # evaluate cancer pileup under best normal genotype
        ll_cancer = log_likelihood_genotype(cancer_bases, best_gt)

        if ll_cancer < LOG_LIKELIHOOD_MUT:
            print(f"Position {pos} has a candidate somatic mutation (Log-likelihood={ll_cancer:.4f})")

    normal_bam.close()
    cancer_bam.close()


# parse arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--normal', required=True) # path to normal BAM file
    parser.add_argument('-c', '--cancer', required=True) # path to cancer BAM file
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    try:
        run(args.normal, args.cancer)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
