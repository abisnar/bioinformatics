#/bin/bash/env python
from Bio import SeqIO, SubsMat

def global_alignment_affinity_gap(i,j, d, e):
    ## stub for algorithm
    return


records = list(SeqIO.parse("assignment1/q4/uniprot-organism.fasta","fasta"))

for record in records:
    print record.id

def main():
    """ Main, Read Fasta Files, Compute Scores, and Output Results"""
    # Gap Penalty
    d = 8
    # Affine Penalty
    e = 3

    # Parse Input Fasta Files
    for record in records:
        print record.id

    # Compute Alignment Scores:
    score = global_alignment_affinity_gap(i, j, d, e)

    print '\n'.join(score)


if __name__ == '__main__':
    main()
