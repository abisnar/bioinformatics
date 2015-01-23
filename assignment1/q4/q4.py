#!/bin/bash/env python
from scripts import BLOSUM62, read_fasta


def global_alignment_affinity_gap(n, m, scoring_matrix, d, e):
    S = [[[0 for j in xrange(len(m)+1)] for i in xrange(len(n)+1)] for k in xrange(3)]

    # Initialize the edges with the given penalties.
    for i in xrange(1, len(n)+1):
        S[0][i][0] = -d - (i-1) * e
        S[1][i][0] = -d - (i-1) * e
        S[2][i][0] = -10000*d
    for j in xrange(1, len(m)+1):
        S[0][0][j] = -10000 * d
        S[1][0][j] = -d - (j-1) * e
        S[2][0][j] = -d - (j-1) * e

    # Fill in the scores for M_ij, Gx, and Gy matrices.
    for j in xrange(1, len(m)+1):
        for i in xrange(1, len(n)+1):
            g_x = [S[0][i-1][j] - e, S[1][i-1][j] - d]
            S[0][i][j] = max(g_x)

            g_y = [S[2][i][j-1] - e, S[1][i][j-1] - d]
            S[2][i][j] = max(g_y)

            m_ij = [S[0][i][j], S[1][i-1][j-1] + scoring_matrix[n[i-1], m[j-1]], S[2][i][j]]
            S[1][i][j] = max(m_ij)

    # Get the maximum score
    matrix_scores = [S[0][i][j], S[1][i][j], S[2][i][j]]
    max_score = max(matrix_scores)

    return str(max_score)


def main():
    # Gap Penalty
    delta = 8
    # Affine Penalty
    eps = 3
    proteins = read_fasta("input/uniprot-organism.fasta")

    results = {}
   # i = recs[0][1]
    with open('input/unknown.txt') as unknown_seq:
        i =unknown_seq.read().strip()
    
    for b in range(0,len(proteins)):
        j = proteins[b][1]

        # Compute Alignment Scores:
        score = global_alignment_affinity_gap(i, j, BLOSUM62(), delta, eps)
        
        results[proteins[b][0]] = int(score)
        print str(b) + ' ' + proteins[b][0] + ' ' + score
    
        with open('output/scores.txt', 'w') as output:
           output.write(str(results))

if __name__ == '__main__':
    main()
