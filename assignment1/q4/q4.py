#!/bin/bash/env python
from scripts import *
import re

def global_alignment_affinity_gap(n, m, scoring_matrix, d, e):
    # Initialize a 3-D Matrix which Stores the 3 separate Matrices Mij, Gx, Gy
    # and their values for a given direction
    S = [[[0 for j in xrange(len(m)+1)] for i in xrange(len(n)+1)] for k in xrange(3)]

    # Initialize the edges with the given penalties.
    for i in xrange(1, len(n)+1):
        S[0][i][0] = -d - (i-1) * e
        S[1][i][0] = -d - (i-1) * e
        S[2][i][0] = -10000 * d
    for j in xrange(1, len(m)+1):
        S[0][0][j] = -10000 * d
        S[1][0][j] = -d - (j-1) * e
        S[2][0][j] = -d - (j-1) * e

    # Compute the scores for M_ij, Gx, and Gy matrices.
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
    print "The alignment score is : " + str(max_score)
    return max_score

def main():
    # Gap Penalty
    delta = 8
    
    # Affine Penalty
    eps = 3

    # List of all 1017 Uniprot Proteins in Tuples by (id,seq)
    proteins = read_fasta("input/uniprot-organism.fasta")

    sequences = get_seqs_from_records(proteins)
    ids = get_seq_id_from_records(proteins)

    # Access the unknown Protein
    with open('input/unknown.txt') as unknown_seq:
        i = unknown_seq.read().strip()

    # Calculating the Global Alignment score of the unknown to all UniProt Proteins
    result_scores = [global_alignment_affinity_gap(i, j, BLOSUM62(),delta,eps) for j in sequences]

    # Make a new list that contains the ID with Corresponding Alignment Score
    result_list = zip(ids, result_scores)

    # Sort the list by decreasing Alignment Scores
    sorted_by_max_scores = sorted(result_list, key=lambda tup:tup[1])[::-1]

    # Store the Top 3 Candidates
    top3 = sorted_by_max_scores[0:3]

    # Print out the Top 3 Candidates in a txt file in /output/top3_scores.txt
    with open('output/top3_scores.txt', 'w') as output:
        output.write("The Top 3 Global Alignment scores for d = 8 and e = 3 are:" + "\n")
        for elem in top3:
            # REGEX pattern to clean up FASTA ID to return only the Protein Name
            filtered_name = re.search(r'(\d?[A-Z])*\w*[_](\d?[A-Z])*', elem[0])
            output.write('\n' + 'Name: ' + filtered_name.group() +' , Score: ' + str(elem[1]))

# Commented out block was to test q5
   # delta = 4
   # epsilon = 1
   # i = 'XX'
   # j = 'YYY'
   # score = global_alignment_affinity_gap(i, j, BLOSUM62(), delta, epsilon)
   # print(score)

if __name__ == '__main__':
    main()
