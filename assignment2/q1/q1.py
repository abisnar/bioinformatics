#!/bin/bash/env python
# Allan Bisnar
from scripts import *
import re

def local_alignment_linear_gap(v, w, scoring_matrix, d):
    # Returns the local alignment score of v and w with constant gap penalty
    # Initialize the matrices.
    S = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
    backtrack = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]

    # Fill in the scores for the lower, middle, upper, and backtrack matrices.
    for i in xrange(1, len(v)+1):
        for j in xrange(1, len(w)+1):
            row_scores = S[i-1][j] - d
            col_scores = S[i][j-1] - d
            dia_scores = S[i-1][j-1] + scoring_matrix[v[i-1], w[j-1]]

            scores = [0, row_scores, col_scores, dia_scores]
            S[i][j] = max(scores)

            ## 0 scores denoted by 0, row scores = 1, col scores, and dia scores = 3
            if S[i][j] == 0:
                backtrack[i][j] = 0
            elif S[i][j] == row_scores:
                backtrack[i][j] = 1
            elif S[i][j] == col_scores:
                backtrack[i][j] = 2
            else:
                backtrack[i][j] = 3

   # Initialize the values of i, j and the aligned sequences.
    i,j = len(v), len(w)
    v_aligned, w_aligned = v, w

    max_score = max([S[i][j]])

    # Lambda function to insert indels.
    insert_indel = lambda word, i: word[:i] + '-' + word[i:]

    # Backtrack to the edge of the matrix starting bottom right.
    while i*j != 0:
        if backtrack[i][j] == 1:
            i -= 1
            w_aligned = insert_indel(w_aligned, j)

        if backtrack[i][j-1] == 2:
            j -= 1
            v_aligned = insert_indel(v_aligned, i)

        if backtrack[i-1][j-1] == 3:
                i -= 1
                j -= 1
        else:
            break


    # Prepend the necessary preceeding indels to get to (0,0).
    for _ in xrange(i):
        w_aligned = insert_indel(w_aligned, 0)
    for _ in xrange(j):
        v_aligned = insert_indel(v_aligned, 0)

    return str(max_score), v_aligned, w_aligned

def main():
    '''Main call. Reads, runs, and saves problem specific data.'''
    # Parse the two input protein strings.
    s = 'PRTEINS'
    t = 'PRTWPSEIN'
    # Get the alignment score.
    score = local_alignment_linear_gap(s, t, BLOSUM62(), 4)

    # Print and save the answer.
    print '\n'.join(score)
    with open('output/q5.txt', 'w') as output_data:
        output_data.write('\n'.join(score))

if __name__ == '__main__':
    main()