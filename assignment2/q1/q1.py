#!/bin/bash/env python
# Allan Bisnar
from scripts import *


def local_alignment_linear_gap(v, w, scoring_matrix, d):
    # Returns the local alignment score of v and w with constant gap penalty

    # Initialize the matrix.
    F = [[0 for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
    backtrack = [["-" for j in xrange(len(w)+1)] for i in xrange(len(v)+1)]
    max_val = 0
    max_ij = (0, 0)

    # Fill in the scores for the lower, middle, upper, and backtrack matrices.
    for j in xrange(1, len(w)+1):
        for i in xrange(1, len(v)+1):
            row_scores = F[i-1][j] - d
            col_scores = F[i][j-1] - d
            dia_scores = F[i-1][j-1] + scoring_matrix[v[i-1], w[j-1]]

            scores = [dia_scores, row_scores, col_scores, 0]
            F[i][j] = max(scores)
            if F[i][j] > max_val:
                max_val = F[i][j]
                max_ij = (i, j)

            ##  Scores denoted by gapX, gapY, match, zero
            if F[i][j] == dia_scores:
                backtrack[i][j] = "match"
            elif F[i][j] == row_scores:
                backtrack[i][j] = "gapY"
            elif F[i][j] == col_scores:
                backtrack[i][j] = "gapX"
            else:
                backtrack[i][j] = "zero"

   # Termination
    i = max_ij[0] # location of the max i value
    j = max_ij[1] # location of the max j value


    # Backtrack to the edge of the matrix starting from the max val
    v_aligned = ''
    w_aligned = ''

    while i*j != 0:

        if backtrack[i][j] == "gapY":
            i -= 1
            v_aligned += v[i]
            w_aligned += '-'

        elif backtrack[i][j] == "gapX":
            j -= 1
            v_aligned += '-'
            w_aligned += w[j]

        elif backtrack[i][j] == "match":
            i -= 1
            j -= 1
            v_aligned += str(v[i])
            w_aligned += str(w[j])

        else:
            break

    #Reverse the reconstructed local alignments
    v_aligned = v_aligned[::-1]
    w_aligned = w_aligned[::-1]

    s_start = i + 1
    t_start = j + 1

    alignment = {"length": max(len(v_aligned), len(w_aligned)),
                 "ref seq": v_aligned[0:60],
                 "ref start": s_start,
                 "unknown seq": w_aligned[0:60],
                 "unknown start": t_start}

    print alignment
    return (max_val, alignment)

def main():
    '''Main call. Reads, runs, and saves problem specific data.'''
    # Parse the two input protein strings.
    proteins = read_fasta("input/uniprot-organism.fasta")
    ref_protein = proteins[999]
    ref_id = ref_protein[0]
    seq1 = ref_protein[1]

    unknown_proteins = proteins[0:998]
    unknown_seq = get_seqs_from_records(unknown_proteins)
    unknown_id = get_seq_id_from_records(unknown_proteins)
    index = list(xrange(998))
    id_index = zip(index, unknown_id)

    result_scores = [local_alignment_linear_gap(seq1, seq2, BLOSUM62(), 4) for seq2 in unknown_seq]

    all_results = zip(id_index, result_scores)

    #sort by max_scores
    sorted_by_max_scores = sorted(all_results, key=lambda tup: tup[1][0])[::-1]

    top3 = sorted_by_max_scores[0:3]

    print(top3)

    #print top 3 alignment information

    for alignment in top3:
        print "\n"
        print "Index= "+str(alignment[0][0])+" Name= "+alignment[0][1]+" Score= "+str(alignment[1][0])
        print "START POS in "+ref_id+": "+str(alignment[1][1]['ref start'])+ " START POS in "+alignment[0][1]+": "+ str(alignment[1][1]['unknown start'])+" LENGTH: "+str(alignment[1][1]['length'])
        print alignment[1][1]["ref seq"]
        print alignment[1][1]["unknown seq"]
        print "\n"

    # Print out the Top 3 Candidates in a txt file in /output/top3_scores.txt
    path_to_output = 'output/top3_results.txt'

    print "printing top 3 to: "+path_to_output
    with open(path_to_output, 'w') as output:
        for alignment in top3:
            output.writelines("\n")
            output.writelines("Index= "+str(alignment[0][0])+" Name= "+alignment[0][1]+" Score= "+str(alignment[1][0]) +"\n")
            output.writelines("START POS in "+ref_id+": "+str(alignment[1][1]['ref start'])+ " START POS in "+alignment[0][1]+": "+ str(alignment[1][1]['unknown start'])+" LENGTH: "+str(alignment[1][1]['length'])+"\n")
            output.writelines(alignment[1][1]["ref seq"] + "\n")
            output.writelines(alignment[1][1]["unknown seq"] +"\n")
            output.writelines("\n")

if __name__ == '__main__':
    main()