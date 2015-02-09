def count_alignments(n, m):
    # Initialize a 3-D Matrix which Stores the 3 separate Matrices Mij, Gx, Gy
    # and their values for a given direction
    S = [[[0 for j in range(m+1)] for i in range(n+1)] for k in range(3)]

    S[1][0][0] = 1

    # Initialize the edges with the given penalties.
    for i in range(1, n+1):
        S[0][i][0] = 1
        
        S[2][i][0] = 0
    for j in range(1, m+1):
        
        S[0][0][j] = 0

        
        S[2][0][j] = 1

    # Compute the Counts for M_ij, Gx, and Gy matrices.
    for j in range(1, m+1):
        for i in range(1, n+1):
            g_x = [S[0][i-1][j], S[1][i-1][j]]
            S[0][i][j] = sum(g_x)

            g_y = [S[2][i][j-1], S[1][i][j-1]]
            S[2][i][j] = sum(g_y)

            m_ij = [S[0][i-1][j-1], S[1][i-1][j-1], S[2][i-1][j-1]]
            S[1][i][j] = sum(m_ij)

    # Get the Total Counts
    matrix_scores = [S[0][i][j], S[1][i][j], S[2][i][j]]
    total_score = sum(matrix_scores)
    print "The alignment count for n =", n," and m =", m," is : " + str(total_score)
    return total_score

count_alignments(1,1) # gap to gap
count_alignments(2,2)
count_alignments(3,3)
count_alignments(4,4)
count_alignments(5,5)
count_alignments(6,6)
count_alignments(7,7)
count_alignments(8,8)
count_alignments(2,3)
