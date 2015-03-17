from scripts import *


def viterbi_local_alignment(seq1, seq2):

	# Initialization:
	m = len(seq1)
	n = len(seq2)

	V = [[[-float("inf") for j in xrange(m + 1)] for i in xrange(n + 1)] for k in xrange(0,13)]
	## V[states][x][y]
	## V[0] = B, V[1] = RX1, V[2] =RY1, V[3] = M, V[4] = X, V[5] = Y, V[6] = RX2 ,
	## V[7] = RY2, V[8] = E, V[9] = S1, V[10] = S2, V[11] = S3, V[12] = S4
	V[0][0][0] = 0

	# # Computation/ Recurrence
	# for j in xrange(0,m + 1):
	# 	for i in xrange(0, n+ 1):
	# 		if i == - 1:
	#
	# 		V[0][i][j] = log(emissions)

	# print qp


	print V

def main():


    #viterbi_local_alignment("A","DE")


if __name__ == '__main__':
    main()
