from scripts import *
from numpy import log


def viterbi_local_alignment(seq1, seq2, delta, epsilon, tao, eta, qp_matrix, q_dict):

	# Initialization:
	m = len(seq1)
	n = len(seq2)

	V = [[[-float("inf") for j in xrange(n+1)] for i in xrange(m+1)] for k in xrange(0,13)]
	#backtrack matrix B
	B = [[[-1 for j in xrange(n+1)] for i in xrange(m+1)] for k in xrange(0,13)]

	## V[states][x][y]
	## V[0] = B, V[1] = RX1, V[2] =RY1, V[3] = M, V[4] = X, V[5] = Y, V[6] = RX2 ,
	## V[7] = RY2, V[8] = S1, V[9] = S2, V[10] = S3, V[11] = S4 V[12] = E
	max_val = -float("inf")
	max_i = 0
	max_j = 0
	max_k = 0
	
	V[0][0][0] = 0.0

	# Computation/ Recurrence
	for i in range(0, m+1):
		for j in range(0,n+1):

				#RX1 Cases
				b_rx1 = V[0][i-1][j] + log(1 - eta)
				rx1_rx1 = V[1][i-1][j] + log(1 - eta)
				rx1_vals = [b_rx1, rx1_rx1]

				V[1][i][j] = log(q_dict[seq1[i-1]]) + max(rx1_vals)
				
				if max(rx1_vals) == b_rx1:
					B[1][i][j] = 0
				elif max(rx1_vals) == rx1_rx1:
					B[1][i][j] = 1				

				#RY1 Cases
				ry1_ry1 = V[2][i][j-1] + log(1 - eta)
				ss1_ry1 = V[8][i][j-1] + log(1 - eta)
				ry1_vals = [ry1_ry1,ss1_ry1]

				V[2][i][j] = log(q_dict[seq2[j-1]]) + max(ry1_vals)
				

				if max(ry1_vals) == ry1_ry1:
					B[2][i][j] = 2
				elif max(ry1_vals) == ss1_ry1:
					B[2][i][j] = 8

				#M Cases
				m_m = V[3][i-1][j-1]+ log(1-2*delta - tao)
				x_m = V[4][i-1][j-1]+ log(1 - epsilon - tao)
				y_m = V[5][i-1][j-1]+ log(1 - epsilon - tao)
				ss2_m = V[9][i-1][j-1]+ log(1 - 2*delta - tao)
				
				qpm_score = QPMATRIX().qp[seq1[i-1]][seq2[j-1]]
				m_scores = [m_m, x_m, y_m, ss2_m]

				V[3][i][j] = log(qpm_score) + max( m_scores )
				
				if max(m_scores) == m_m:
					B[3][i][j] = 3
				elif max(m_scores) == x_m:
					B[3][i][j] = 4

				elif max(m_scores) == y_m:
					B[3][i][j] = 5
				elif max(m_scores) == ss2_m:
					B[3][i][j] = 9

				#X Cases
				m_x = V[3][i-1][j]+ log(delta)
				x_x = V[4][i-1][j]+ log(epsilon)
				ss2_x = V[9][i-1][j]+ log(delta)
				x_vals = [ss2_x, m_x, x_x]

				V[4][i][j] = log (q_dict[seq1[i-1]]) + max (x_vals)
				
				if max(x_vals) == m_x:
					B[4][i][j] = 3
				elif max(x_vals) == x_x:
					B[4][i][j] = 4
				elif max(x_vals) == ss2_x:
					B[4][i][j] = 9

				#Y Cases

				m_y = V[3][i][j-1]+ log(delta)
				y_y = V[5][i][j-1]+ log(epsilon)
				ss2_y = V[9][i][j-1]+ log(delta)
				y_vals = [ss2_y,m_y,y_y]

				V[5][i][j] = log (q_dict[seq2[j-1]]) + max (y_vals)
				if max(y_vals) == m_y:
					B[5][i][j] = 3
				elif max(y_vals) == y_y:
					B[5][i][j] = 5
				elif max(y_vals) == ss2_y:
					B[5][i][j] = 9

				#RX2 Cases
				rx2_rx2 = V[6][i-1][j] + log(1 - eta)
				ss3_rx2 = V[10][i-1][j] + log(1 - eta)
				rx2_vals = [ss3_rx2, rx2_rx2]

				V[6][i][j] = log(q_dict[seq1[i-1]]) + max(rx2_vals)
				

				if max(rx2_vals) == rx2_rx2:
					B[6][i][j] = 6
				elif max(rx2_vals) == ss3_rx2:
					B[6][i][j] = 10

				#RY2 Cases
				ry2_ry2 = V[7][i][j-1] + log(1 - eta)
				ss4_ry2 = V[11][i][j-1] + log(1 - eta)
				ry2_vals = [ss4_ry2, ry2_ry2]

				V[7][i][j] = log(q_dict[seq2[j-1]]) + max(ry2_vals)
				
				if max(ry2_vals) == ry2_ry2:
					B[7][i][j] = 7
				else:
					B[7][i][j] = 11

				#V[8] = S1 cases
				b_ss1 = V[0][i][j] + log(eta)
				rx1_ss1 = V[1][i][j] + log(eta)
				s1_vals = [b_ss1,rx1_ss1]

				V[8][i][j] = max(s1_vals)
				
				if max(s1_vals) == b_ss1:
					B[8][i][j] = 0
				elif max(s1_vals) == rx1_ss1:
					B[8][i][j] = 1

				#V[9] = S2,
				ry1_ss2 = V[2][i][j] + log(eta)
				ss1_ss2 = V[8][i][j] + log(eta)
				s2_vals = [ry1_ss2, ss1_ss2]

				V[9][i][j] = max(s2_vals)
				
				if max(s2_vals) == ry1_ss2:
					B[9][i][j] = 2
				elif max(s2_vals) == ss1_ss2:
					B[9][i][j] = 8				

				#V[10] = S3
				m_ss3 = V[3][i][j] + log(tao)
				x_ss3 = V[4][i][j] + log(tao)
				y_ss3 = V[5][i][j] + log(tao)
				ss2_ss3 = V[9][i][j] + log(tao)
				s3_vals = [m_ss3,x_ss3,y_ss3,ss2_ss3]

				V[10][i][j] = max(s3_vals)
				
				if max(s3_vals) == m_ss3:
					B[10][i][j] = 3
				elif max(s3_vals) == x_ss3:
					B[10][i][j] = 4
				elif max(s3_vals) == y_ss3:
					B[10][i][j] = 5
				elif max(s3_vals) == ss2_ss3:
					B[10][i][j] = 9
				#V[11] = S4
				ry2_ss4 = V[6][i][j] + log(eta)
				ss3_ss4 = V[10][i][j] + log(eta)
				s4_vals = [ry2_ss4, ss3_ss4]

				V[11][i][j] = max(s4_vals)
				if max(s4_vals) == ry2_ss4:
					B[11][i][j] = 6
				elif max(s4_vals) == ss3_ss4:
					B[11][i][j] = 10

				#V[12] = E
				ry2_e = V[7][i][j] + log(eta)
				ss4_e = V[11][i][j] + log(eta)
				e_vals = [ry2_e, ss4_e]

				V[12][i][j] = max(e_vals)
				if max(e_vals) == ry2_e:
					B[12][i][j] = 7
				elif max(e_vals) == ss4_e:
					B[12][i][j] = 11

	#Termination
	log_prob = V[12][m][n]

	#print log_prob


	v_aligned = ''
	w_aligned = ''

	k = B[12][m][n]
	i = m
	j = n
	
	while k != 0:

		#RX1
		if k == 1:
			k = B[k][m][n]
			m -= 1
		#RY1
		elif k == 2:
			k = B[k][m][n]
			n -=1
		#M
		elif k == 3:
			k = B[k][m][n]
			i -= 1
			j -= 1
			m -= 1
			n -= 1
			v_aligned += seq1[i]
			w_aligned += seq2[j]
		#X
		elif k == 4:
			k = B[k][m][n]
			i -= 1
			m -= 1
			v_aligned += seq1[i]
			w_aligned += '-'
		#Y
		elif k == 5:
			k = B[k][m][n]
			j -= 1
			n -= 1
			v_aligned += '-'
			w_aligned += seq2[j]
		#RX1
		elif k == 6:
			k = B[k][m][n]
			i -= 1
			m -= 1

		elif k == 7:
			k = B[k][m][n]
			j -= 1
			n -= 1

		elif k == 8:
			k = B[k][m][n]

		elif k == 9:
			k = B[k][m][n]

		elif k == 10:
			k = B[k][m][n]

		elif k == 11:
			k = B[k][m][n]
		
		elif k == 12:
			k = B[k][m][n]
		else:
			break

	v_aligned = v_aligned[::-1]
	w_aligned = w_aligned[::-1]

	seq1_start = i + 1
	seq2_start = j + 1
	t_aligned = w_aligned[0:60]
	b_aligned = v_aligned[0:60]
	len_aligned = max(len(v_aligned), len(w_aligned))

	alignment = {"length": len_aligned, 
					"top_seq" : t_aligned,
					"bot_seq" : b_aligned, 
					"ref_start" : seq1_start, 
					"unk_start" :seq2_start }

	print alignment
	return (log_prob, alignment)

	# ## V[0] = B, V[1] = RX1, V[2] =RY1, V[3] = M, V[4] = X, V[5] = Y, V[6] = RX2 ,
	# ## V[7] = RY2, V[8] = S1, V[9] = S2, V[10] = S3, V[11] = S4 V[12] = E


def main():
	sig = 0.08
	e = 0.35
	t = 0.002
	n = 0.12

	proteins = read_fasta("input/uniprot-organism.fasta")
	ref_protein = proteins[999]
	ref_id = ref_protein[0]
	seq1 = ref_protein[1]

	unknown = read_fasta("input/test.fasta")

	unknown_proteins = unknown[0:998]
	unknown_seq = get_seqs_from_records(unknown_proteins)
	unknown_id = get_seq_id_from_records(unknown_proteins)
	index = list(xrange(998))
	id_index = zip(index, unknown_id)


	result_scores = [viterbi_local_alignment(seq1,seq2, sig, e, t, n, QPMATRIX(),QDICT()) for seq2 in unknown_seq]
	
	all_results = zip(id_index, result_scores)

	# # Sort by max scores
	sorted_by_max_scores = sorted(all_results, key=lambda tup: tup[1][0])[::-1]

	top3 = sorted_by_max_scores[0:3]

	print(top3)

	path_to_output = 'output/top3_results.txt'

	for alignment in top3:
		print "\n"
		print "Index= "+str(alignment[0][0])+" Name= "+alignment[0][1]+"ln Prob= "+str(alignment[1][0])
		print "START POS in "+ref_id+": "+str(alignment[1][1]['ref_start'])+ " START POS in "+alignment[0][1]+": "+ str(alignment[1][1]['ref_start'])+" LENGTH: "+str(alignment[1][1]['length'])
		print alignment[1][1]["top_seq"]
		print alignment[1][1]["bot_seq"]
		print "\n"

	print "printing top 3 to: "+path_to_output

	with open(path_to_output, 'w') as output:
		for alignment in top3:
			output.writelines("\n")
			output.writelines("Index= "+str(alignment[0][0])+" Name= "+alignment[0][1]+"ln Prob= "+str(alignment[1][0]) +"\n")
			output.writelines("START POS in "+ref_id+": "+str(alignment[1][1]['ref_start'])+ " START POS in "+alignment[0][1]+": "+ str(alignment[1][1]['ref_start'])+" LENGTH: "+str(alignment[1][1]['length']) +"\n")
			output.writelines(alignment[1][1]["top_seq"]+"\n")
			output.writelines(alignment[1][1]["bot_seq"]+"\n")
			output.writelines("\n")

if __name__ == '__main__':
    main()
