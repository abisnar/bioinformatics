	double score_consensus(char i, double** weighted_seqs, int j_pos){
		if (i == 'A') return weighted_seqs[0][j_pos-1];
		if (i == 'C') return weighted_seqs[1][j_pos-1];
		if (i == 'G') return weighted_seqs[2][j_pos-1];
		if (i == 'T') return weighted_seqs[3][j_pos-1];
	}
