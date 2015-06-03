	// Backwards search  -- 
	// Takes in the start and end of SA indexes and calculates the next indexes for the given char
	struct S_E backwards_search(struct SA_BWT indexes, char c, struct S_E range){
		
		int s = range.s;
		int e = range.e;
		// C['c']
		int CT = 0;
		if (c == '$') CT = indexes.C_table[0];
		if (c == 'A') CT = indexes.C_table[1];
		if (c == 'C') CT = indexes.C_table[2];
		if (c == 'G') CT = indexes.C_table[3];
		if (c == 'T') CT = indexes.C_table[4];

		int Rank_s = 0;
		int Rank_e = 0;
		if (s > 1){
			if (c == '$') Rank_s = indexes.BWT_RANK[0][s-2];
			if (c == 'A') Rank_s = indexes.BWT_RANK[1][s-2];
			if (c == 'C') Rank_s = indexes.BWT_RANK[2][s-2];
			if (c == 'G') Rank_s = indexes.BWT_RANK[3][s-2];
			if (c == 'T') Rank_s = indexes.BWT_RANK[4][s-2];
		}
		if (c == '$') Rank_e = indexes.BWT_RANK[0][e-1];
		if (c == 'A') Rank_e = indexes.BWT_RANK[1][e-1]; 
		if (c == 'C') Rank_e = indexes.BWT_RANK[2][e-1];
		if (c == 'G') Rank_e = indexes.BWT_RANK[3][e-1];
		if (c == 'T') Rank_e = indexes.BWT_RANK[4][e-1];

		s = CT + Rank_s + 1;
		e = CT + Rank_e;
		range.s = s;
		range.e = e;
		return range;
	}
