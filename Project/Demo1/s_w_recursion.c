
	void smith_waterman_recursion(struct node* parent,struct SA_BWT indexes, char c, struct S_E range){

		// Local Variables
		int length = indexes.pattern_length;
		int j, jj;
		double max = 0;
		int check_pos = 0;	
		S_E rangeA, rangeC, rangeG, rangeT;	

		// Initilize node
		struct node* node = (struct node*) malloc(sizeof(struct node));
		node->i = parent->i+1;
		node->prev = parent;
		node->c = c;
		node->match = (double*) malloc(sizeof(double)*length+1);
		node->gapX = (double*) malloc(sizeof(double)*length+1);
		node->gapY = (double*) malloc(sizeof(double)*length+1);
		node->max = (double*) malloc(sizeof(double)*length+1);
		
		// Initilize Row (j=0)
		node->match[0] = -D - (node->i * E);
		node->gapX[0] = M_INF;
		node->gapY[0] = M_INF;
		node->max[0] = node->match[0];

		// DP calculations
		for (j=1; j<length+1; j++){
			// Match
			node->match[j] = parent->max[j-1] + score(c, seq1[j-1]);
			// X gap
			max = parent->max[j] - D; 
			if (max < parent->gapX[j] - E) max = parent->gapX[j] - E;
			node->gapX[j] = max;
			// Y gap
			max = node->max[j-1] - D;
			if (max < node->gapY[j-1] - E) max = node->gapY[j-1] -E;
			node->gapY[j] = max;
			// Max
			if ( max < node->gapX[j]) max = node->gapX[j];
			if ( max < node->match[j]) max = node->match[j];
			node->max[j] = max;
			if (node->max[j] > 0) check_pos = 1;
		}
		if (check_pos == 1){
			int check_score, r, pos, temp, j_pos;
			for (j=1; j<length+1; j++){
				check_score = 0;
				for (r = 0; r < R_SIZE; r++){
					if (check_score == 0){
						if (node->max[j] > results[r].score){
							check_score = 1;
							pos = r;
							j_pos = j;
						}
					} 
				}
					if (check_score == 1){ 
						for (r = R_SIZE -1; r > pos; r--) results[r] = results[r-1];	
						traceback(node,indexes,j_pos,pos,range.s);
					}
				}

			// Recursion step
			rangeA = backwards_search(indexes, 'A', range);
			if (rangeA.s <= rangeA.e) smith_waterman_recursion(node, indexes, 'A', rangeA);
			
			rangeC = backwards_search(indexes, 'C', range);
			if (rangeC.s <= rangeC.e) smith_waterman_recursion(node, indexes, 'C', rangeC);
			rangeG = backwards_search(indexes, 'G', range);
			if (rangeG.s <= rangeG.e) smith_waterman_recursion(node, indexes, 'G', rangeG);

			rangeT = backwards_search(indexes, 'T', range);
			if (rangeT.s <= rangeT.e) smith_waterman_recursion(node, indexes, 'T', rangeT);
		}

		// Free memory 
		free(node->match);
		free(node->gapX);
		free(node->gapY);
		free(node->max);
		free(node);
	}
