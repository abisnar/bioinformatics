
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>    
#include "kseq.h"  
#include "cs445_project.h"
#include <math.h>



	int score(char i, char j){
	if (i == 'A'){
		if (j == 'A') return 10;
		if (j == 'C') return -3;
		if (j == 'G') return -1;
		if (j == 'T') return -4;
	}
	if (i == 'C'){
		if (j == 'A') return -2;
		if (j == 'C') return 4;
		if (j == 'G') return -1;
		if (j == 'T') return -2;
	}
	if (i == 'G'){
		if (j == 'A') return -1;
		if (j == 'C') return -5;
		if (j == 'G') return 7;
		if (j == 'T') return -3;
	}
	if (i == 'T'){
		if (j == 'A') return -4;
		if (j == 'C') return 0;
		if (j == 'G') return -9;
		if (j == 'T') return 8;
	}
}


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

		// rank(s-1, 'c')
		int Rank_s = 0;
		// rank(e, 'c')
		int Rank_e = 0;
		if (s == 1) Rank_s = 0;
		else{
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

	// TODO not distingushing between same score different sequence.  Ignoring traceback after only 1 of i or j hits 0
	void traceback(struct node* node, int j, int r_pos){
		results[r_pos].score = node->max[j];
		int x = 0;
		int y = 0;


		while (node->i*j != 0){
			if (node->max[j] == node->match[j]){

				results[r_pos].n_seq[x] = node->c;
				results[r_pos].p_seq[y] = seq1[j-1];
				
				j--;
				node = node->prev;
			}
			else {
				if (node->max[j] == node->gapX[j]){

					results[r_pos].n_seq[x] = node->c;
					results[r_pos].p_seq[y] = '-';
					node = node->prev;
				}
				else {


					if (node->max[j] == node->gapY[j]){
						results[r_pos].n_seq[x] = '-';
						results[r_pos].p_seq[y] = seq1[j-1];
						j--;
					}
				}
			}
			x++;
			y++;
		}

		results[r_pos].n_seq[x] = '\0';
		results[r_pos].p_seq[y] = '\0';
		//printf("%s\n", results[r_pos].n_seq);
		//printf("%s\n", results[r_pos].p_seq);
	}

	void smith_waterman_recursion(struct node* parent,struct SA_BWT indexes, char c, S_E range){
	//	printf("p-i-%d\n,", parent->i+1);
		// Local Variables
		int length = indexes.pattern_length;
		int j;
		int max = 0;
		int check_pos = 0;	
		S_E rangeA, rangeC, rangeG, rangeT;	

		// Initilize node
		struct node* node = (struct node*) malloc(sizeof(struct node));
	
		node->i = parent->i+1;
		node->prev = parent;
		//printf("%d ", node->i);
		node->c = c;
		node->match = (int*) malloc(sizeof(int)*length);
		node->gapX = (int*) malloc(sizeof(int)*length);
		node->gapY = (int*) malloc(sizeof(int)*length);
		node->max = (int*) malloc(sizeof(int)*length);
		


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
			int check_score, r, pos, temp, j_pos, jj;
			for (j=1; j<length+1; j++){
			//	if (node->i == 2 && j == 2) printf("node->max-%d\n",node->max[j]);
				check_score = 0;
				// results[0] contains top score
				for (r = 0; r < R_SIZE; r++){
					if (check_score == 0){
						
						if (node->max[j] > results[r].score){
				//		printf("node->max-%d results.score-%d\n",node->max[j], results[r].score);
							check_score = 1;
							pos = r;
							j_pos = j;
		
						}
					} 
				}
					if (check_score == 1){ 
					//for (jj=1; jj<length+1; jj++){
					//	printf("Match: %c-%d: %d     ",c,node->i, node->match[jj]); 
					//	printf("XGAP: %c-%d: %d     ",c,node->i, node->gapX[jj]); 
					//	printf("YGAP: %c-%d: %d     ",c,node->i, node->gapY[jj]); 
					//	printf("MAX: %c-%d: %d\n",c,node->i, node->max[jj]); 
					//}
					//	printf("Pos-%d  i-%d  j-%d char-%c\n", pos, node->i, j_pos, node->c);
						for (r = R_SIZE -1; r > pos; r--) results[r] = results[r-1];	
						//  Make function call to enter new alignment_result	
						traceback(node,j_pos,pos);
					//	for (r = 0; r < R_SIZE; r++) printf("Rscore-%d: %d\n", r, results[r].score);
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
		
			

			
			

			// 	How Ideas: 	Remember sequence of chars in each node? 
			//		S_E gives index into suffix array that we can use to get position in original sequence
			// 		traverse parent nodes to recover sequence when a new top score is reached
		

		// Free memory 
		free(node->match);
		free(node->gapX);
		free(node->gapY);
		free(node->max);
		free(node);
	}

	void init_smith_waterman(struct node* root, SA_BWT indexes){
		struct S_E range, rangeA, rangeC, rangeG, rangeT;
		range.s = 1;
		range.e = indexes.length;		
		int i; 
		for (i = 0; i < R_SIZE; i++){
			results[i].score = 0;
			//results[i].n_seq = (char*) malloc(2*sizeof(char)*indexes.pattern_length);
			//results[i].p_seq = (char*) malloc(2*sizeof(char)*indexes.pattern_length);
		}

		rangeA = backwards_search(indexes, 'A', range);
		if (rangeA.s <= rangeA.e) smith_waterman_recursion(root, indexes, 'A', rangeA);

		rangeC = backwards_search(indexes, 'C', range);
		if (rangeC.s <= rangeC.e) smith_waterman_recursion(root, indexes, 'C', rangeC);

		rangeG = backwards_search(indexes, 'G', range);
		if (rangeG.s <= rangeG.e) smith_waterman_recursion(root, indexes, 'G', rangeG);

		rangeT = backwards_search(indexes, 'T', range);
		if (rangeT.s <= rangeT.e) smith_waterman_recursion(root, indexes, 'T', rangeT);

	}





	// Reverse any string
	void strrev(char* s){
		char* t = s;
		while (t && *t) ++t;
		for (--t; s < t; ++s, --t){
			*s = *s ^ *t;
			*t = *s ^ *t;
			*s = *s ^ *t;
		}
	}

	// Compare two strings for quick sort
	static int cmpfunc(const void* a,const void* b){

		const char **ia = (const char **) a;
		const char **ib = (const char **) b;
		return strcmp(*ia, *ib);
	}



	// create BWT and SA
	struct SA_BWT create_bwt(char* seq){
		int length = strlen(seq)+1;
		char* appended_seq = (char*) malloc(sizeof(char)*(length+1)); 
		appended_seq = strcpy(appended_seq, seq);
		appended_seq[length-1] = '$';
		struct permutation permutations[length];
		int i, j, k;
		
		// ROTATIONS OF SEQUENCE

		for (i = 0; i < length; i++){
			permutations[i].string = (char*) malloc(sizeof(char)*length);
			permutations[i].number = i;
			for (j = 0; j < length; j++){
				k = i+j;
				if (k >= length) k = k - length;
				permutations[i].string[j] = appended_seq[k];
			}
	//	printf("%s\n", permutations[i].string);  // To print rotations
		}
		
		// SORT ROTATIONS
		qsort(permutations, length, sizeof(struct permutation), cmpfunc);
		//for (i = 0; i < length; i++) printf("%s\n", permutations[i].string);  //To print sorted strings


		// Calculate C lookup table		(Here we will assume DNA) 
		int C_table[5]; 	// $, A, C, G, T 
		C_table[0] = 0;
		C_table[1] = 1;
		for (i = length-1; i > 0; i--){
			if (permutations[i].string[0] == 'C') C_table[2] = i;
			if (permutations[i].string[0] == 'G') C_table[3] = i;
			if (permutations[i].string[0] == 'T') C_table[4] = i;
		}
		//printf("C_table: $-%d A-%d C-%d G-%d T-%d", C_table[0], C_table[1], C_table[2], C_table[3], C_table[4]);		
		// To print the C_table 
		


		// MAKE BWT SEQUENCE && SUFFIX ARRAY
		char* bwt_seq = (char*) malloc(sizeof(char)*(length+1));
		int* suffix_array = (int*) malloc(sizeof(int)*(length+1));
		for (i = 0; i < length; i++){
			bwt_seq[i] = permutations[i].string[length-1];
			suffix_array[i] = permutations[i].number;

		}
		printf("\nBWT Sequence:\n%s\n", bwt_seq);  			//To print BWT sequence
		//for (i = 0; i < length; i++) printf("%d, ", suffix_array[i]);  	// To print Suffix Array
		//  Free memory 		TODO check for more memory leaks (there will be a couple)
		free ((void*)appended_seq);
		for (i = 0; i < length; i++){
			free(permutations[i].string);
		}

		// RANK Matrix
		int* BWT_Rank[5];
		for (i = 0; i < 5; i++) BWT_Rank[i] = (int*) malloc(sizeof(int)*(length+1));
		int count$ = 0, countA = 0, countC = 0, countG = 0, countT = 0;
		for (i = 0; i < length; i++){
			if (bwt_seq[i] == '$') count$++;
			if (bwt_seq[i] == 'A') countA++;
			if (bwt_seq[i] == 'C') countC++;
			if (bwt_seq[i] == 'G') countG++;
			if (bwt_seq[i] == 'T') countT++;
			BWT_Rank[0][i] = count$;
			BWT_Rank[1][i] = countA;
			BWT_Rank[2][i] = countC;
			BWT_Rank[3][i] = countG;
			BWT_Rank[4][i] = countT;			
		}
		printf("\n\nRank Computed\n");
											//To Print the BWT_Rank Matrix
		//printf("\n\nBWT_Rank matrix:\n");
		//for (i = 0; i < length; i++){
		//	printf("$-%d, A-%d, C-%d, G-%d, T-%d\n", BWT_Rank[0][i], BWT_Rank[1][i], BWT_Rank[2][i], BWT_Rank[3][i], BWT_Rank[4][i]);
		//}
		

		
		
		// Return Struct
		struct SA_BWT seq_SA_BWT;
		seq_SA_BWT.length = length;
		for (i=1; i < 5; i++) seq_SA_BWT.C_table[i] = C_table[i];
		seq_SA_BWT.bwt_seq = bwt_seq;
		seq_SA_BWT.suffix_array = suffix_array;
		for (i=1; i < 5; i++) seq_SA_BWT.BWT_RANK[i] = BWT_Rank[i];

		return seq_SA_BWT;

	}

/* Initial code for reading a fasta file in the main function 
was provided by the kseq package available by searching kseq.h  */

// STEP 1: declare the type of file handler and the read() function  
    KSEQ_INIT(gzFile, gzread)  
      
    int main(int argc, char *argv[])  
    {  
        gzFile fp;  
        kseq_t *seq;  
        int l;  
        if (argc == 1) {  
            fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);  
            return 1;  
        }  
        fp = gzopen(argv[1], "r"); // STEP 2: open the file handler  
        seq = kseq_init(fp); // STEP 3: initialize seq  

	//  START OF INDIVIDUAL CODE

		int i,j;
		struct node* root = (struct node*) malloc(sizeof(struct node));
		root->i = 0;
		root->c = '\0';
		root->match = (int*) malloc(sizeof(int)*1000); // TODO reSIZE
		root->gapX = (int*) malloc(sizeof(int)*1000); // TODO reSIZE
		root->gapY = (int*) malloc(sizeof(int)*1000); // TODO reSIZE
		root->max = (int*) malloc(sizeof(int)*1000); // TODO reSIZE
		root->prev = NULL;

		for (j = 0; j < 1000; j++){		//TODO reSIZE
			root->max[j] = 0;
			root->match[j] = 0;
			if (j > 0) root->match[j] = -D - (j * E);
			root->gapX[j] = M_INF;
			root->gapY[j] = M_INF;
		}
  
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
		strrev(seq->seq.s);
		struct SA_BWT indexes = create_bwt(seq->seq.s);

		indexes.pattern_length = 4;
		init_smith_waterman(root, indexes);
		//printf("\n\n%s\n",results[0].p_seq+1);
		for (i = 0; i < R_SIZE; i++) printf("\nAlignment %d: Score=%d\n%s\n%s\n\n",i, results[i].score, results[i].n_seq, results[i].p_seq);

		

	
	}

	// END OF INDIVIDUAL CODE

        kseq_destroy(seq); // STEP 5: destroy seq  
        gzclose(fp); // STEP 6: close the file handler  
        return 0;  
    }  


