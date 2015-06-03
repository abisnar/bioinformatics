
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>    
#include "kseq.h"  
#include "cs445_project_consensus.h"
#include <math.h>
#include "make_bwt.c"
#include "backwards_search.c"
#include "score_consensus.c"
#include "combined_consensus.h"



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


	void traceback(struct node* node, struct SA_BWT indexes, int j, int r_pos, int s){
		results[r_pos].score = node->max[j];
		int x = 0;
		int y = 0;
		int cx;
		int jj;
		for (cx=0; cx < SEQ_SIZE; cx++){ 
			results[r_pos].n_seq[cx] = '\0';
			results[r_pos].p_seq[cx] = '\0';
		}
		while (node->i*j > 0){
			if (node->max[j] == node->gapY[j]){
				results[r_pos].n_seq[x] = '-';
				results[r_pos].p_seq[y] = indexes.pattern[j-1];
				j--;
			}

			else {
				if (node->max[j] == node->gapX[j]){
					results[r_pos].n_seq[x] = node->c;
					results[r_pos].p_seq[y] = '-';
					node = node->prev;
				}
				else {
					if (node->max[j] == node->match[j]){
						results[r_pos].n_seq[x] = node->c;
						results[r_pos].p_seq[y] = indexes.pattern[j-1];
						j--;
						node = node->prev;
					}
				}
			}
			x++;
			y++;
		}
		strrev(results[r_pos].n_seq);
		strrev(results[r_pos].p_seq);
		results[r_pos].s = indexes.suffix_array[s-1];
	}

	void smith_waterman_recursion(struct node* parent,struct SA_BWT indexes, char c, struct S_E range){

		// Local Variables
		int length = indexes.pattern_length;
		int j;
		double max = 0;
		int check_pos = 0;	
		S_E rangeA, rangeC, rangeG, rangeT;	

		// Initilize node
		struct node* node = (struct node*) malloc(sizeof(struct node));
	
		node->i = parent->i+1;
		node->prev = parent;
		node->c = c;
		node->match = (double*) malloc(sizeof(double)*(length+1));
		node->gapX = (double*) malloc(sizeof(double)*(length+1));
		node->gapY = (double*) malloc(sizeof(double)*(length+1));
		node->max = (double*) malloc(sizeof(double)*(length+1));
		


		// Initilize Row (j=0)
		node->match[0] = -D - (node->i * E);
		for (j=1; j < length+1; j++) node->gapX[0] = M_INF;
		for (j=1; j < length+1; j++) node->gapY[0] = M_INF;
		node->max[0] = node->match[0];

		double test;
		// DP calculations
		for (j=1; j<length+1; j++){

			// Switch for Gaps:
			if (indexes.pattern[j-1] == 'Q') indexes.gaps_allowed = 1;
			if (indexes.pattern[j-1] == 'Z') indexes.gaps_allowed = 0;
			if (indexes.gaps_allowed){
				// Match
				test = score_consensus(c, indexes.weighted_seqs, j);
				node->match[j] = parent->max[j-1] + test;
				// X gap
				max = parent->max[j] - D; 
				if (max < parent->gapX[j] - E) max = parent->gapX[j] - E;
				node->gapX[j] = max;
				// Max
				if (max < node->match[j]) max = node->match[j];
				node->max[j] = max;
				if (node->max[j] > 0) check_pos = 1;
			}
			else{
				// Match
				test = score_consensus(c, indexes.weighted_seqs, j);
				node->match[j] = parent->max[j-1] + test;
				max = node->match[j];
				node->max[j] = max;
				if (node->max[j] > 0) check_pos = 1;
			}
		}

	
	

		if (check_pos == 1){
			int check_score, r, pos, temp, j_pos, jj;
			for (j=1; j<length+1; j++){
				check_score = 0;
				// results[0] contains top score
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
						//  Make function call to enter new alignment_result	
						traceback(node,indexes,j_pos,pos, range.s);
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

	void init_smith_waterman(struct node* root, SA_BWT indexes){
		struct S_E range, rangeA, rangeC, rangeG, rangeT;
		range.s = 1;
		range.e = indexes.length;
		indexes.gaps_allowed = 0;
		
		int i; 
		for (i = 0; i < R_SIZE; i++){
			results[i].score = 0;
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
		root->match = (double*) malloc(sizeof(double)*1000); 
		root->gapX = (double*) malloc(sizeof(double)*1000); 
		root->gapY = (double*) malloc(sizeof(double)*1000); 
		root->max = (double*) malloc(sizeof(double)*1000); 
		root->prev = NULL;

		for (j = 0; j < 1000; j++){		
			root->max[j] = 0;
			root->match[j] = 0;
			if (j > 0) root->match[j] = -D - (j * E);
			root->gapX[j] = M_INF;
			root->gapY[j] = M_INF;
		}
  
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
		strrev(seq->seq.s);
		struct SA_BWT indexes = create_bwt(seq->seq.s);
							printf("\n\nBWT Constructed\n\nPress Enter to Search\n");  
							getchar();

		int col = 0;
		int ln = 0;
		indexes.pattern_length = strlen(seq1);
		char* rev_seq1 = (char*) malloc(sizeof(char)*(1+indexes.pattern_length));	
		strcpy(rev_seq1, seq1);
		strrev(rev_seq1);
		indexes.rev_pattern = rev_seq1;
		indexes.pattern = seq1;
		indexes.weighted_seqs[0] = seq1_weightA;
		indexes.weighted_seqs[1] = seq1_weightC;
		indexes.weighted_seqs[2] = seq1_weightG;
		indexes.weighted_seqs[3] = seq1_weightT;
		init_smith_waterman(root, indexes);
		col = (indexes.length - results[0].s)%60 - strlen(results[0].p_seq);
		ln = (indexes.length -results[0].s)/60 + 2;
		if (col < 0){
		col = col + 60;
		ln = ln -1;
		}
		printf("\nAligned pattern: %s\n", indexes.pattern);
		printf("\nTop Alignment: Score=%f Genome_Pos: %d-%d\n", results[0].score, ln, col);
		printf("%s\n%s\n\n",results[0].n_seq, results[0].p_seq);
		results[0].score = 0;

		indexes.pattern_length = strlen(seq2);
		char* rev_seq2 = (char*) malloc(sizeof(char)*(1+indexes.pattern_length));	
		strcpy(rev_seq2, seq2);
		strrev(rev_seq2);
		indexes.rev_pattern = rev_seq2;
		indexes.pattern = seq2;
		indexes.weighted_seqs[0] = seq2_weightA;
		indexes.weighted_seqs[1] = seq2_weightC;
		indexes.weighted_seqs[2] = seq2_weightG;
		indexes.weighted_seqs[3] = seq2_weightT;
		init_smith_waterman(root, indexes);
		col = (indexes.length - results[0].s)%60 - strlen(results[0].p_seq);
		ln = (indexes.length -results[0].s)/60 + 2;
		if (col < 0){
		col = col + 60;
		ln = ln -1;
		}
		printf("\nAligned pattern: %s\n", indexes.pattern);
		printf("\nAlignment: Score=%f Genome_Pos: %d-%d\n", results[0].score, ln, col);
		printf("%s\n%s\n\n",results[0].n_seq, results[0].p_seq);
		results[0].score = 0;

		indexes.pattern_length = strlen(seq3);
		char* rev_seq3 = (char*) malloc(sizeof(char)*(1+indexes.pattern_length));	
		strcpy(rev_seq3, seq3);
		strrev(rev_seq3);
		indexes.rev_pattern = rev_seq3;
		indexes.pattern = seq3;
		indexes.weighted_seqs[0] = seq3_weightA;
		indexes.weighted_seqs[1] = seq3_weightC;
		indexes.weighted_seqs[2] = seq3_weightG;
		indexes.weighted_seqs[3] = seq3_weightT;
		init_smith_waterman(root, indexes);
		col = (indexes.length - results[0].s)%60 - strlen(results[0].p_seq);
		ln = (indexes.length -results[0].s)/60 + 2;
		if (col < 0){
		col = col + 60;
		ln = ln -1;
		}
		printf("\nAligned pattern: %s\n", indexes.pattern);
		printf("\nAlignment: Score=%f Genome_Pos: %d-%d\n", results[0].score, ln, col);
		printf("%s\n%s\n\n",results[0].n_seq, results[0].p_seq);
		results[0].score = 0;

		indexes.pattern_length = strlen(seq4);
		char* rev_seq4 = (char*) malloc(sizeof(char)*(1+indexes.pattern_length));
		strcpy(rev_seq4, seq4);
		strrev(rev_seq4);
		indexes.rev_pattern = rev_seq4;
		indexes.pattern = seq4;
		indexes.weighted_seqs[0] = seq4_weightA;
		indexes.weighted_seqs[1] = seq4_weightC;
		indexes.weighted_seqs[2] = seq4_weightG;
		indexes.weighted_seqs[3] = seq4_weightT;
		init_smith_waterman(root, indexes);
		col = (indexes.length - results[0].s)%60 - strlen(results[0].p_seq);
		ln = (indexes.length -results[0].s)/60 + 2;
		if (col < 0){
		col = col + 60;
		ln = ln -1;
		}
		printf("\nAligned pattern: %s\n", indexes.pattern);
		printf("\nTop Alignment: Score=%f Genome_Pos: %d-%d\n", results[0].score, ln, col);
		printf("%s\n%s\n\n",results[0].n_seq, results[0].p_seq);
		results[0].score = 0;

		indexes.pattern_length = strlen(seq5);
		char* rev_seq5 = (char*) malloc(sizeof(char)*(1+indexes.pattern_length));
		strcpy(rev_seq5, seq5);
		strrev(rev_seq5);
		indexes.rev_pattern = rev_seq5;
		indexes.pattern = seq5;
		indexes.weighted_seqs[0] = seq5_weightA;
		indexes.weighted_seqs[1] = seq5_weightC;
		indexes.weighted_seqs[2] = seq5_weightG;
		indexes.weighted_seqs[3] = seq5_weightT;
		init_smith_waterman(root, indexes);
		col = (indexes.length - results[0].s)%60 - strlen(results[0].p_seq);
		ln = (indexes.length -results[0].s)/60 + 2;
		if (col < 0){
		col = col + 60;
		ln = ln -1;
		}
		printf("\nAligned pattern: %s\n", indexes.pattern);
		printf("\nTop Alignment: Score=%f Genome_Pos: %d-%d\n", results[0].score, ln, col);
		printf("%s\n%s\n\n",results[0].n_seq, results[0].p_seq);
		results[0].score = 0;


/*
		for (i = 0; i < R_SIZE; i++){ 
			col = (indexes.length - results[i].s)%60 - strlen(results[i].p_seq);
			ln = (indexes.length -results[i].s)/60 + 2;
			if (col < 0){
			col = col + 60;
			ln = ln -1;
			}
			printf("Aligned pattern: %s\n\n", indexes.pattern);
			printf("\nAlignment %d: Score=%f Genome_Pos: %d-%d\n",i, results[i].score, ln, col);
			printf("%s\n%s\n\n",results[i].n_seq, results[i].p_seq);
		}

	*/
	}

	// END OF INDIVIDUAL CODE

        kseq_destroy(seq); // STEP 5: destroy seq  
        gzclose(fp); // STEP 6: close the file handler  
        return 0;  
    }  


