
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>    
#include "kseq.h"  
#include "cs445_project.h"
#include <math.h>
#include "score.c"
#include "backwards_search.c"
#include "traceback.c"
#include "s_w_recursion.c"
#include "make_bwt.c"



	void init_smith_waterman(struct node* root, SA_BWT indexes){
		struct S_E range, rangeA, rangeC, rangeG, rangeT;
		range.s = 1;
		range.e = indexes.length;		
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

					/*  START OF  CODE	*/
		// Initialize root
		int i,j;
		struct node* root = (struct node*) malloc(sizeof(struct node));
		root->i = 0;
		root->c = '\0';
		root->match = (double*) malloc(sizeof(double)*SEQ_SIZE*2); 
		root->gapX = (double*) malloc(sizeof(double)*SEQ_SIZE*2); 
		root->gapY = (double*) malloc(sizeof(double)*SEQ_SIZE*2); 
		root->max = (double*) malloc(sizeof(double)*SEQ_SIZE*2); 
		root->prev = NULL;

		for (j = 0; j < SEQ_SIZE*2; j++){		
			root->max[j] = 0;
			root->match[j] = 0;
			if (j > 0) root->match[j] = -D - (j * E);
			root->gapX[j] = M_INF;
			root->gapY[j] = M_INF;
		}
  
        while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 

				printf("\n\nExample sequence (Length %d):\n%s\n\n",strlen(seq->seq.s),seq->seq.s);
				printf("Press ENTER\n");  
				getchar();

		strrev(seq->seq.s);

		struct SA_BWT indexes = create_bwt(seq->seq.s);

		indexes.pattern_length = strlen(seq1);
		char* rev_seq1 = (char*) malloc(sizeof(char)*(1+indexes.pattern_length));
		strcpy(rev_seq1, seq1);
		strrev(rev_seq1);
		indexes.rev_pattern = rev_seq1;
		indexes.pattern = seq1;

		
		init_smith_waterman(root, indexes);


		int col = 0;
		int ln = 0;
			strrev(seq->seq.s);
			printf("Original Sequence:\n%s\n\n", seq->seq.s);
			printf("Sequence being aligned:\n%s\n", indexes.pattern);
			printf("\nTop Alignments:\n");
		for (i = 0; i < R_SIZE; i++){ 
			col = (indexes.length - results[i].s)%60 - strlen(results[i].p_seq)+1;
			ln = (indexes.length -results[i].s)/60 + 1;
			if (col < 0){
			col = col + 60;
			ln = ln -1;
			}
			printf("\nAlignment %d: Score=%f Genome_Pos: Line-%d, Pos-%d\n",i, results[i].score, ln, col);
			printf("%s\n%s\n\n",results[i].n_seq, results[i].p_seq);
		}
	
	}

	// END OF INDIVIDUAL CODE

        kseq_destroy(seq); // STEP 5: destroy seq  
        gzclose(fp); // STEP 6: close the file handler  
        return 0;  
    }  


