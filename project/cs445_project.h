/*
 * File: cs445_project.h
 */

#define M_INF -9999999
#define D 5
#define E 2
#define R_SIZE 5

 typedef struct permutation {
	char* string;
	int number;
} permutation;

typedef struct SA_BWT {
	int length;
	int C_table[5];
	char* bwt_seq;
	int* suffix_array;
	int* BWT_RANK[5];
	int pattern_length;
} SA_BWT;

typedef struct S_E {
	int s;
	int e;
} S_E;

 struct node {
	int i;
	char c;
	int* match;
	int* gapX;
	int* gapY;
	int* max;
	struct node * prev;
};
typedef struct alignment_results {
	int score;
	char n_seq[500];
	char p_seq[500];
} alignment_results; 

char seq1[] = "AAGAGCTAGATCACAAA";
//int seq1_weight[] = {1,1,1,2,3,4,3,2,4,3,3};
int topscore =0;
struct alignment_results results[R_SIZE];




