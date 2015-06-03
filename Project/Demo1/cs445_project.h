/*
 * File: cs445_project.h
 */

#define M_INF -9999.9
#define D 5
#define E 2
#define R_SIZE 5
#define SEQ_SIZE 500

char seq1[] = "AGGACCATT";

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
	double* weighted_seqs[4];
	int pattern_length;
	char* pattern;
	char* rev_pattern;

} SA_BWT;

typedef struct S_E {
	int s;
	int e;
} S_E;

 struct node {
	int i;
	char c;
	double* match;
	double* gapX;
	double* gapY;
	double* max;
	struct node * prev;
};
typedef struct alignment_results {
	int s;
	double score;
	char n_seq[SEQ_SIZE];
	char p_seq[SEQ_SIZE];
} alignment_results; 

struct alignment_results results[R_SIZE];




