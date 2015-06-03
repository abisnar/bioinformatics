
#include "strrev.c"

	// TODO not distingushing between same score different sequence.  Ignoring traceback after only 1 of i or j hits 0
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
		results[r_pos].n_seq[x] = '\0';
		results[r_pos].s = indexes.suffix_array[s-1];
	}
