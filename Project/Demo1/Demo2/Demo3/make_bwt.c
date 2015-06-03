	

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
		//printf("\nBWT Sequence:\n%s\n", bwt_seq);  			//To print BWT sequence
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
