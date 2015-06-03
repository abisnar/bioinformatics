	double score(char i, char j){
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

