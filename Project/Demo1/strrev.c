
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
