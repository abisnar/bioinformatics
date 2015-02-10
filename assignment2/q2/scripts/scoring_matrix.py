#!/usr/bin/env python
import os


class BLOSUM62(object):
    #BLOSUM62 class

    def __init__(self):
        #Initialize BLOSUM scoring matrix
        with open(os.path.join(os.path.dirname(__file__), 'data/BLOSUM62.txt')) as data:
            libs = [line.strip().split() for line in data.readlines()]
            self.scoring_matrix = {(lib[0], lib[1]): int(lib[2]) for lib in libs}

    def __getitem__(self, pair):
        #return the score of the Paired Proteins
        return self.scoring_matrix[pair[0], pair[1]]