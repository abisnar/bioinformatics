#!/usr/bin/env python
import os


class QDICT(object):
    #QDICT class

    def __init__(self):
        #Initialize Qscoring dictionary
        with open(os.path.join(os.path.dirname(__file__), 'data/q_data.txt')) as data:
            libs = [line.strip().split() for line in data.readlines()]
            self.q_val = { lib[0]: float(lib[1]) for lib in libs}

    def __getitem__(self, val):
        #return the score of the Paired Proteins
        return self.q_val[val]
