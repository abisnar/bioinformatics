#!/usr/bin/env python
import os


class QPMATRIX(object):
    #QPMATRIX class

    def __init__(self):
        #Initialize QPMATRIX scoring matrix
        with open(os.path.join(os.path.dirname(__file__), 'data/qp_data.txt')) as data:
            self.qp = {}
            letters = []
            first = True
            for line in data.readlines():
              splitted = line.split()
              if first:
                for a in splitted:
                  self.qp[a] = {}
                  letters.append(a)
                first = False
              else:
                a = splitted[0]
                for i in range(1, len(splitted)):
                  b = letters[i-1]
                  self.qp[a][b] = splitted[i]


    def __getitem__(self, p1,p2):
        #return the score of the Paired Proteins
        return self.qp[p1][p2]
