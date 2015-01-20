#/usr/bin/env python
# Allan Bisnar
from Bio.Seq import Seq
from difflib import SequenceMatcher
from math import log10

## Compute similarity matrix for RNA based on the following alignment of
## 2 rRNA sequences

rna1 = "-----CCACCCGGCGAUAGUGAGCGGGCAACACCCGGACUCAUCUCGAACCCGGAAGUAAAG-UCCCCUACGUUGGUAAG-GCA--GUGGGAUCCGCAAGGGCCUGCAGCCUUGCCAAGCUGGGAUGGACAUU"
rna2 = "GAUGGGUGCACGGUCAUAGCGGUGGAGUU-UACCCGGUCUCAUCCCGAACCCGGAAGUCAAGCCCUCCUGCGUC-UGUCCC-AAUACUGUGGUACGAGAGUCCACGGGAACGGCGGUCACUGUG-C-------"

seq1 = Seq(rna1)
seq2 = Seq(rna2)
length_seq1 = len(rna1)
length_seq2 = len(rna2)

print str(seq1.count("A")) + " " + str(seq1.count("C")) + " " + str(seq1.count("U")) + " " + str(seq1.count("G"))
print str(seq2.count("A")) + " " + str(seq2.count("C")) + " " + str(seq2.count("U")) + " " + str(seq2.count("G"))
print length_seq1

pa1 = float(seq1.count("A")) / length_seq1
pa2 = float(seq2.count("A")) / length_seq2

pc1 = float(seq1.count("C")) / length_seq1
pc2 = float(seq2.count("C")) / length_seq2

pu1 = float(seq1.count("U")) / length_seq1
pu2 = float(seq2.count("U")) / length_seq2

pg1 = float(seq1.count("G")) / length_seq1
pg2 = float(seq2.count("G")) / length_seq2

random_model = pa1 * pa2 * pc1 * pc2 * pu1 * pu2 * pg1 * pg2

print random_model

def similar(a,b):
    return SequenceMatcher(None,a,b).ratio()

match_model = similar(rna1,rna2)

print match_model

log_likelihood = log10(match_model / random_model)
print log_likelihood
