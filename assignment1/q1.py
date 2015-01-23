#/usr/bin/env python
# Allan Bisnar
from Bio.Seq import Seq
from difflib import SequenceMatcher
from math import log10

## Compute similarity matrix for RNA based on the following alignment of
## 2 rRNA sequences

rna1 = "-----CCACCCGGCGAUAGUGAGCGGGCAACACCCGGACUCAUCUCGAACCCGGAAGUAAAG-UCCCCUACGUUGGUAAG-GCA--GUGGGAUCCGCAAGGGCCUGCAGCCUUGCCAAGCUGGGAUGGACAUU"
rna2 = "GAUGGGUGCACGGUCAUAGCGGUGGAGUU-UACCCGGUCUCAUCCCGAACCCGGAAGUCAAGCCCUCCUGCGUC-UGUCCC-AAUACUGUGGUACGAGAGUCCACGGGAACGGCGGUCACUGUG-C-------"

#concat both strings
random_rna = rna1 + rna2

length_both = len(random_rna)

seq1 = Seq(rna1)
seq2 = Seq(rna2)
length_seq1 = len(rna1) - seq1.count("-")
length_seq2 = len(rna2) - seq2.count("-")

print str(seq1.count("A")) + " " + str(seq1.count("C")) + " " + str(seq1.count("U")) + " " + str(seq1.count("G"))
print str(seq2.count("A")) + " " + str(seq2.count("C")) + " " + str(seq2.count("U")) + " " + str(seq2.count("G"))
print length_seq1
print length_seq2

pa1 = float(seq1.count("A")) / length_seq1
pa2 = float(seq2.count("A")) / length_seq2

pc1 = float(seq1.count("C")) / length_seq1
pc2 = float(seq2.count("C")) / length_seq2

pu1 = float(seq1.count("U")) / length_seq1
pu2 = float(seq2.count("U")) / length_seq2

pg1 = float(seq1.count("G")) / length_seq1
pg2 = float(seq2.count("G")) / length_seq2

pa = float(random_rna.count("A")) / length_both
pc = float(random_rna.count("C")) / length_both
pu = float(random_rna.count("U")) / length_both
pg = float(random_rna.count("G")) / length_both

#random_model = pa1 * pa2 * pc1 * pc2 * pu1 * pu2 * pg1 * pg2
random_model = pa * pc * pu * pg

print "Pa is:" + str(pa)
print "Pc is:" + str(pc)
print "Pu is:" + str(pu)
print "Pg is:" + str(pg)
print "Random score is: " + str(random_model)

def similar(a,b):
    return SequenceMatcher(None,a,b).ratio()

match_model = similar(rna1,rna2)

print match_model

log_likelihood = log10(match_model / random_model)
print log_likelihood
