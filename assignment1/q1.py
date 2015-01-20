#/usr/bin/env python
# Allan Bisnar
from Bio.Seq import Seq

## Compute similarity matrix for RNA based on the following alignment of
## 2 rRNA sequences

rna1 = "-----CCACCCGGCGAUAGUGAGCGGGCAACACCCGGACUCAUCUCGAACCCGGAAGUAAAG-UCCCCUACGUUGGUAAG-GCA--GUGGGAUCCGCAAGGGCCUGCAGCCUUGCCAAGCUGGGAUGGACAUU"
rna2 = "GAUGGGUGCACGGUCAUAGCGGUGGAGUU-UACCCGGUCUCAUCCCGAACCCGGAAGUCAAGCCCUCCUGCGUC-UGUCCC-AAUACUGUGGUACGAGAGUCCACGGGAACGGCGGUCACUGUG-C-------"

seq1 = Seq(rna1)
seq2 = Seq(rna2)
length_seq = len(rna1)

print str(seq1.count("A")) + " " + str(seq1.count("C")) + " " + str(seq1.count("U")) + " " + str(seq1.count("G"))
print str(seq2.count("A")) + " " + str(seq2.count("C")) + " " + str(seq2.count("U")) + " " + str(seq2.count("G"))
print length_seq
