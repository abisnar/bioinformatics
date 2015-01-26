#/usr/bin/env python
# Allan Bisnar
from math import log

## Compute similarity matrix for RNA based on the following alignment of
## 2 rRNA rnauences

rna1 = "-----CCACCCGGCGAUAGUGAGCGGGCAACACCCGGACUCAUCUCGAACCCGGAAGUAAAG-UCCCCUACGUUGGUAAG-GCA--GUGGGAUCCGCAAGGGCCUGCAGCCUUGCCAAGCUGGGAUGGACAUU"
rna2 = "GAUGGGUGCACGGUCAUAGCGGUGGAGUU-UACCCGGUCUCAUCCCGAACCCGGAAGUCAAGCCCUCCUGCGUC-UGUCCC-AAUACUGUGGUACGAGAGUCCACGGGAACGGCGGUCACUGUG-C-------"

match = zip(rna1, rna2)
no_gap_match = []

for elem in match:
	if not '-' in elem:
		no_gap_match.append(elem)

length_match = len(no_gap_match)

AA = 0
AU = 0
AG = 0
AC = 0
CU = 0
CG = 0
CC = 0
GG = 0
GU = 0
UU = 0

for elem in no_gap_match:
	if elem[0] == 'A' and elem[1] == 'A':
		AA +=1
	elif (elem[0] == 'A' and elem[1] == 'U') or (elem[0] == 'U' and elem[1] == 'A'):
		AU +=1
	elif (elem[0] == 'A' and elem[1] == 'G') or (elem[0] == 'G' and elem[1] == 'A'):
		AG +=1
	elif (elem[0] == 'A' and elem[1] == 'C') or (elem[0] == 'C' and elem[1] == 'A'):
		AC +=1
	elif (elem[0] == 'C' and elem[1] == 'U') or (elem[0] == 'U' and elem[1] == 'C'):
		CU +=1
	elif (elem[0] == 'C' and elem[1] == 'G') or (elem[0] == 'G' and elem[1] == 'C'):
		CG +=1
	elif (elem[0] == 'C' and elem[1] == 'C'):
		CC +=1
	elif (elem[0] == 'G' and elem[1] == 'G'):
		GG +=1
	elif (elem[0] == 'G' and elem[1] == 'U') or (elem[0] == 'U' and elem[1] == 'G'):
		GU +=1
	else:
		UU +=1

pAA = float(AA)/length_match
pAU = float(AU)/length_match
pAC = float(AC)/length_match
pAG = float(AG)/length_match
pCU = float(CU)/length_match
pCG = float(CG)/length_match
pCC = float(CC)/length_match
pGG = float(GG)/length_match
pGU = float(GU)/length_match
pUU = float(UU)/length_match

print "========================"
print "pAA is : " +str(pAA)
print "pAU is : " +str(pAU) 
print "pAC is : " +str(pAC) 
print "pAG is : " +str(pAG) 
print "pCU is : " +str(pCU)
print "pCG is : " +str(pCG)
print "pCC is : " +str(pCC)
print "pGG is : " +str(pGG)
print "pGU is : " +str(pGU)
print "pUU is : " +str(pUU)

#concat both strings
random_rna = rna1 + rna2

length_both = len(random_rna) - random_rna.count("-")

length_rna1 = len(rna1) - rna1.count("-")
length_rna2 = len(rna2) - rna2.count("-")

qA = float(random_rna.count("A")) / length_both
qC = float(random_rna.count("C")) / length_both
qU = float(random_rna.count("U")) / length_both
qG = float(random_rna.count("G")) / length_both

print "========================"
print "qA is: " + str(qA)
print "qC is: " + str(qC)
print "qU is: " + str(qU)
print "qG is: " + str(qG)
print "Total sum of q is : "+str(qA + qC + qU + qG)

def match_scores(pab,qa,qb):
	numerator = float(pab)
	denom = float(qa*qb)
	return log(numerator/denom,2)

sAA = match_scores(pAA, qA, qA)
sAU = match_scores(pAU, qA, qU)
sAC = match_scores(pAC, qA, qC)
sAG = match_scores(pAG, qA, qG)
sCU = match_scores(pCU, qC, qU)
sCG = match_scores(pCG, qC, qG)
sCC = match_scores(pCC, qC, qC)
sGG = match_scores(pGG, qG, qG)
sGU = match_scores(pGU, qG, qU)
sUU = match_scores(pUU, qU, qU)

#match model = sAA * sAU * sAC * sAG * sCU * sCG * sCC * sGG * sGU * sUU
print "========================"
print "sAA is : " +str(sAA)
print "sAU is : " +str(sAU) 
print "sAC is : " +str(sAC) 
print "sAG is : " +str(sAG) 
print "sCU is : " +str(sCU)
print "sCG is : " +str(sCG)
print "sCC is : " +str(sCC)
print "sGG is : " +str(sGG)
print "sGU is : " +str(sGU)
print "sUU is : " +str(sUU)
print "========================="