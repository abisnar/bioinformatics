from scripts import *

def find_words(sequence):
    word_list = []
    for i in xrange(0,len(sequence)):
        word = sequence[i:i+3]
        if len(word) == 3:
            word_list.append(word)
        i +=1

    return word_list

seq = 'MAAALIRRLLRG'
print find_words(seq)

alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z']

def find_neighborhood_words(word, threshold):
    tally_scores = []

    for letter in alphabet:
        for letter2 in alphabet:
            for letter3 in alphabet[::-1]:
                first_score = BLOSUM62().scoring_matrix[letter,word[0]]
                print letter
                print first_score
                second_score = BLOSUM62().scoring_matrix[letter2,word[1]]
                print letter2
                print second_score
                third_score = BLOSUM62().scoring_matrix[letter3,word[2]]
                print letter3
                print third_score

                total_score = first_score + second_score + third_score
                print "------------"
                print total_score
                print "------------"
                if total_score >= threshold:
                    tally_scores.append((total_score,letter+letter2+letter3))
    return tally_scores

words = find_words(seq)
neighborhood_words = []

for word in words:
    neighbors = find_neighborhood_words(word,13)
    neighborhood_words.append((word, neighbors))

print neighborhood_words

