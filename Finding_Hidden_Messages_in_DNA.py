### AUTHOR: Zeynep KOKER --> zeyneepkkr@gmail.com
### Coursera: Finding Hidden Messages in DNA (Bioinformatics I)

import numpy as np
import matplotlib.pyplot as plt

# HammingDistance: Compare two equal length string, return count their difference.
def HammingDistance(Seq1,Seq2):
    count = 0
    for i in range(0,len(Seq1)):
        if Seq1[i] != Seq2[i]:
            count += 1
    return count

# input_file = open('dataset_30278_3.txt', 'r')
# file = input_file.read().splitlines()
# Seq1 = file[0]
# Seq2 = file[1]
# input_file.close()
# print(HammingDistance(Seq1,Seq2))

def ApproximatePatternMatching(Pattern,Text,d):
    appears = []
    m = len(Text)
    k = len(Pattern)
    for i in range(0,m-k+1):
        if HammingDistance(Text[i:i+k],Pattern) <= d:
            appears.append(i)
    return appears

# input_file = open('dataset_30278_4.txt', 'r')
# file = input_file.read().splitlines()
# Pattern = file[0]
# Text = file[1]
# d = int(file[2])
# result = str(ApproximatePatternMatching(Pattern,Text,d))
# result = result.replace(","," ")
# result = result.replace("[","")
# result = result.replace("]","")
# print(result)

def ApproximatePatternCount(Pattern,Text,d):
    count = 0
    m = len(Text)
    k = len(Pattern)
    for i in range(0,m-k+1):
        if HammingDistance(Text[i:i+k],Pattern) <= d:
            count += 1
    return count

# input_file = open('dataset_30278_6.txt', 'r')
# file = input_file.read().splitlines()
# Pattern = file[0]
# Text = file[1]
# d = int(file[2])
# print(ApproximatePatternCount(Pattern,Text,d))


def ImmediateNeighbors(Pattern):
    neighborhood = []
    for i in range(0, len(Pattern)):
        symbol = Pattern[i]
        for j in "ACGT":
            if j != symbol:
                neighbor = Pattern[:i] + j + Pattern[i+1:]
                neighborhood.append(neighbor)
    return neighborhood

def FrequentWordsMismatch(Text, k, d):
    freqMap = {}
    patterns = []
    maxPatterns = []
    n = len(Text)
    for i in range(0,n-k+1):
        Pattern = Text[i:i+k]
        patterns.append(Pattern)
    for i in range(0,len(patterns)):
        count = 0
        for j in range(0,len(patterns)):
            if (HammingDistance(patterns[i],patterns[j]) <= d):
                count += 1
        freqMap[patterns[i]] = count
    m = max(freqMap.values())
    print(freqMap)
    for key in freqMap:
        if freqMap[key] == m:
            maxPatterns.append(key)
    return maxPatterns


all_possible_seq = []

def printAllKLength(set, k):
    n = len(set)
    printAllKLengthRec(set, "", n, k)

def printAllKLengthRec(set, prefix, n, k):
    global all_possible_seq
    if (k == 0) :
        all_possible_seq.append(prefix)
        return
    else:
        for i in range(n):
            newPrefix = prefix + set[i]
            printAllKLengthRec(set, newPrefix, n, k - 1)


# Solve the Frequent Words with Mismatches Problem

def FrequentWordsMismatchAllPossible(Text, k, d):
    freqMap = {}
    patterns = []
    maxPatterns = []
    set = ['A', 'C','G', 'T']
    n = len(Text)
    for i in range(0,n-k+1):
        Pattern = Text[i:i+k]
        patterns.append(Pattern)
    printAllKLength(set, k)
    for i in range(0,len(all_possible_seq)):
        count = 0
        for j in range(0,len(patterns)):
            if (HammingDistance(all_possible_seq[i],patterns[j]) <= d):
                count += 1
        freqMap[all_possible_seq[i]] = count
    m = max(freqMap.values())
    for key in freqMap:
        if freqMap[key] == m:
            maxPatterns.append(key)
    return maxPatterns


# input_file = open('dataset_30278_9.txt', 'r')
# file = input_file.read().splitlines()
# Text = file[0]
# values = file[1].split()
# k = int(values[0])
# d = int(values[1])
# print(FrequentWordsMismatchAllPossible(Text,k,d))


def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern


def Reverse(Pattern):
    return Pattern[::-1]



def Complement(Pattern):
    comp = ""
    for i in range(len(Pattern)):
        if Pattern[i] == 'A':
            comp += 'T'
        elif Pattern[i] == 'T':
            comp += 'A'
        elif Pattern[i] == 'G':
            comp += 'C'
        elif Pattern[i] == 'C':
            comp += 'G'
    return comp




def FrequentWordswithMismatchesandReverseComplements(Text, k, d):
    freqMap = {}
    patterns = []
    maxPatterns = []
    set = ['A', 'C','G', 'T']
    n = len(Text)
    for i in range(0,n-k+1):
        Pattern = Text[i:i+k]
        patterns.append(Pattern)
        patterns.append(ReverseComplement(Pattern))
    printAllKLength(set, k)
    for i in range(0,len(all_possible_seq)):
        count = 0
        for j in range(0,len(patterns)):
            if (HammingDistance(all_possible_seq[i],patterns[j]) <= d):
                count += 1
        freqMap[all_possible_seq[i]] = count
    m = max(freqMap.values())
    for key in freqMap:
        if freqMap[key] == m:
            maxPatterns.append(key)
    return maxPatterns

# input_file = open('dataset_30278_10.txt', 'r')
# file = input_file.read().splitlines()
# Text = file[0]
# values = file[1].split()
# k = int(values[0])
# d = int(values[1])
# print(FrequentWordswithMismatchesandReverseComplements(Text,k,d))

def SkewArray(Genome):
    skew = []
    count = 0
    skew.append(count)
    for i in Genome:
        if i == 'G':
            count = count + 1
            skew.append(count)
        elif i == 'C':
            count = count - 1
            skew.append(count)
        else:
            skew.append(count)
    return skew


def MinimumSkew(Genome):
    positions  = []
    skew_array = SkewArray(Genome)
    min_value = min(skew_array)
    for i in range(len(skew_array)):
        if min_value == skew_array[i]:
            positions.append(i)
    return positions

# input_file = open('Salmonella_enterica.txt', 'r')
# file = input_file.read()
# Text = file
# skew_plt = SkewArray(Text)
#
# y_axis = []
# for i in range(0,len(skew_plt)):
#     y_axis.append(i)
#
# plt.plot(y_axis, skew_plt)
# plt.show()
#
# print(MinimumSkew(Text))


def PossibleMismatch(Text, d):
    patterns = []
    set = ['A', 'C', 'G', 'T']
    n = len(Text)
    printAllKLength(set, n)
    for i in range(0,len(all_possible_seq)):
        if (HammingDistance(all_possible_seq[i],Text) == d):
            patterns.append(all_possible_seq[i])
    return patterns

# input_file = open('dataset_30282_4.txt', 'r')
# file = input_file.read().splitlines()
# Text = file[0]
# d = int(file[1])
#
# result  = str(PossibleMismatch(Text,d))
# result = result.replace(","," ")
# result = result.replace("[","")
# result = result.replace("]","")
# result = result.replace("'","")
# print(result)
# with open("final.txt", "w") as out:
#     out.write(result)
# out.close()

# Week 3
# 1.2 Motif Finding Is More Difficult Than You Think

def MotifEnumeration(Dna, k, d):
    All_Patterns = []
    Final_Array = {}
    maxPatterns = []
    for seq in Dna:
        Seq_Patterns = []
        for i in range(0,len(seq)-k+1):
            Pattern = seq[i:i+k]
            global all_possible_seq
            all_possible_seq = []
            kmer_neighbors = PossibleMismatch(Pattern, d)
            kmer_neighbors.append(Pattern)
            for kmer in kmer_neighbors:
                if kmer not in Seq_Patterns:
                    Seq_Patterns.append(kmer)
        All_Patterns.append(Seq_Patterns)
    for seq_patterns in All_Patterns:
        for item in seq_patterns:
            if item not in Final_Array:
                Final_Array[item] = 1
            else:
                Final_Array[item] += 1
    m = max(Final_Array.values())
    for key in Final_Array:
        if Final_Array[key] == m:
            maxPatterns.append(key)
    return maxPatterns

# final = MotifEnumeration(['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT'], 3, 1)
# print(final)

# input_file = open('dataset_30302_8.txt', 'r')
# file = input_file.read().splitlines()
# k = int((file[0]).split()[0])
# d = int((file[0]).split()[1])
# Dna = (file[1]).split()
# print(k,d)
# print(Dna)
#
# result = str(MotifEnumeration(Dna,k,d))
# result = result.replace(",","")
# result = result.replace("[","")
# result = result.replace("]","")
# result = result.replace("'","")
# print(result)
# with open("final.txt", "w") as out:
#     out.write(result)
# out.close()

def countMotifPercent(motifs):
    count = {}
    columns = []
    for i in range(len(motifs[0])):
        columns.append([motif[i] for motif in motifs])
    for i in range(len(columns)):
        count[i] = {'A': columns[i].count('A')/len(columns[i]), 'C': columns[i].count('C')/len(columns[i]), 'G': columns[i].count('G')/len(columns[i]), 'T': columns[i].count('T')/len(columns[i])}
    return count

import math
def motifEntropy(motifs):
    entropy = 0
    percents = countMotifPercent(motifs)
    for i in range(len(percents)):
        for nucleotide in percents[i]:
            if percents[i][nucleotide] != 0:
                entropy += percents[i][nucleotide] * math.log2(percents[i][nucleotide])
    return -entropy

#motifEntropy(["TCGGGGGTTTTT","CCGGTGACTTAC","ACGGGGATTTTC","TTGGGGACTTTT","AAGGGGACTTCC","TTGGGGACTTCC","TCGGGGATTCAT","TCGGGGATTCCT","TAGGGGAACTAC","TCGGGTATAACC"])
# if u run this function you get ===> 9.916290005356972

def Pr(Text, Profile):
    n = len(Text)
    multiply  = 1
    for i in range(n):
        multiply = Profile[Text[i]][i] * multiply
    return multiply

def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def Score(Motifs):
    cons_seq = Consensus(Motifs)
    count = 0
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            if cons_seq[j] != Motifs[i][j]:
                count += 1
    return count

def Profile(Motifs):
    count = Count(Motifs)
    profile = {}
    symbol_str = "ACGT"
    for i in range(len(symbol_str)):
        profile[symbol_str[i]] = []
        for j in range(len(count[symbol_str[i]])):
            profile[symbol_str[i]].append(0)
            profile[symbol_str[i]][j] = count[symbol_str[i]][j]/len(Motifs)
    return profile

def ProfileMostProbableKmer(text, k, profile):
    max_value = -0.1
    max_str = ""
    for i in range(len(text)-k+1):
        if max_value < Pr(text[i:i+k],profile):
            max_value = Pr(text[i:i+k],profile)
            max_str = text[i:i+k]
    return max_str


# file = open("dataset_30305_3.txt", "r")
#
# text = file.readline().strip()
# k = int(file.readline().strip())
# profile = {}
# for i in range(4):
#     profile['ACGT'[i]] = [float(item) for item in file.readline().strip().split()]
# print(ProfileMostProbableKmer(text, k, profile))

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# print(GreedyMotifSearch(["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"],3,5))

# input_file = open('dataset_30306_9.txt', 'r')
# file = input_file.read().splitlines()
# k = int((file[0]).split()[0])
# d = int((file[0]).split()[1])
# Dna = (file[1]).split()
# print(k,d)
# print(Dna)
#
# result = str(GreedyMotifSearch(Dna,k,d))
# result = result.replace(",","")
# result = result.replace("[","")
# result = result.replace("]","")
# result = result.replace("'","")
# print(result)
# with open("final.txt", "w") as out:
#     out.write(result)
# out.close()

def AllPossibleDna(n, curr, ways):
    e = ['A', 'C', 'G', 'T']
    if len(curr) == n:
        ways.append(''.join(curr))
        return
    for element in e:
        AllPossibleDna(n, list(curr) + [element], ways)

def Mutation_Probability(d):
    all_possible = []
    AllPossibleDna(3,[],all_possible)
    file = open('Codon_Table.txt', 'r')
    codon_table = file.readlines()
    results = []
    codon_table_array = []
    for c_t in codon_table:
        line = c_t.split()
        codon_table_array.append(line)

    for c_t in codon_table_array:
        each_codon = []
        each_codon.append("Wild_Type:" + " " + c_t[0] + " " + "=" + " " + c_t[2])
        for each_possible in all_possible:
            if d == HammingDistance(each_possible,c_t[0]):
                for mut_control in codon_table_array:
                    if mut_control[0] == each_possible:
                        each_codon.append(each_possible + " " + "=" + " " + mut_control[2])
        results.append(each_codon)
    return results

# results_single_mutation = Mutation_Probability(1)
# print(results_single_mutation)
#
# results_double_mutation = Mutation_Probability(2)
# print(results_double_mutation)
#
# results_triple_mutation = Mutation_Probability(3)
# print(results_triple_mutation)


# Each Sequences in Dna Array compared with Pattern and returns distance of pattern with windows
def DistanceBetweenPatternAndStrings(Pattern, DnaArray):
    k = len(Pattern)
    Distance = 0
    for each_seq in DnaArray:
        hamming_distance = float("inf")
        for i in range(0,len(each_seq)-k+1):
            if hamming_distance > HammingDistance(Pattern, each_seq[i:i+k]):
                hamming_distance = HammingDistance(Pattern, each_seq[i:i+k])
        Distance = Distance + hamming_distance
    return Distance

# print(DistanceBetweenPatternAndStrings("AAC", ["AAAAA", "AAAAC", "AAAAG"]))
# print(DistanceBetweenPatternAndStrings("AAA", ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"]))
# print(DistanceBetweenPatternAndStrings("GAC", ["AAATTGAC", "GCGCAAAA", "CCGTAGCG"]))

# file = open("dataset_30312_1.txt", "r")
# lines = file.read().splitlines()
# Pattern = lines[0]
# Dna = (lines[1].split())
# result = str(DistanceBetweenPatternAndStrings(Pattern, Dna))
# print(result)
# with open("final.txt", "w") as out:
#     out.write(result)
# out.close()


def MedianString(k, DnaArray):
    Distance = float("inf")
    all_possible = []
    AllPossibleDna(k,[],all_possible)
    median_seq = []
    for each_possible in all_possible:
        if Distance >= DistanceBetweenPatternAndStrings(each_possible, DnaArray):
            Distance = DistanceBetweenPatternAndStrings(each_possible, DnaArray)
            median_seq.append(each_possible)
    return median_seq

# print(MedianString(7,["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC", "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC", "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"]))


# file = open("dataset_30304_9.txt", "r")
# lines = file.read().splitlines() # read in the input from STDIN
# k = int(lines[0])
# Dna = []
# for i in range(1, len(lines)):
#     Dna.append(lines[i])
# print(MedianString(k, Dna))

def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = Count(Motifs)
    for symbol in "ACGT":
        for i in range(k):
            count[symbol][i] = count[symbol][i] + 1
    return count

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = CountWithPseudocounts(Motifs)
    sum_value = sum(profile[item][0] for item in profile)
    for symbol in "ACGT":
        for i in range(k):
            profile[symbol][i] = profile[symbol][i]/sum_value
    return profile

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# file = open("dataset_30306_9.txt", "r")
# lines = file.read().splitlines()
# k = int((lines[0].split())[0])
# t = int((lines[0].split())[1])
# Dna = (lines[1].split())
# result = str(GreedyMotifSearchWithPseudocounts(Dna, k, t))
# result = result.replace(",","")
# result = result.replace("[","")
# result = result.replace("]","")
# result = result.replace("'","")
# print(result)
# with open("final.txt", "w") as out:
#     out.write(result)
# out.close()


def Motifs(Profile, Dna):
    motifs = []
    for i in Dna:
        motifs.append(ProfileMostProbableKmer(i, len(Profile['A']), Profile))
    return motifs

# print(Motifs({'A': [0, 0.25, 0.5], 'C': [0, 0, 0], 'G': [0.5, 0.5, 0], 'T':[0.5, 0.25, 0.5]},["TGACGTTC","TAAGAGTT","GGACGAAA","CTGTTCGC"]))

import random

def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(t):
        start = random.randint(0, len(Dna[i]) - k)
        motifs.append(Dna[i][start:start + k])
    return motifs

def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    # while True:
    for i in range(1000):
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        # else:
        #     return BestMotifs
    return BestMotifs

# print(RandomizedMotifSearch(["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"], 8,5))

# file = open("dataset_30307_5.txt", "r")
# lines = file.read().splitlines()
# k = int((lines[0].split())[0])
# t = int((lines[0].split())[1])
# Dna = (lines[1].split())
# result = str(RandomizedMotifSearch(Dna, k, t))
# result = result.replace(",","")
# result = result.replace("[","")
# result = result.replace("]","")
# result = result.replace("'","")
# print(result)
# with open("final.txt", "w") as out:
#     out.write(result)
# out.close()

def GibbsSampler(Dna, k, t, N):
    BestMotifs = []
    Motifs = RandomMotifs(Dna, k ,t)
    BestMotifs = Motifs.copy()
    for j in range(0,N):
        i = random.randint(0,t-1)
        del Motifs[i]
        profile = ProfileWithPseudocounts(Motifs)
        Motifs.insert(i, ProfileMostProbableKmer(Dna[i], k , profile))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


# print(GibbsSampler(["CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"],8,5,100))

# file = open("dataset_30307_5.txt", "r")
# lines = file.read().splitlines()
# k = int((lines[0].split())[0])
# t = int((lines[0].split())[1])
# Dna = (lines[1].split())
# result = str(RandomizedMotifSearch(Dna, k, t))
# result = result.replace(",","")
# result = result.replace("[","")
# result = result.replace("]","")
# result = result.replace("'","")
# print(result)
# with open("final.txt", "w") as out:
#     out.write(result)
# out.close()



