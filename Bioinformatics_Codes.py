
# Bioinformatics


import random


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for i in range(n-k+1):
        Pattern=Text[i:i+k]
        freq[Pattern]=freq[Pattern]+1
    return freq


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words

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

def PatternMatching(Pattern, Genome):
    positions = []
    k = len(Pattern)
    for i in range(len(Genome)-k+1):
        if Pattern == Genome[i:i+k]:
            positions.append(i)

    return positions



def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count


def SymbolArray_list(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i+(n//2)],symbol)
    return array


def SymbolArray(Genome, symbol):
    array = []
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        count = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
        array.append(count)
    return array


def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        array[i] = array[i-1]

        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array


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


def HammingDistance(p, q):
    count = 0
    for i in range(0,len(p)):
        if p[i] != q[i]:
            count += 1
    return count


def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    counts_mutations = []
    for i in range(len(Text)-len(Pattern)+1):
        temp = HammingDistance(Text[i:i+len(Pattern)], Pattern)
        counts_mutations.append(temp)
        if temp <= d:
            positions.append(i)
    return positions

def ApproximatePatternCount(Pattern, Text, d):
    positions = []
    counts_mutations = []
    for i in range(len(Text)-len(Pattern)+1):
        temp = HammingDistance(Text[i:i+len(Pattern)], Pattern)
        counts_mutations.append(temp)
        if temp <= d:
            positions.append(i)
    return len(positions)


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


def CountAminoacids(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ARNDCQEGHILKMFPSTWYV-":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def ConsensusAminoacids(Motifs):
    k = len(Motifs[0])
    count = CountAminoacids(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ARNDCQEGHILKMFPSTWYV-":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def CountAminoacidsInConsensus(Seq):
    count = {}
    k = len(Seq)
    for symbol in "ARNDCQEGHILKMFPSTWYV-":
        count[symbol] = []
        count[symbol].append(0)
        for i in range(k):
            if symbol == Seq[i]:
                count[symbol][0] += 1
    return count


def CountAminoacidsInConsensusFrequency(Seq):
    count = {}
    k = len(Seq)
    for symbol in "ARNDCQEGHILKMFPSTWYV-":
        count[symbol] = []
        count[symbol].append(0)
        for i in range(k):
            if symbol == Seq[i]:
                count[symbol][0] += 1
    for symbol in "ARNDCQEGHILKMFPSTWYV-":
        for j in range(len(count[symbol])):
            count[symbol][j] = count[symbol][j]/len(Seq)
    return count

def Score(Motifs):
    cons_seq = Consensus(Motifs)
    count = 0
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            if cons_seq[j] != Motifs[i][j]:
                count += 1
    return count

def Pr(Text, Profile):
    n = len(Text)
    multiply  = 1
    for i in range(n):
        multiply = Profile[Text[i]][i] * multiply
    return multiply

def ProfileMostProbableKmer(text, k, profile):
    max_value = -0.1
    max_str = ""
    for i in range(len(text)-k+1):
        if max_value < Pr(text[i:i+k],profile):
            max_value = Pr(text[i:i+k],profile)
            max_str = text[i:i+k]
    return max_str

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



def Motifs(Profile, Dna):
    motifs = []
    for i in Dna:
        motifs.append(ProfileMostProbableKmer(i, len(Profile['A']), Profile))
    return motifs


def RandomMotifs(Dna, k, t):
    motifs = []
    for i in range(t):
        start = random.randint(0, len(Dna[i]) - k)
        motifs.append(Dna[i][start:start + k])
    return motifs


def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs


def Normalize(Probabilities):
    total_probability = sum(Probabilities.values())
    normalized_probabilities = {key: value / total_probability for key, value in Probabilities.items()}

    return normalized_probabilities


def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    cumulative_prob = 0
    for kmer, prob in Probabilities.items():
        cumulative_prob += prob
        if p < cumulative_prob:
            return kmer
    return None

def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

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


