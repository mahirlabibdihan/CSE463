import random
from functools import reduce


def Profile(Motifs):
    count = Count(Motifs)
    number_of_motifs = len(Motifs)
    return {symbol: motif_profile_for_symbol(count[symbol], number_of_motifs) for symbol in "AGCT"}


def motif_profile_for_symbol(count_array_for_symbol, number_of_motifs):
    return list(map(lambda x: x / number_of_motifs, count_array_for_symbol))


def Motifs(Profile, Dna):
    return list(map(lambda text: ProfileMostProbableKmer(text, len(Profile[next(iter(Profile))]), Profile), Dna))


def RandomizedMotifSearch(Dna, k, t):
    random_motifs = RandomMotifs(Dna, k, t)
    return converge_to_optimum_motifs(random_motifs, Dna)


def HammingDistance(p, q):
    return sum(list(map(lambda x, y: int(x != y), list(p), list(q))))


def Score(Motifs):
    return sum(map(lambda motif: HammingDistance(motif, Consensus(Motifs)), Motifs))


def probability_of_generation(motif, profile_matrix):
    return reduce(lambda x, y: x * y, [profile_matrix[motif[i]][i] for i in range(len(motif))], 1)


def ProfileMostProbableKmer(text, k, profile):
    kmers = [text[iterator:iterator + k]
             for iterator in range(len(text) - k + 1)]
    probabilities = [probability_of_generation(
        kmer, profile) for kmer in kmers]
    return kmers[probabilities.index(max(probabilities))]


def Consensus(Motifs):
    count_matrix = Count(Motifs)
    symbols = "ACGT"
    string_length = len(Motifs[0])
    consensus_list = [consensus_symbol_at_index(
        count_matrix, i, symbols) for i in range(string_length)]
    return ''.join(consensus_list)


def consensus_symbol_at_index(count_matrix, index, symbols):
    count_to_symbols_tuple = {
        count_matrix[symbol][index]: symbol for symbol in symbols}
    max_count = max(count_to_symbols_tuple.keys())
    return count_to_symbols_tuple[max_count]


def add_pseudocount_toarray(motif_count_array):
    return list(map(lambda x: x + 1, motif_count_array))


def symbol_count_array_for_single_motif(symbol, motif):
    return [1 if current_symbol == symbol else 0 for current_symbol in motif]


def count_array(symbol, motifs):
    individual_count_arrays = [symbol_count_array_for_single_motif(
        symbol, motif) for motif in motifs]
    return [sum(elements) for elements in zip(*individual_count_arrays)]


def Count(Motifs):
    return {symbol: count_array(symbol, Motifs) for symbol in "ACGT"}


def CountWithPseudocounts(Motifs):
    motifs_count = Count(Motifs)
    return {key: add_pseudocount_toarray(value) for key, value in motifs_count.items()}


def ProfileWithPseudocounts(Motifs):
    motifs_pseudocounts = CountWithPseudocounts(Motifs)
    divisor = len(Motifs) + 4
    return {key: list(map(lambda x: x / divisor, value)) for key, value in motifs_pseudocounts.items()}


def converge_to_optimum_motifs(current_motifs, dna):
    newly_computed_motifs = Motifs(
        ProfileWithPseudocounts(current_motifs), dna)
    if Score(newly_computed_motifs) == Score(current_motifs):
        return current_motifs
    return converge_to_optimum_motifs(newly_computed_motifs, dna)


def RandomMotifs(Dna, k, t):
    return list(map(lambda text: select_random_motif(text, k), Dna))


def select_random_motif(dna_text, motif_length):
    position = random.randint(0, len(dna_text) - motif_length)
    return dna_text[position: position + motif_length]


def best_randomised_motifs(dna_array, motif_length, runs):
    best_of_current_run = RandomizedMotifSearch(
        dna_array, motif_length, len(dna_array))
    if runs == 1:
        return best_of_current_run
    best_of_other_runs = best_randomised_motifs(
        dna_array, motif_length, runs - 1)
    return best_of_current_run if Score(best_of_current_run) < Score(best_of_other_runs) else best_of_other_runs


def main():
    Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
           "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
           "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
           "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
           "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
    k = 8
    t = 5
    N = 100
    result = best_randomised_motifs(Dna, k, N)
    print(result)


if __name__ == "__main__":
    main()
