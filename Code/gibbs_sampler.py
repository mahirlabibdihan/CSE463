import random
from functools import reduce


def Normalize(Probabilities):
    sum_of_probabilities = sum(Probabilities.values())
    return {k: v/sum_of_probabilities for k, v in Probabilities.items()}


def key_vs_range(Probabilities):
    key_list = Probabilities.keys()
    value_list = map(lambda key: Probabilities[key], key_list)
    upper_bound_list = reduce(
        lambda result, element: result + [result[-1] + element], value_list, [0])[1:]
    lower_bound_list = [0] + upper_bound_list[:-1]
    lower_and_upper_bounds = list(
        map(lambda lower, upper: {'lower': lower, 'upper': upper}, lower_bound_list, upper_bound_list))
    return dict(zip(key_list, lower_and_upper_bounds))

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

def probability_of_generation(motif, profile_matrix):
    return reduce(lambda x, y: x * y, [profile_matrix[motif[i]][i] for i in range(len(motif))], 1)

def Pr(Text, Profile):
    return probability_of_generation(Text, Profile)

def WeightedDie(Probabilities):
    key_to_range = key_vs_range(Probabilities)
    random_fraction = random.uniform(0, 1)
    return next(key for key, value in key_to_range.items() if value['lower'] < random_fraction <= value['upper'])

def ProfileGeneratedString(Text, profile, k):
    text_length = len(Text)
    probabilities = {Text[index:index+k]: Pr(Text[index:index+k], profile)
                     for index in range(text_length - k + 1)}
    weighted_probabilities = Normalize(probabilities)
    return WeightedDie(weighted_probabilities)

def generate_new_motifs(dna_strings, kmer_length, current_motifs, index):
    temp_motifs = current_motifs[:index] + current_motifs[index+1:]
    profile_matrix = ProfileWithPseudocounts(temp_motifs)
    new_motif = ProfileGeneratedString(
        dna_strings[index], profile_matrix, kmer_length)
    return current_motifs[:index] + [new_motif] + current_motifs[index+1:]

def select_random_motif(dna_text, motif_length):
    position = random.randint(0, len(dna_text) - motif_length)
    return dna_text[position: position + motif_length]

def RandomMotifs(Dna, k, t):
    return list(map(lambda text: select_random_motif(text, k), Dna))

def GibbsSampler(Dna, k, t, N):
    random_indices = [random.randint(0, t-1) for i in range(N)]
    initial_value = RandomMotifs(Dna, k, t)
    output = list(reduce(
        lambda current_motifs, index: generate_new_motifs(
            Dna, k, current_motifs, index),
        random_indices,
        initial_value))
    return output


def main():
    Dna = [
        "GCCTTGAATGCTGAACATCAGATCTGGGGCTATATTCTTACCTTGATACATTTCAGAAGCAACTGAAATCGTAGGACCTTCCTTGCTTCTCTATTGGGTGAATGTTTCTCAGTCTTGGTGTGAGTCTCAGTGCCTACGTAGTTAAAGCTTACTGAAATGTTCCCTTTACAATTCTAGAGAGATATGTCCTTTATGTTGACATGTTCATGTTGACAGACTGCATCTGATTAAACAGCTGCCTGTGCAATGCCTCCAAGTGTGGATAAAAGAAAAATTAAACTCATAATCTTGGACAGCCATGTGTAGACTAGTTACATTGATCAAAGGGCAATAGAAATGATCCAGTGAGGATTTGTCTGAATTTCCCACAATTATTTAAAATCTACCTCAAATACCTGTTCATCTATAATGCCTCCCCTGAGGCCTTCATTCTGAATAGTACCTCTGTCTCTGTCCCCAAAGCACTAACTGATCCCTGTGATAGCGCACTTCCCAGCCAGGCTGATATGTAGACTTGGCTGCCTGTGTATCTTTTCCCCATAGACTGTGAGCTTCCTTTTATGAATAATAATTGTAGCTAGCATTTAGTAGGGTGCTCCTACCTGTTAAACTCTATGATGAGTGCTTTACATAGATTATATCATTTATTCACTAAACAGTCCTTTAAAATGGTGCTATATTCACTAAACAGTCCTTTAAAATGGTGCTATATTCACTAAACAGTCATTTAAAATGGTATTATTCTTCTTCATCTTACAGGTAAACAAACTAAGGCAAAAAAAAAAGTGAAATAATAAGTGCCAGTACACAGAGCTAGTAAGGAATAGGGTCTGCCAGGTCCCAAAAAGCATGCCATCACCTTTGCCCCATACTGCCTCTGGTACAGATAGAGGTAATGTCTTATTTATCACTGCCATCCACTGGACCCAGCTTAGTGCCTGACACACAGAGGGGCTCAGTCAATGCTGATTGGTTTGAGGTGGAGCAAAAATGCTTAGCAGGGTGAGCACCTTTGCTGTGATTGAGTATCTGATTCTCTATGAAGAGAAGGGGAGTCCTGAGCCAAACACATTCCTCTGGCTCCTGGCTGTCATCTTTATTTGCCCGGCTTCTTTGCTCTTCCTCCTTCCTAACTGCACCGTTTGGATTCAAAGCTGGAGCTTAATGCAGATAAAGGGAAAACAGAACTTTGAATGACCACTGTGGGACTAAGAGAGGAGAAACAAGAAATTTGACAGATGAGGAATAAAGTGAGGAGAAGAGAAAATGATTAAGCTTTATCACTTTAACTTAATATTTAACCTAATGAAAACAAAATCTTATTTGAAATTGGAAAAATCAATGTATTGATTGCTGGTTCATTGCCCTCTTCTTTATGATTTGACAGTCTGTGAATAATCTAATGGGTGTGGCTTAAAGACCTAGATCATGTGTGGAACTGGAATCGGGTGTTATTCAAGCAAAAAAAATAAATAAATACCTATGCAATACACCTGCTTT",
        "CACAGAGATACATGCTTCTAAAGAGCATAGGCCCTGTATTGGAGAATGGTGAGGAATAGTTGGATGAGAAGGCTAAAGCCAGGTTATCTGGGGCTATCTTCTTTCCGTCATACCTTTCCAAACTAGAAGCAACTAAAATCATAGGACCTCCTTTGCTTTCTATTGGGTGAATGTTTCTCAGTCTTGGTGTGAGTCTCAGTGCCTACATAATTAAAGCTTACTGAAATGTTCCCTTGGCAATTTTAGAGAAATATGTCCTCTATGCTAACATGTTCATATTGACAGACTATATCTAATTAAATAGTTGCCTGTATAGTTTTATGACGCCTCCAAGTGTGGATAAAAGAAAAATCAAACTCCTATGCTTGGAAAGCCATGAGGATACTAGCTACATTGATCAAAGGGTAATGGAAATGATCCAGTGAGGATTTTGCTCCATTTCTCACAATTATTTAAAATCTGCCACAAATGCCTGTTCCTCTACAATGCCTCCCCTGAGTCCTTCATTCTGAATAGCAACTCTGTCTCTGTTCCCAAAGCATTAACTGACCCCTGTGATAGGGCACTTCCCAGCCAGGCTGATATGTAGACTTGGTTTCCTGCGTGTCTTTTCACTATAAACTGTGAGCTTTCTTTCATGAATAATAATTGTAGCTAACATTTAGGAGGGTGCTCCTATCTGCAAAACTCCATGATGAGTGCTTTACAAGGATTATTTCATTTATTCACAAACAGTCCTTTAAAATGACGCCATTATTCATCTTCATCTCACAGATAAGGAAATTAAGGCAAAAAAAATTCAGTGAAATAATAAGTCCAAGTACACAGAACAAGTAAGGAACTGGGTTGGCCTGGTCCCACAAGCCATGCCATCACCGCTGCCCCATACTGCCTCTTGGCACAGGCAGAGGTCATGTCTTATTCATCACTGCTACTCCACTGGACCCAGCTTAGTGCCTGACACAGGGAGGTGTTCAGTCCATACATACTAGTTTTGAGTGGAGCAGAATTGCCAGGCACAGTGATCCCCTTTGCTGTGACTAATTGGGTGTCTGATTCTCTGCTGTAAAGACAAAGAGGCTGCTCTCTCTATCCTACTTCGTTTTTCCTGTTTTCTTACTCTTCTTCCTTCCTCGCCACCCCATGTGCATTCAAAGTTGCAGCTTAGTGCAGACAAAGGGAATACAAGGCAGAACGACCTATGTGGGATTGAAAGAAGAGAAACAAGGAATTCCAGAGGCCGGGGGGGGGTGGGAAGTGAGGAAAAGAGAAAGTGATTACAATTTATCACTTTAACTTAATATTTAAACTAATGAAAACAAAATCTTATCTAGAATTTGGAAGTCAATATTTTGATTGCTGGTTCAGTACCCTTTTATCTGTTTTGACAGTCTGGGAATAATCCAGTGGGTGTGGCTTAAAGACATAGATCACGTGTGGAATTGGAATTGGATGTTACACAAGCAAACAAAATAAATATCTGTGCAATATATCTGCTTT",
        "AGGTAGGGACTTGGAGAACTTTTCTGTCTAGCTAAAGGATTGTAAACACACCAATCAGTGCTCTGTGTCTAGCTAAAAGTTTGTAAACACAACAATCAGCACTCTGTAAAAATGCACTAATCAGTGCTCTGTGTCTAGCTAAAGGTTTGTAAACGCACCAATCAGCACTCTGTAAAAATGCACCAATCAGTGCTCTGTGTCTAGCTAAAGGTTTGTAAATGCACCAATCAGCACTCTGTAAAAACGGACCAATCAGCACTCTGTAAAATGGACCAATCAGCGCTCTGTAAAATGGACCAATCAGCAGGACATGGGGTGGGGGGGGCCAAATAAGGGAATAAAAGCTGGCCACCCAAGACAGCAGCAGCAACCCACTCAGGTCCCCTTCCATGCCATGGAAGCTTTCTTCTTTTGCTCTTCGCAATAAATCTTGCTGCTGCCCACTCTTTGGGTCCATGCTACCTTTATGAGCTGTAACACTCATGGCAAAGGTCTGCAGCTTCACTCCTGAAGCCAGCAAGGCCATGAACCCACCAGGAGGAACAAACAACTCTGGACGTGCCATGTTTAAGAGCTGTAACACTCACTGTGAAGGTCTGCAGCTTCACTCCTGAAGTCAGTGAGACCATGAAGCCACTGGGAGGAATGAACAACTCTGGACATGTCACCTTTAAGAGCTCTGACACTCACTGCGAAGGTCTGCAGCTTCTGGACACAATACTATTCATTCTCACCTTAAAGACGAGGAAACTAAGGCAAAGAACAGTCAACTAATAAGTCCAAGTATACAGAGCTGCTAAGGAATAGTCTGTCTGATCCCAAAGGCTGTGTCATAACCGCTTCCCTATACTGCCTCTCAGCAGAGGTAAGAGTCAAGTTTTATTTATCACTGCCACCCCATCAGCCCCAGCTTAGTGCCTGACACAGGGAGATGCTCAATCAATGCTGATTGTTATTGAGTGGACTAGAAATGCAAGGCACAGTGAGCCCCTTTGCTGTGACTGATGGGGTGTCTGATTTTCTGCTATAAAGAGGAGAGTGCTGTATCAAACACACTCCTCTGGCTCCTAGCTCTCTCTGTTCCACTTTGTTTATCCAATTTCCCTACTCCTCCTTCGTAACTGCACCATGTGGATTCAAAATTGCAGCTTAGTGCAGACAAAGGGAAAACGGAATTCTGAATGACCCCAAAGGGAAAACTGAACTCTGAATGACCCCTGTGGGTTTGAGAGAAGAGAAGCAGGAACTTGAGAGAGGAGGAAGAGAGAAAGTAATTAAAATGTATCGTTTTAACTTAATATTTAACCGAATGATAGCAAAATCTTATCTGAAATTGGAAAAGTCAAGGTTTTGAGTGCTGGTTCGGTGCCCATTTCTTTATGATTTGATAGTCTGAGAAGAATACGACGGGTGTGGCTTAAAAACCTAGATCACGTGTGTAGTTGGAATTGGGTGTTATATGAGCAAACAAAATAAATACCTGTGCAACATACCTGCTTT",
        "TAGTCAGAAAGACAAGTCTGGTATTTCCTTTCAGGACTCCCTTGAGTCATTAAAAAAAATCTTCCTATCTATCTATGTATCTATCATCCATCTAGCTTTGATTTTTTCCTCTTCTGTGCTTTATTAGTTAATTAGTACCCATTTCTGAAGAAGAAATAACATAAGATTATAGAAAATAATTTCTTTCATTGTAAGACTGAATAGAAAAAATTTTCTTTCATTATAAGACTGAGTAGAAAAAATAATACTTTGTTAGTCTCTGTGCCTCTATGTGCCATGAGGAAATTTGACTACTGGTTTTGACTGACTGAGTTATTTAATTAAGTAAAATAACTGGCTTAGTACTAATTATTGTTCTGTAGTATCAGAGAAAGTTGTTCTTCCTACTGGTTGAGCTCAGTAGTTCTTCATATTCTGAGCAAAAGGGCAGAGGTAGGATAGCTTTTCTGAGGTAGAGATAAGAACCTTGGGTAGGGAAGGAAGATTTATGAAATATTTAAAAAATTATTCTTCCTTCGCTTTGTTTTTAGACATAATGTTAAATTTATTTTGAAATTTAAAGCAACATAAAAGAACATGTGATTTTTCTACTTATTGAAAGAGAGAAAGGAAAAAAATATGAAACAGGGATGGAAAGAATCCTATGCCTGGTGAAGGTCAAGGGTTCTCATAACCTACAGAGAATTTGGGGTCAGCCTGTCCTATTGTATATTATGGCAAAGATAATCATCATCTCATTTGGGTCCATTTTCCTCTCCATCTCTGCTTAACTGAAGATCCCATGAGATATACTCACACTGAATCTAAATAGCCTATCTCAGGGCTTGAATCACATGTGGGCCACAGCAGGAATGGGAACATGGAATTTCTAAGTCCTATCTTACTTGTTATTGTTGCTATGTCTTTTTCTTAGTTTGCATCTGAGGCAACATCAGCTTTTTCAGACAGAATGGCTTTGGAATAGTAAAAAAGACACAGAAGCCCTAAAATATGTATGTATGTATATGTGTGTGTGCATGCGTGAGTACTTGTGTGTAAATTTTTCATTATCTATAGGTAAAAGCACACTTGGAATTAGCAATAGATGCAATTTGGGACTTAACTCTTTCAGTATGTCTTATTTCTAAGCAAAGTATTTAGTTTGGTTAGTAATTACTAAACACTGAGAACTAAATTGCAAACACCAAGAACTAAAATGTTCAAGTGGGAAATTACAGTTAAATACCATGGTAATGAATAAAAGGTACAAATCGTTTAAACTCTTATGTAAAATTTGATAAGATGTTTTACACAACTTTAATACATTGACAAGGTCTTGTGGAGAAAACAGTTCCAGATGGTAAATATACACAAGGGATTTAGTCAAACAATTTTTTGGCAAGAATATTATGAATTTTGTAATCGGTTGGCAGCCAATGAAATACAAAGATGAGTCTAGTTAATAATCTACAATTATTGGTTAAAGAAGTATATTAGTGCTAATTTCCCTCCGTTTGTCCT",
        "CTGTGGATGAGCAACTGATCACTGGAGGGAGTTTAGCTGCCCATAGGAGTTCATGGCTAATGACAATATCTGAATAAGGACAGGTGTGGAGCCCAGGTGCAGGAAGCAGGCGAAGGTCTTTCTGTGAGTCTCCTCTGAGGGAACTGGGTCTTTATACATAGTTACTGTTTCAGAATTGATCCTTCTGGAATCATCAGTCTTCACCAGTAGCTTGTTACATCTGGGGTTATCTCATAATTCAAACAAAGCTGACAAGTTGTAACAATGAGCACACACTGACTTCTGCAACAGGCGCTGTCCACTTCCCATCCGCACTCTACCGGCTTGCTCCTGGCCGCCTCCCACTCGCCTTCCTGGGTGGTCCCCCAGCAGTTATACCTACCTGGTTGTCGCCCCCTCTATCCTACCACAATTGCTCACTAGCGGTTTCCTGCGTACACAGCTTGTCTCCCTAACCAGAGTGGAGGTGCCTTGGGGACACAGCCAGGCTCAGACATTCACTCAGCTCATCATAGTGCCATCCCATCAATAACCCCTTCTGAGTGATCCTGGGTTAGTAAACCGAGTGTCCCTGAAATTCCACTACCGCTGATTCCCTCCAGCTGGGCAGAGGCAGCGAGCGCTGGCTGAAGCTTCCGGTGGGAAATGGGCAGTGCCTAGAAGAGAAGGAAACGATGCATGAGAAGGTTCCAGATGTCTATGAGGAACATGACGTGTCCTGTCCACTACTCTGCTTTTCCTCGTCCGCCTCCCCACCACTGGAGGAAACCTAGAAGCTGGTGCAGGAAATCCTCCTCTCAACAACCCAAGAACACTTTGCACAAGAGGGGTGCGCCCTCGGAGGTTGCTCTTCCCCAGAGGCCTCTCCTCGCTGGGGTTTCTTGAAGACAGATACTTGGACTCCTGCTGGGACCAGGCAGGCCACCCATCCTCAGGGGCAGTGACTGGTCACTCACCAGACCTCCCTGCATCCCCCTTCTCTCTCCTCCCCCAGCACGGGCTGAACCCCGCAGCCACAGATTCTGATCAGGATTAGGGTGTGGGTGCAAATCCAAGGTCCACCAAAATGGAAAAGAAGTAACCGATGGGAACACGTCTCCACCAAGACAGCGCTCAGGACTGGTTCTCCTCGTGGCTCCCAATTCAGTCCAGGAGAAGCAGAGATTTTGTCCCCATGGTGGGTCATCTGAAGAAGGCACCCCTGGTCAGGGCAGGCTTCTCAGACCCTGAGGCGCTGGCCATGGCCCCACTGAGACACAGGAAGGGCCGCGCCAGAGCACTGAAGACGCTTGGGGAAGGGAACCCACCTGGGACCCAGCCCCTGGTGGCTGCGGCTGCATCCCAGGTGGGCCCCCTCCCCGAGGCTCTTCAAGGCTCAAAGAGAAGCCAGTGTAGAAAAGCAAACAGGTCAGGCCCGGGAGGCGCCCTTTGGACCTTTTGCAATCCTGGCGCTCTTGCAGCCTGGGCTTCCTATAAATGGGGTGCGGGCGCCGGCCGCGC",
        "GCCACTGCCCTCCAGCCTGGGTGACAGGGCAAGACTCTAAAAAAAAAACACTCAAACAAACAAAATATCCCCCAAAAAGTAGGAGGCTGGTTACTTTCTCACAATATAACAAGAGGCCTGTAACCTGTAAGAATGAGGCAGTTCTTTGCTCACTGAGGTGAAATAGCCTCTGAGGTATATTGTTCATGAAAAAACGAAACAAAACGAAACCCAAGATTTAACTGAAGAGACCAGGAAGAATAGTATGTGCTATGTGCTGTCCACAGGGCACAGTAGTTCACACCAGCACTTTGTGAGGCTGCTGCGGGAGGATCACTTGAGCCCAGGAGTTCAAGACTGGACTGGGCAACATCGTGGGACCCCCATCTCCACAAAAATAAAAAAATTATCCGGGCATGGTGGCGGCCACCCGTAGTCCCGGCTACTTGGGTGGTTGAGCCAGGATGATCACTTGACCCCAGGAGGTTGAGGCTGCAGTGAGCTGTGATTGCACCACTGCAATTTAGCCTGAGTGACAGAATGAAAAAAAAATTTTTTTAAAGGAAAACACAAAAAGAATATGCTGTCAACAGGGATGGGAGGAAGACCACCTTTACTGCTATACACATTTGTACCTTTTAGATGTTGATCAATATGAATATATTATACACACAGACACACACACAGACACACACACACACACAAACAATACAATTTAATATCCTAAGAGGATATTGACATTAGACAGGTACAAAAGCTCTAGAAATGAGGACTTTCCTCAGTGATGACTTTTTTCACCACCAAAGTCACTCAGGCATCCTGACAAGGGTAAGTGAGGGGAGCCTCCTTGGAAAATAAACTCACTTGGATAGTGAACTCCTGCACATACCTCAAAGCCCATCTGAAATGTCCCCTCCTACAGGAAGTTTTCCCTGACCCTCCAAGAAGCAGAGTTCTATTTCACTGGGGAAAACATTTCTTCTTCTTCTTTTTTTTCCCTGCCCTGCACATGAGCTAGAAAACATTTCATGAAACTGGGAGTTTCTGTGCTGGGCTCTGTCCCTCCCCCATTCTACTTCCCCTCCCTCAGCATGGAAGCCTCTGGAAGTGGGGCTCTGACTCCCAGCCTACAGAGAGATTCCTAGGAAGTGTTCGACTGATAAACGCATGGCCAAAAGTGAACTGGGGATGAGGTCCAAGACATCTGCGGTGGGGGGTTCTCCAGACCTTAGTGTTCTTCCACTACAAAGTGGGTCCAACAGAGAAAGGTCTGTGTTCACCAGGTGGCCCTGACCCTGGGAGAGTCCAGGGCAGGGTGCAGCTGCATTCATGCTGCTGGGGAACATGCCCTCAGGTTACTCACCCCATGGACATGTTGGCCCCAGGGACTGAAAAGCTTAGGAAATGGTATTGAGAAATCTGGGGCAGCCCCAAAAGGGGAGAGGCCATGGGGAGAAGGGGGGGCTGAGTGGGGGAAAGGCAGGAGCCAGATAAAAAGCCAGCTCCAGCAGGCGCTGCTCA",
        "CAGATATAGAATTTGCAAATATTTTCTCCCATTCTGTAGACTGTCTTCACTTTTTTTTTTGAGACAGAGTCTTGTTCTGTCACCCAGGCTGGAGTGCAGTGGTGCGATCTAGGCTCACTGCAACCTCTGCCTCCCAGGTTCAAGCAATTCTCCTGCCTCAGCCTCCGAGTAGCTGGGATTAAAGCATGCGCCACCATGCCTGGCTAATTTTTGTATTTTTAGTAGAGATGGAATTTCACCACATTGGTCAGGCTGGTCTCAAACTCCTGACCTCATGATCCGCCCACCTCAGCCTCCCGAAGTGCTGGGATTACAGGCGTGAGCCACCACGCCTGGCCTGTCTTCATTTTATTAATACTGTCTTTTGATACACAAATGTTCATTCCGATGAAGCCCAATTTATCTAATTTTTCTTTGATTGCTGGTGATTTTGGTGTCATTATCTAAGAATCCATTGCCAAGTTTTACTTGGTTTTGATCTCTAGAAGACTGCCATCATACTGAAGGTAATCTTCTGAAACTTGCGGTTTTCTTTCAAAAATTATGTTTATAAGATGATCCATGTTCTTGCAGAGTTTATTCATTTTATGCTGTATAATATTCCATTATATCCACATACAATGCAGTATTGACCCTTCCTCCTGTTGATGGGCATTTGTCTTGTTTCTAGTTACTTTGCTATTATATCAGTGTCACCATGATTATCCAAAAGTAATTCTTTTGTACACTCTAATTTAAGAACAACTAACCCTTTTTAATGAATAAATCAACCTTGTATTGAGTTGCTACTAAGTTTCAGTTGACTAGTACCTGGGATACACACAGGTGCAGACATTTGACTGAGACATATTGATTTTTCTCATCTGCCTATTTAGGCTAATCACCAGACTATAAAACCATGAGAACCACTGCCATTGAGTATAGTCTGTGTCAGTCTACACTATAGCTTTAACTAGTTGTGTGATTTCTTGCAAAGAGCAATCAGAGAAGACACAATAAACACATTTACTGATTTCAGGCTGGAGAGCTTTTAAGCAATAGGGAGATGGCCACACACAAGGTGGAGAAAATTACTGTGAAAAGGAAGTACTTTCTTTAGAGCCCCACCTAAGCTAGGCTGCAGAAATGTCTACAATGGGTTTGAAAAAACTCAAAATGAGCCTTTCTGCAGTGTGAAAATCCTCCAAGATAAAGAGACAGATTGATGGTTCCTGCCGCCGCCCTGTCCTGCCCAGTTGCTGATTTCAGGAAATACTTTGGCAGGTTTGTGGGTCATAGAGTTGCCAGGTTTCTTGGGATTTGTAATAGAACATCACAAGAAAATCAAGTGTGAAGCAAGAGCTCAACTCTTAACAGGGGTATTGTTTGTGGTTTTGTTACTGGAAAAGATAGTGACCTTACCAGGGCCAAAGTTTGTAGACACAGGAATTACGAAATGGAGAAGGGGGAGAAGTGAGCTAGTGGCAGCATAAAAAGACCAGCAGATGCCCCACAGCACTG",
        "AGGGGTCTGGAGCCATCTGTGAGGGATCAGGGCCCTTTCAGCCTTGGCTAGGGAGCAGGGGTCCTGGAACTTCATCCTGGCCCATAGCTGAGTCTGCCCATAATTCTTTTCTGACTCACTAGGCAAATCTCACACAGAAATGGGGCAGCTTTGGGAGTGGGCCCAGGAAGTACTGAGGATAGCAGGTGAGATCCCAGGAAGAGATGGATGTGGGGCCGAGACACTGGAGAGAGAAACAGGACTGTCAGATAAAGGGCGTCTGTGACTCCTAGATCTCATTATGCCTACTACCATAACCTACCCCCAATTCCTAATATTCTCCTACCCTAGAGGGGGGGAAATTGTCAGAAATTTGGCTGCAACACTAGCAACACTACTCAGTACTTGAAATGCATTTTTGCATTTTTTTCATTCAACAAATATTTCTGGAACAACTCTTATATGCCAGGCACTATTTTAGGAGTCAGGGATATATAATGGTAAACAAGACAGGCAAAACAAAGCAAAGCAACAACAACCATCACCAGATAAGTAGACAGATGAAAGAATTTCAAGTTTTAGTAAGTAAAATAAAACAAGCAAGGGTCTGAAATGGCTAGATAAGGTGGTCAAGAAAGGCTTCATTGAGAAGGTAGCATTTAAGCAGGAGTCAGCTAGAAATATTGTGAAATTCCAGTTACAGTTCTATTTGTTCTGGGTTGGTTAAATAAAGCTTTTTCCCCCAAGGTGGAAACTACCAAGAAAGACTAATTACTAGTAGTGGTGGTGCTCTCTGGAAGAGAGACACCTCCTGTTTCTGCCTCATTACTGTCAACCCTTCACTTCCAGGCACTTTTTGCAAAGCCCTTTGCCAGTCAGGGAAGGCGAGAGGCTGGGCATGGGGCTTGGACATTTGACAACAGTGAGACATTATTGTCCCCAGACTCACTAGCCCAAGGGTAAAGCTGAAGAGGCTTGGGCATGCCCCAGAAAGGCCCCTGATGAAGCTTGGAAAAAGCTGTTCTCTGAGTATTTCTAAGTAAGTTTATCTGTGTGTGTGGTTACTAAAAGTAGTAAGTATTGCTGTCTCTAGCTGCCTTAGAGCAGGGCTTGACACAGTACACAGCAATATTAGTTCCCTCCTTTTCTCACCTCCCCCATTGTGGAGATAAACTCAATCACAAAAGGTGATCCTCAGTCTACTCACTTCCCTGACTTATGGATGCCTGGACCCATTGCCAGTGTGAGAGTCACAGCTGGACGTCAGCAGTGTAGCCCAGTTACTGCTTGAAAATTGCTGAAGGGGGTTGGGGGGCAGCTGCCGGGAAAAAGGAGTCTTGGATTCAGATTTCTGTCCAGACCCTGACCTTATTTGCAGTGATGTAATCAGCCAATATTGGCTTAGTCCTGGGAGACAGCACATTCCCAGTAGAGTTGGAGGTGGGGGTGGTGCTGCTGCCAACTCTATATAGGGAGTTCAACTGGTCACCCAGAGCTGTCCTGTGGCCTCTGCAGCTCAGC",
        "AGGCAGGAATAGGTTTGGCCTGTTGCATGAACAGTGGGTCCAGCTCCTAGCAAACTGTTTATTGAATGAAAGAAGAATGAATGCCTTGGGTCTAGGGTTGTGCTGGGCGCTTTCTTAAGTTTTCTTTCCCGGGTACCTCCCCAGAACTGGCATGCAGGTATTATTAAACCCATTACACAAGTGAAACTGGCCCAGAGACAGAAAAGTCCCTGGTCCAAGACCACACAGGAGTGAGGGGTGGAGGAACCCTCCTCCCATTGAGTTCTGGCTTTCCTATACTGAAAGCCCCTTCCTCTCCTGCAGTAAGGTAGGTGGAACCGCTGTCCCGCCTTGTTGGTGAATGTCGTTGCTAGACTTCAGACACATACAGGCTGGTCTGCTGAAAATCAGAGATGTCCACCTGCGCCCTATTCGAGGTCTCCGGCGTCTTCTTTGGCGTCGTCTTTGCCCTTTCAGAAGCGTCTGCACATTTTTCCAGGTGTCATTTCTCCAACTTGAACACAGGGAGCGCACTGGGCACGCGGGCACGTGGCTGTCCCCAGGGGCCTGGCTTGGGTCTCGCCCCTGGGCCGGGGCGCACTCGCGGGCGGGACATCTGGGGGCGCCCACGCGCTCTGGGACGAGTGTCGCTGGCCAGGCCCGGACTGAGGAAAGGCGAGTGAGACACTACTCGCCTGGGGTGCAAAATTTAAGGGAGTGAAAAAAAAAAAAAAAGAAAGAAACCAAAACCACCTCGAGTCACCAAAATAAACATTTTAATGCAGTATTTTTTAAAAAATCAACAGGAATCCTCCAAAGCCCACTATGAACAAAATAGCAAAATGGTAGAGAAAGGATCTGTGCCGCTGCGTCGGGCCTGTGGGGCGCCTCCGGGGGTCTGAAACTGGAGGAGACTCGGGGCTGTAGGGCGCGCGGATCTGGGGCGCGCCCTCGGTCCCGGCGCGCCCAGGGCCTCCCGCGCGGGGCCCGGCACAGGGAGGCGGGGAGGCGGGCGGGGCGGGGCGGGGCCGGGCGGCACCTCCCTCCCCTGCAAGCTTTCCCTCCCTCTCCTGGGCCTCTCCCGGGCGCAGAGTCCCTTCCTAGGCCAGATCCGCGCCGCCTTTTCCCGCGGCCCGCACGGGGCCCAGCTGACGGGCCGCGTTGTTTACGGGCCCGAGCAGCCCTCTCTCCCGCCGCCCGCCCGCCACCCGCCAGCCCAGGTGCCCGCCCGCCAGTCAGCTAGTCCGTCGGTCCGCGCGTCCCTCTGTCCCGGAGCCCGCAGATCGCGACCCAGAGCGCGCGGGGCCGAGAGCCGAGAGACAGTCCCGGGCGCAGCGCGGAGCTCCGGGCCCCGAGATCCTGGGACGGGGCCCGGGCCGCAGCGGCCGGGGGGTCGGGGCCACCACCGCAAGGGCCTCCGCTCAGTATTTGTAGCTGGCGAAGCCGCGCGCGCCCTTCCCGGGGCTGCCTCTGGGCCCTCCCCGGCAGGGGGGCTGCGGCCCGCGGGTCGCGGGCGTGGAA",
        "CTCCCACCTGAGTCTGCAGTCATGACAGATGGAGCCTTTCAGGTGTCCTAAGGGAAGATGGGGTTCACCTCACAATGCTAAAGGTGAAAAGGGGAGTGCCGACCCTGCCCTTTGTTTTGGCATTTTTATTTCCATACTTTGTGGTTTGTATGTGCCCAGTTAAACGGATGTACAGATTAGATAATCTGGCTTTAACTCCAACCTTCACGAGATGCACACATGCCTTAAAAAGATTCTTACAAAGAAACTAGAATTAGGAGACAATTCTGACAACACAAAAACAAGTGAACAAAAATCAAATGACTGACCTGACGTGTCCACTAGGCTTGATATGAGCTTGGTCAACTAAAACATGAGACAGAATGAATGACCATCTGGTCAGGTTGGATAAGAAGTAACTTTTGGCAAAATCTCTGAAGTCTGGATTTACACTTCGTGGAAGGAATCTTGACCTTTCCCCAAGTATATTATGCCTTGAGATAATAATCCTGTATAAGGCCGCATAAGCAAAACGACATTTGAAAATGAGTTTCATCTTTAGCTCATCCTTTTACAAATTTCTATTCCTTGTGCATGATTGAAGTGTACACATTTTTTAGTTAAAATATAAGCATATAGAAGGCACAAGTTCAATATTTTTCATGATAAGGATGTGTAATTAAACTGGTTAGGAGACTCATCGACTAGAGGGGGACATGGTGGCCCCAGGCTGTAAGAACAGGCCACACCGTCCACTGGGCCGCTTGCTTTGTGCTAAGATGACACTTTGTTCTGAGCCTCACAGTGTCTTGACCATGTTCCTGGAACCTTCTTGTTGGAGGGAGTTCATCTTCCCCTATGACTCTGTCCCTAGTCTAAGGTGTCCCACAGGAAGCTTGAGGGCGGGAAGTTTTCCAGCCCAGGAGCCTGAGCTCAGCGGGGCAGGAAGAGGGAGCAGCTCCTCCGTGGGGGACCTTTGAGAGCCCAGGAGCAGGATTTCGAGGGACACCTGGTGGGGAGCAAAAGGTGCTGAGTCTGTCTTTGACCTTGAGCCCAGCTTGTTTCTCCTGCATCCTCCCCCAAAAGGGGCTTTGCCTGTCATTCTGCAGTTCTAGTGTGGGGTCTGGGCGCAGTTCTTTTCCCTCTCCAGCCTCGGAGTCTTCCTCTGTGGACTGCGCAGATAGGACTGGTGGCACGGACCAGCTCTGCAGCCCTGGAGTCAGGAGCAGAGCCCCCCGGCTCCCAGCCCGCCGTAGCCGCTCCTGGCACCGAGCGAGCCGCGATGACAATGGCTGCATTGTGCTTCATGTCCCTTCCCATCAACATTTCTGTGCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAAAATGAGGTGATCAGTGGGACGAGTAAGGAAGGGGGGTTGGGAGAGGGGCGATTGGGCAACCCGGCTGCACAAACACGGGAGGTCAAAGATTGCGCCCAGCCCGCCCAGGCCGGGAATGGAATAAAGGGACGCGGGGCGCCGGAGGCT"]

    k = 8
    t = 10
    
    N = 100
    result = GibbsSampler(Dna, k, t, N)
    print(result)


if __name__ == "__main__":
    main()