import random
import string


def generate_dna_string(length):
    return ''.join(random.choices("ATCG", k=length))


def mutate_sequence(sequence, mutation_rate):
    mutated_sequence = ""
    for base in sequence:
        if random.random() < mutation_rate:
            mutated_base = random.choice("ATCG".replace(base, ""))
            mutated_sequence += mutated_base
        else:
            mutated_sequence += base
    return mutated_sequence


def generate_multi_sequence_sets(num_sets, sequences_per_set, sequence_length, motif_length, mutation_rate):
    sets = []
    for _ in range(num_sets):
        sequence_set = []
        # Generate a random motif
        motif = generate_dna_string(motif_length)
        for _ in range(sequences_per_set):
            # Generate a random DNA string of length sequence_length
            random_sequence = generate_dna_string(sequence_length)
            # Randomly plant the motif in the random string
            start_index = random.randint(0, sequence_length - motif_length)
            end_index = start_index + motif_length
            modified_sequence = random_sequence[:start_index] + \
                motif + random_sequence[end_index:]
            # Mutate the sequence
            mutated_sequence = mutate_sequence(
                modified_sequence, mutation_rate)
            sequence_set.append(mutated_sequence)
        sets.append((sequence_set, motif))  # Include the motif in the set
    return sets


if __name__ == "__main__":
    t = 20
    n = 600
    motif_length = 8

    sets = generate_multi_sequence_sets(20, t, n, motif_length, 0.1)

    # Save each set in a FASTA file
    for i, (sequences, motif) in enumerate(sets):
        with open(f'set_{i+1}.fasta', 'w') as f:
            for j, sequence in enumerate(sequences):
                f.write(f'>Sequence_{j+1}\n')
                f.write(sequence + '\n')
            f.write(f'>Motif\n')
            f.write(motif + '\n')
