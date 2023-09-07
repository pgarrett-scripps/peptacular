import timeit
import random

from peptacular import constants
from peptacular.fragment import *
from peptacular.sequence import *


# Function to generate a random sequence of amino acids
def generate_random_sequence(length: int) -> str:
    amino_acids = list(constants.AMINO_ACIDS)
    return ''.join(random.choice(amino_acids) for _ in range(length))


# Function to generate random modifications for a given sequence
def generate_random_modifications(sequence: str, num_mods: int) -> Dict[int, float]:
    mod_positions = random.sample(range(calculate_sequence_length(sequence)), num_mods)
    return {pos: round(random.uniform(0.5, 1.5), 4) for pos in mod_positions}


# Generate sample data for benchmarking
sample_sequence = generate_random_sequence(50)
sample_modifications = generate_random_modifications(sample_sequence, 5)
sample_span = (10, 20, 0)


# Benchmarking script
def benchmark_func(func, *args):
    timer = timeit.Timer(lambda: func(*args))
    n_loops, time_taken = timer.autorange()
    calls_per_second = n_loops / time_taken
    return calls_per_second


functions_to_benchmark = [
    (get_modifications, "PE(3.1415)PTIDE"),
    (add_modifications, "PEPTIDE", {1: 3.1415}),
    (strip_modifications, "PE(3.1415)PTIDE"),
    (build_fragments, "PE(3.1415)PTIDE", ['a', 'b', 'c', 'x', 'y', 'z'], [1, 2, 3], True, False),
    (fragment, "PE(3.1415)PTIDE", ['a', 'b', 'c', 'x', 'y', 'z'], [1, 2, 3], True, False),
    (split_sequence, "PE(3.1415)PTIDE"),
]

# Benchmark and store results
benchmark_results = {}
for func, *params in functions_to_benchmark:
    calls_per_sec = benchmark_func(func, *params)
    benchmark_results[func.__name__] = round(calls_per_sec)

print(benchmark_results)
