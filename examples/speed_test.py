import timeit
import random
from typing import Dict

import peptacular as pt


# Function to generate a random sequence of amino acids
def generate_random_sequence(length: int) -> str:
    amino_acids = list(pt.AMINO_ACIDS)
    return ''.join(random.choice(amino_acids) for _ in range(length))


# Function to generate random modifications for a given sequence
def generate_random_modifications(sequence: str, num_mods: int) -> Dict[int, float]:
    mod_positions = random.sample(range(pt.sequence_length(sequence)), num_mods)
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
    (pt.get_mods, "<13C>[Acetyl]-PE[3.1415][1]PTIDE-[Amide]"),
    (pt.add_mods, "<13C>[Acetyl]-PEPTIDE-[Amide]", {1: [pt.Mod(3.1415, 1)]}),
    (pt.strip_mods, "<13C>[Acetyl]-PE[3.1415]PTIDE-[Amide]"),
    (pt.mass, "[Acetyl]-PE[3.1415]PTIDE-[100]"),
    (pt.fragment, "[Acetyl]-PE[3.1415]PTIDE-[100]", ['a', 'b', 'c', 'x', 'y', 'z'], [1, 2, 3]),
    (pt.digest, "MASFRLFLLCLAGLVFVSEAGSVGAGEPKCPLMVKVLDAVRGSPAANVGVKVFK"*25, 'Trypsin'),
    (pt.split, "<13C>[Acetyl]-PE[3.1415]PTIDE-[Amide]"),
    (pt.span_to_sequence, "<13C>[Acetyl]-PE[3.1415]PTIDE-[Amide]", (3, 5, 0)),
    (pt.comp_mass, "<13C>[Acetyl]-PE[3.1415]PTIDE-[Oxidation]"),
    (pt.apply_static_mods, "<13C>[Acetyl]-PE[3.1415]PTIDE-[Amide]", {'P': 57.021464}),
]

# Benchmark and store results
benchmark_results = {}
for func, *params in functions_to_benchmark:
    calls_per_sec = benchmark_func(func, *params)
    benchmark_results[func.__name__] = round(calls_per_sec)

print(benchmark_results)



