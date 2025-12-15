"""
Mass Calculation Benchmark: Peptacular vs Pyteomics
Testing with random sequences
"""

import random
import statistics
import time

from pyteomics import mass

import peptacular as pt


def generate_random_proforma(length: int = 20, mod_probability: float = 0.2) -> str:
    """Generate a random ProForma sequence with modifications."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    mods = ["[+15.995]", "[+57.021]", "[+79.966]", "[+42.011]"]
    
    sequence_parts: list[str] = []
    
    # Random N-term modification (20% chance)
    if random.random() < mod_probability:
        sequence_parts.append(random.choice(mods))
        sequence_parts.append("-")
    
    # Build sequence with random modifications
    for _ in range(length):
        aa = random.choice(amino_acids)
        sequence_parts.append(aa)
        if random.random() < mod_probability:
            sequence_parts.append(random.choice(mods))
    
    # Random C-term modification (20% chance)
    if random.random() < mod_probability:
        sequence_parts.append("-")
        sequence_parts.append(random.choice(mods))
        
    
    return "".join(sequence_parts)


def benchmark_mass_calculation():
    """Benchmark: Calculate mass of random ProForma peptides"""
    
    iterations = 10
    num_sequences = 10000
    
    # Generate random sequences
    print(f"Generating {num_sequences} random ProForma sequences...")
    sequences = [generate_random_proforma(length=random.randint(10, 30)) for _ in range(num_sequences)]
    
    print(f"Example sequences:")
    for i, seq in enumerate(sequences[:3]):
        print(f"  {i+1}. {seq}")
    
    # Benchmark Peptacular
    peptacular_times = []
    all_masses_pept = []
    
    print(f"\nBenchmarking Peptacular (sequential - {iterations} iterations per sequence)...")
    for seq in sequences:
        for _ in range(iterations):
            start = time.perf_counter()
            result_pept = pt.mass(seq)
            end = time.perf_counter()
            peptacular_times.append((end - start) * 1_000_000)  # microseconds
        
        all_masses_pept.append(float(pt.mass(seq)))
    
    # Benchmark Peptacular with parallel processing
    peptacular_parallel_times = []
    all_masses_pept_parallel = []
    
    _ = pt.mass(sequences) 

    print(f"Benchmarking Peptacular (parallel - {iterations} iterations)...")
    for _ in range(iterations):
        start = time.perf_counter()
        results_parallel = pt.mass(sequences)  # Pass all sequences at once
        end = time.perf_counter()
        peptacular_parallel_times.append((end - start) * 1_000_000)  # microseconds
    
    all_masses_pept_parallel = [float(m) for m in pt.mass(sequences)]
        
    
    # Benchmark Pyteomics
    pyteomics_times = []
    all_masses_pyto = []
    
    print(f"Benchmarking Pyteomics ({iterations} iterations per sequence)...")
    for seq in sequences:
        for _ in range(iterations):
            start = time.perf_counter()
            result_pyto = mass.calculate_mass(proforma=seq)
            end = time.perf_counter()
            pyteomics_times.append((end - start) * 1_000_000)  # microseconds
        
        all_masses_pyto.append(float(mass.calculate_mass(proforma=seq)))
    
    # Results
    pept_mean = statistics.mean(peptacular_times)
    pept_median = statistics.median(peptacular_times)
    pept_stdev = statistics.stdev(peptacular_times)
    
    pept_parallel_mean = statistics.mean(peptacular_parallel_times)
    pept_parallel_median = statistics.median(peptacular_parallel_times)
    pept_parallel_stdev = statistics.stdev(peptacular_parallel_times)
    
    pyto_mean = statistics.mean(pyteomics_times)
    pyto_median = statistics.median(pyteomics_times)
    pyto_stdev = statistics.stdev(pyteomics_times)
    
    speedup = pyto_mean / pept_mean
    speedup_parallel = (pyto_mean * num_sequences) / pept_parallel_mean  # Compare total time
    
    # Check mass agreement
    mass_differences = [abs(m1 - m2) for m1, m2 in zip(all_masses_pept, all_masses_pyto)]
    mass_differences_parallel = [abs(m1 - m2) for m1, m2 in zip(all_masses_pept_parallel, all_masses_pyto)]
    max_diff = max(mass_differences)
    avg_diff = statistics.mean(mass_differences)
    max_diff_parallel = max(mass_differences_parallel)
    avg_diff_parallel = statistics.mean(mass_differences_parallel)
    
    print(f"\n{'='*70}")
    print(f"Mass Calculation Benchmark - Random ProForma Sequences")
    print(f"Sequences: {num_sequences}")
    print(f"Iterations per sequence: {iterations:,}")
    print(f"Total calculations: {num_sequences * iterations:,}")
    print(f"{'='*70}")
    print(f"\n{'Method':<25} {'Mean (µs)':<15} {'Median (µs)':<15} {'StdDev (µs)':<15}")
    print(f"{'-'*70}")
    print(f"{'Peptacular (sequential)':<25} {pept_mean:>13.2f}   {pept_median:>13.2f}   {pept_stdev:>13.2f}")
    print(f"{'Peptacular (parallel)':<25} {pept_parallel_mean:>13.2f}   {pept_parallel_median:>13.2f}   {pept_parallel_stdev:>13.2f}")
    print(f"{'Pyteomics':<25} {pyto_mean:>13.2f}   {pyto_median:>13.2f}   {pyto_stdev:>13.2f}")
    print(f"{'-'*70}")
    print(f"\nMass Agreement (vs Pyteomics):")
    print(f"  Sequential - Max diff: {max_diff:.6f} Da, Avg diff: {avg_diff:.6f} Da")
    print(f"  Parallel   - Max diff: {max_diff_parallel:.6f} Da, Avg diff: {avg_diff_parallel:.6f} Da")
    
    print(f"\nPerformance Comparison:")
    if speedup > 1:
        print(f"  Sequential: Peptacular is {speedup:.2f}x FASTER than Pyteomics ✅")
    else:
        print(f"  Sequential: Pyteomics is {1/speedup:.2f}x FASTER than Peptacular ⚠️")
    
    if speedup_parallel > 1:
        print(f"  Parallel:   Peptacular is {speedup_parallel:.2f}x FASTER than Pyteomics (sequential) ✅")
    else:
        print(f"  Parallel:   Pyteomics is {1/speedup_parallel:.2f}x FASTER than Peptacular (parallel) ⚠️")
    
    # Compare sequential vs parallel within Peptacular
    total_seq_time = pept_mean * num_sequences
    parallel_speedup = total_seq_time / pept_parallel_mean
    print(f"\n  Peptacular parallel is {parallel_speedup:.2f}x faster than sequential 🚀")
    
    print(f"{'='*70}\n")


if __name__ == "__main__":
    benchmark_mass_calculation()