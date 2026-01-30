"""
Simple parallel scaling test
"""

import multiprocessing as mp
import random
import time

import peptacular as pt

if __name__ == "__main__":
    pt.set_start_method("fork")

    def generate_random_proforma(length: int = 20) -> str:
        """Generate a random ProForma sequence."""
        return pt.ProFormaAnnotation.random().serialize()

    # Generate test data
    NUM_SEQUENCES = 10_000
    sequences = [generate_random_proforma(random.randint(10, 30)) for _ in range(NUM_SEQUENCES)]

    print(f"Testing {NUM_SEQUENCES} sequences")
    print(f"System CPUs: {mp.cpu_count()}")

    # Sequential baseline
    start = time.perf_counter()
    _ = [pt.mass(seq) for seq in sequences]
    baseline_time = time.perf_counter() - start

    print(f"Sequential: {baseline_time:.3f}s\n")

    # Test different worker counts
    worker_counts = [1, 2, 4, 8, 16, 20, 24, 32]
    max_workers = mp.cpu_count()
    worker_counts = [w for w in worker_counts if w <= max_workers * 2]

    print(f"{'Workers':<10} {'Sequential Time':<15} {'Speedup':<12} {'Process Time':<15} {'Speedup':<12}")
    print("-" * 70)

    for n in worker_counts:
        # Thread
        start = time.perf_counter()
        _ = pt.mass(sequences, n_workers=n, method="sequential")
        thread_time = time.perf_counter() - start
        thread_speedup = baseline_time / thread_time

        # Process
        start = time.perf_counter()
        _ = pt.mass(sequences, n_workers=n, method="process")
        process_time = time.perf_counter() - start
        process_speedup = baseline_time / process_time

        print(f"{n:<10} {thread_time:>10.3f}s     {thread_speedup:>8.2f}x     {process_time:>10.3f}s     {process_speedup:>8.2f}x")
