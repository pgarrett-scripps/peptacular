import random
import time

import peptacular as pt


def generate_random_protein(length: int = 500) -> str:
    """Generate a random protein sequence."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(amino_acids) for _ in range(length))


def benchmark_digest_methods():
    """Benchmark all three parallelization methods."""

    # Generate test data
    proteins = [generate_random_protein(1_000) for _ in range(1000)]

    methods = ["sequential", "thread", "process"]
    results = {}

    print("=" * 80)
    print("BENCHMARKING DIGEST METHODS")
    print("=" * 80)
    print(f"Test data: {len(proteins)} proteins, ~1000 AA each")
    print("Parameters: cleave_on='KR', missed_cleavages=2, min_len=6, max_len=50")
    print("=" * 80)

    for method in methods:
        print(f"\nTesting method: {method.upper()}")
        print("-" * 40)

        # Warm-up run (not timed)
        _ = list(
            pt.digest(
                proteins[:10],
                cleave_on="KR",
                missed_cleavages=2,
                semi=False,
                min_len=6,
                max_len=50,
                method=method,
            )
        )

        # Actual benchmark
        start_time = time.time()
        peptides = list(
            pt.digest(
                proteins,
                cleave_on="KR",
                missed_cleavages=2,
                semi=False,
                min_len=6,
                max_len=50,
                method=method,
            )
        )
        elapsed = time.time() - start_time

        total_peptides = len(peptides)
        results[method] = {
            "time": elapsed,
            "peptides": total_peptides,
            "ms_per_protein": elapsed / len(proteins) * 1000,
        }

        print(f"  Time: {elapsed:.3f} seconds")
        print(f"  Peptides generated: {total_peptides:,}")
        print(f"  Speed: {elapsed / len(proteins) * 1000:.2f} ms per protein")
        print(f"  Throughput: {len(proteins) / elapsed:.1f} proteins/sec")

    # Summary comparison
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    fastest_method = min(results, key=lambda x: results[x]["time"])
    fastest_time = results[fastest_method]["time"]

    print(f"{'Method':<12} {'Time (s)':<12} {'Speedup':<12} {'ms/protein':<15}")
    print("-" * 80)

    for method in methods:
        time_val = results[method]["time"]
        speedup = results["sequential"]["time"] / time_val if method != "sequential" else 1.0
        ms_per_protein = results[method]["ms_per_protein"]

        marker = " ← FASTEST" if method == fastest_method else ""
        print(f"{method:<12} {time_val:<12.3f} {speedup:<12.2f}x {ms_per_protein:<15.2f}{marker}")

    print("=" * 80)

    return results


if __name__ == "__main__":
    benchmark_digest_methods()
