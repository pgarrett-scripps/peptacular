"""
Benchmark script for ProForma parser performance.

Tests various ProForma sequence complexity levels:
- Simple unmodified peptides
- Modified peptides with various modification types
- Complex sequences with intervals, ambiguity, etc.
- Chimeric and crosslinked sequences
"""

import random
import time
from collections import defaultdict

import peptacular as pt


def generate_random_sequence(length: int = 10) -> str:
    """Generate a random amino acid sequence."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(amino_acids) for _ in range(length))


def generate_test_sequences(count: int = 1000) -> dict[str, list[str]]:
    """Generate test sequences of various complexity levels."""
    sequences = defaultdict(list)

    # 1. Simple unmodified sequences
    for _ in range(count):
        sequences["simple"].append(generate_random_sequence(random.randint(7, 30)))

    # 2. Single modification
    for _ in range(count):
        seq = generate_random_sequence(random.randint(7, 30))
        pos = random.randint(0, len(seq) - 1)
        sequences["single_mod"].append(seq[:pos] + seq[pos] + "[Oxidation]" + seq[pos + 1 :])

    # 3. Multiple modifications
    for _ in range(count):
        seq = generate_random_sequence(random.randint(10, 30))
        mod_count = random.randint(2, 5)
        positions = random.sample(range(len(seq)), min(mod_count, len(seq)))
        parts = []
        last = 0
        for pos in sorted(positions):
            parts.append(seq[last:pos])
            parts.append(seq[pos] + "[Oxidation]")
            last = pos + 1
        parts.append(seq[last:])
        sequences["multi_mod"].append("".join(parts))

    # 4. Terminal modifications
    for _ in range(count):
        seq = generate_random_sequence(random.randint(7, 30))
        sequences["term_mod"].append(f"[Acetyl]-{seq}-[Amidated]")

    # 5. Global modifications
    for _ in range(count):
        seq = generate_random_sequence(random.randint(7, 30))
        sequences["global_mod"].append(f"<[Oxidation]@M>{seq}")

    # 6. Labile modifications
    for _ in range(count):
        seq = generate_random_sequence(random.randint(7, 30))
        sequences["labile_mod"].append(f"{{Glycan:Hex}}{seq}")

    # 7. Intervals
    for _ in range(count):
        seq1 = generate_random_sequence(random.randint(3, 5))
        seq2 = generate_random_sequence(random.randint(4, 8))
        seq3 = generate_random_sequence(random.randint(3, 5))
        sequences["interval"].append(f"{seq1}({seq2}){seq3}")

    # 8. Ambiguous positions
    for _ in range(count):
        seq = generate_random_sequence(random.randint(7, 30))
        sequences["ambiguous"].append(f"[Phospho]?{seq}")

    # 9. Charge states
    for _ in range(count):
        seq = generate_random_sequence(random.randint(7, 30))
        charge = random.randint(1, 4)
        sequences["charge"].append(f"{seq}/{charge}")

    # 10. Complex (combination of features)
    for _ in range(count):
        seq = generate_random_sequence(random.randint(10, 20))
        sequences["complex"].append(f"<13C>[Acetyl]-{seq[:3]}[Oxidation]({seq[3:8]})[Phospho]{seq[8:]}-[Amidated]/2")

    return sequences


def benchmark_parsing(sequences: dict[str, list[str]], iterations: int = 1):
    """Benchmark parsing for each sequence type."""
    results = {}

    print("=" * 80)
    print("BENCHMARKING PROFORMA PARSER")
    print("=" * 80)
    print(f"Sequences per category: {len(next(iter(sequences.values()))):,}")
    print(f"Iterations: {iterations}")
    print("=" * 80)

    for seq_type, seq_list in sequences.items():
        print(f"\nTesting: {seq_type.upper().replace('_', ' ')}")
        print("-" * 40)

        # Warm-up
        for seq in seq_list[:10]:
            try:
                _ = pt.ProFormaAnnotation.parse(seq)
            except Exception:
                pass

        # Benchmark
        start_time = time.perf_counter()
        successful = 0
        failed = 0

        for _ in range(iterations):
            for seq in seq_list:
                try:
                    _ = pt.ProFormaAnnotation.parse(seq)
                    successful += 1
                except Exception:
                    failed += 1

        elapsed = time.perf_counter() - start_time

        total = successful + failed
        results[seq_type] = {
            "time": elapsed,
            "count": len(seq_list) * iterations,
            "successful": successful,
            "failed": failed,
            "us_per_seq": elapsed / total * 1_000_000,
            "sequences_per_sec": total / elapsed,
        }

        print(f"  Time: {elapsed:.3f} seconds")
        print(f"  Sequences parsed: {successful:,}/{total:,}")
        if failed > 0:
            print(f"  Failed: {failed:,}")
        print(f"  Speed: {elapsed / total * 1_000_000:.2f} μs per sequence")
        print(f"  Throughput: {total / elapsed:,.0f} sequences/sec")

    return results


def benchmark_parse_methods():
    """Compare parse vs from_proforma methods."""
    sequences = [generate_random_sequence(15) for _ in range(1000)]

    print("\n" + "=" * 80)
    print("COMPARING PARSING METHODS")
    print("=" * 80)

    methods = {
        "parse": lambda s: pt.ProFormaAnnotation.parse(s),
        "from_proforma": lambda s: pt.ProFormaAnnotation.from_proforma(s),
    }

    for method_name, method_func in methods.items():
        print(f"\nTesting: {method_name}")
        print("-" * 40)

        # Warm-up
        for seq in sequences[:10]:
            _ = method_func(seq)

        # Benchmark
        start_time = time.perf_counter()
        for seq in sequences:
            _ = method_func(seq)
        elapsed = time.perf_counter() - start_time

        print(f"  Time: {elapsed:.3f} seconds")
        print(f"  Speed: {elapsed / len(sequences) * 1_000_000:.2f} μs per sequence")
        print(f"  Throughput: {len(sequences) / elapsed:,.0f} sequences/sec")


def profile_parser_components():
    """Profile individual parser components."""
    import cProfile
    import pstats
    from io import StringIO

    sequences = generate_test_sequences(count=100)

    # Combine all sequences
    all_sequences = []
    for seq_list in sequences.values():
        all_sequences.extend(seq_list)

    print("\n" + "=" * 80)
    print("PROFILING PARSER (Top 30 functions)")
    print("=" * 80)

    profiler = cProfile.Profile()
    profiler.enable()

    for seq in all_sequences:
        try:
            _ = pt.ProFormaAnnotation.parse(seq)
        except Exception:
            pass

    profiler.disable()

    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
    ps.print_stats(30)
    print(s.getvalue())


def benchmark_sequence_lengths():
    """Benchmark parsing performance vs sequence length."""
    lengths = [5, 10, 20, 50, 100, 200, 500]

    print("\n" + "=" * 80)
    print("PARSING PERFORMANCE VS SEQUENCE LENGTH")
    print("=" * 80)

    print(f"{'Length':<10} {'Time (ms)':<12} {'μs/seq':<12} {'Seqs/sec':<15}")
    print("-" * 80)

    for length in lengths:
        sequences = [generate_random_sequence(length) for _ in range(100)]

        start_time = time.perf_counter()
        for seq in sequences:
            _ = pt.ProFormaAnnotation.parse(seq)
        elapsed = time.perf_counter() - start_time

        print(f"{length:<10} {elapsed * 1000:<12.3f} {elapsed / len(sequences) * 1_000_000:<12.2f} {len(sequences) / elapsed:<15,.0f}")


def main():
    """Run all benchmarks."""
    # Generate test sequences
    print("Generating test sequences...")
    sequences = generate_test_sequences(count=1000)

    # Run benchmarks
    results = benchmark_parsing(sequences, iterations=1)

    # Compare methods
    benchmark_parse_methods()

    # Benchmark vs length
    benchmark_sequence_lengths()

    # Profile
    profile_parser_components()

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    fastest = min(results, key=lambda x: results[x]["time"])
    slowest = max(results, key=lambda x: results[x]["time"])

    print(f"\nFastest: {fastest} ({results[fastest]['us_per_seq']:.2f} μs/seq)")
    print(f"Slowest: {slowest} ({results[slowest]['us_per_seq']:.2f} μs/seq)")

    avg_time = sum(r["time"] for r in results.values()) / len(results)
    avg_throughput = sum(r["sequences_per_sec"] for r in results.values()) / len(results)

    print(f"\nAverage throughput: {avg_throughput:,.0f} sequences/sec")

    print("=" * 80)


if __name__ == "__main__":
    main()
