"""
Profiling script for ProForma parser with detailed analysis.

Uses cProfile to identify hotspots in the parsing code.
"""

import cProfile
import pstats
import random
from io import StringIO

import peptacular as pt


def generate_random_sequence(length: int = 10) -> str:
    """Generate a random amino acid sequence."""
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(amino_acids) for _ in range(length))


def generate_complex_sequences(count: int = 100) -> list[str]:
    """Generate complex ProForma sequences for profiling."""
    sequences = []

    # Mix of different complexity levels
    for _ in range(count // 10):
        # Simple
        sequences.append(generate_random_sequence(20))

    for _ in range(count // 10):
        # Single mod
        seq = generate_random_sequence(20)
        pos = random.randint(0, len(seq) - 1)
        sequences.append(seq[:pos] + seq[pos] + "[Oxidation]" + seq[pos + 1 :])

    for _ in range(count // 10):
        # Multiple mods
        seq = generate_random_sequence(30)
        sequences.append(f"{seq[:5]}[Phospho]{seq[5:10]}[Oxidation]{seq[10:15]}[Acetyl]{seq[15:]}")

    for _ in range(count // 10):
        # Terminal mods
        seq = generate_random_sequence(20)
        sequences.append(f"[Acetyl]-{seq}-[Amidated]")

    for _ in range(count // 10):
        # Global mods
        seq = generate_random_sequence(20)
        sequences.append(f"<[Oxidation]@M>{seq}")

    for _ in range(count // 10):
        # Intervals
        seq = generate_random_sequence(20)
        sequences.append(f"{seq[:5]}({seq[5:12]}){seq[12:]}")

    for _ in range(count // 10):
        # Labile
        seq = generate_random_sequence(20)
        sequences.append(f"{{Glycan:Hex}}{seq}")

    for _ in range(count // 10):
        # Charge states
        seq = generate_random_sequence(20)
        sequences.append(f"{seq}/2")

    for _ in range(count // 10):
        # Ambiguous
        seq = generate_random_sequence(20)
        sequences.append(f"[Phospho]?{seq}")

    for _ in range(count // 10):
        # Complex combination
        seq = generate_random_sequence(25)
        sequences.append(f"<13C>[Acetyl]-{seq[:5]}[Oxidation]({seq[5:12]})[Phospho]{seq[12:]}-[Amidated]/2")

    return sequences


def profile_parsing(iterations: int = 10):
    """Profile the parsing function."""
    sequences = generate_complex_sequences(count=1000)

    print("=" * 80)
    print("PROFILING PROFORMA PARSER")
    print("=" * 80)
    print(f"Sequences: {len(sequences):,}")
    print(f"Iterations: {iterations}")
    print("=" * 80)

    profiler = cProfile.Profile()
    profiler.enable()

    for _ in range(iterations):
        for seq in sequences:
            try:
                annotation = pt.ProFormaAnnotation.parse(seq)
                # Access some properties to ensure full parsing
                _ = annotation.sequence
                _ = annotation.charge
            except Exception as e:
                print(f"Error parsing '{seq}': {e}")

    profiler.disable()

    # Print stats
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
    ps.print_stats(50)  # Top 50 functions
    print(s.getvalue())

    print("\n" + "=" * 80)
    print("BY TIME SPENT")
    print("=" * 80)
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("time")
    ps.print_stats(30)
    print(s.getvalue())

    print("\n" + "=" * 80)
    print("BY NUMBER OF CALLS")
    print("=" * 80)
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats("calls")
    ps.print_stats(30)
    print(s.getvalue())


def profile_specific_patterns():
    """Profile specific parsing patterns to identify bottlenecks."""
    patterns = {
        "simple": [generate_random_sequence(20) for _ in range(1000)],
        "single_mod": [generate_random_sequence(10) + "[Oxidation]" + generate_random_sequence(10) for _ in range(1000)],
        "multi_mod": [f"{generate_random_sequence(5)}[Oxidation]{generate_random_sequence(5)}[Phospho]{generate_random_sequence(5)}" for _ in range(1000)],
        "intervals": [f"{generate_random_sequence(5)}({generate_random_sequence(8)}){generate_random_sequence(5)}" for _ in range(1000)],
        "complex": [
            f"<13C>[Acetyl]-{generate_random_sequence(5)}[Oxidation]({generate_random_sequence(8)})[Phospho]{generate_random_sequence(5)}-[Amidated]/2"
            for _ in range(1000)
        ],
    }

    for pattern_name, sequences in patterns.items():
        print("\n" + "=" * 80)
        print(f"PROFILING: {pattern_name.upper()}")
        print("=" * 80)

        profiler = cProfile.Profile()
        profiler.enable()

        for seq in sequences:
            try:
                _ = pt.ProFormaAnnotation.parse(seq)
            except Exception:
                pass

        profiler.disable()

        s = StringIO()
        ps = pstats.Stats(profiler, stream=s).sort_stats("cumulative")
        ps.print_stats(20)
        print(s.getvalue())


if __name__ == "__main__":
    # Main profiling
    profile_parsing(iterations=10)

    # Pattern-specific profiling
    # profile_specific_patterns()
