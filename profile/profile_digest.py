import cProfile
import pstats
import random
from io import StringIO
import peptacular as pt
import time

def generate_random_protein(length: int = 500) -> str:
    """Generate a random protein sequence."""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(random.choice(amino_acids) for _ in range(length))

def run_digest():
    """Run digest on random proteins."""
    # Generate 500 DIFFERENT proteins instead of 500 copies of the same one
    proteins = [generate_random_protein(1_000) for _ in range(100)]
    
    start_time = time.time()
    peptides = list(pt.digest(
        proteins,
        cleave_on="KR",
        missed_cleavages=2,
        semi=False,
        min_len=6,
        max_len=50
    ))
    elapsed = time.time() - start_time
    
    total_peptides = len(peptides)
    print(f"Generated {total_peptides} peptides from {len(proteins)} proteins")
    print(f"Took {elapsed:.2f} seconds ({elapsed/len(proteins)*1000:.2f} ms per protein)")
    return peptides

def run_digest_single():
    """Run digest on a single protein to profile without multiprocessing."""
    protein = generate_random_protein(1_000)
    
    start_time = time.time()
    peptides = pt.digest(
        protein,  # Single protein
        cleave_on="KR",
        missed_cleavages=2,
        semi=False,
        min_len=6,
        max_len=50
    )
    elapsed = time.time() - start_time

    total_peptides = len(list(peptides))
    print(f"Generated {total_peptides} peptides from 1 protein")
    print(f"Took {elapsed*1000:.2f} ms")
    return peptides

def profile_digest_single():
    """Profile single digest (no multiprocessing overhead)."""
    profiler = cProfile.Profile()
    profiler.enable()
    result = run_digest_single()
    profiler.disable()
    
    # Print stats
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
    ps.print_stats(30)  # Top 30 functions
    print(s.getvalue())
    return result

def profile_digest_multi():
    """Profile multi-protein digest (with multiprocessing)."""
    profiler = cProfile.Profile()
    profiler.enable()
    result = run_digest()
    profiler.disable()
    
    # Print stats
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
    ps.print_stats(30)  # Top 30 functions
    print(s.getvalue())
    return result

if __name__ == "__main__":
    print("=" * 60)
    print("SINGLE PROTEIN DIGEST (no multiprocessing)")
    print("=" * 60)
    profile_digest_single()
    
    print("\n" + "=" * 60)
    print("MULTI-PROTEIN DIGEST (with multiprocessing)")
    print("=" * 60)
    #profile_digest_multi()