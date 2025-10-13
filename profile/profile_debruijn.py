import cProfile
import pstats
from io import StringIO
import peptacular as pt

# Test sequence - a longer protein for better profiling
test_sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
test_sequence = test_sequence * 5  # Make it longer
print(f"Sequence length: {len(test_sequence)}")
print(f"Profiling debruijn_sequence with k=3, seed=20...\n")

# Profile the function
profiler = cProfile.Profile()
profiler.enable()

result = pt.debruijin_sequence(test_sequence, k=3, static_residues="P", seed=20)

profiler.disable()

# Print stats
s = StringIO()
ps = pstats.Stats(profiler, stream=s).sort_stats('cumulative')
ps.print_stats(30)  # Top 30 functions
print(s.getvalue())

print(f"\nOriginal: {test_sequence[:50]}...")
print(f"Result:   {result[:50]}...")