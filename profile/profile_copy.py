import time
import statistics
from typing import Callable

# Assuming your imports work like this - adjust as needed
from peptacular.proforma.dclasses.modlist import ModList
from peptacular.mod import Mod


def benchmark(func: Callable, iterations: int = 10000, warmup: int = 1000) -> dict:
    """Benchmark a function with warmup and multiple iterations."""
    # Warmup
    for _ in range(warmup):
        func()
    
    # Actual benchmark
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        func()
        end = time.perf_counter()
        times.append(end - start)
    
    return {
        'mean': statistics.mean(times),
        'median': statistics.median(times),
        'stdev': statistics.stdev(times) if len(times) > 1 else 0,
        'min': min(times),
        'max': max(times),
        'total': sum(times)
    }


def format_time(seconds: float) -> str:
    """Format time in appropriate units."""
    if seconds < 1e-6:
        return f"{seconds * 1e9:.2f} ns"
    elif seconds < 1e-3:
        return f"{seconds * 1e6:.2f} µs"
    elif seconds < 1:
        return f"{seconds * 1e3:.2f} ms"
    else:
        return f"{seconds:.2f} s"


def print_results(name: str, results: dict):
    """Print benchmark results nicely."""
    print(f"\n{name}:")
    print(f"  Mean:   {format_time(results['mean'])}")
    print(f"  Median: {format_time(results['median'])}")
    print(f"  Stdev:  {format_time(results['stdev'])}")
    print(f"  Min:    {format_time(results['min'])}")
    print(f"  Max:    {format_time(results['max'])}")


print("=" * 70)
print("ModList Initialization Benchmark")
print("=" * 70)

# Test 1: Create new ModList from scratch
def test_init_empty():
    return ModList(allow_dups=True, stackable=False)

results_init = benchmark(test_init_empty)
print_results("1. Init empty ModList()", results_init)


# Test 2: Copy a pre-created empty ModList
empty_modlist = ModList(allow_dups=True, stackable=False)

def test_copy_empty():
    return empty_modlist.copy()

results_copy = benchmark(test_copy_empty)
print_results("2. Copy empty ModList", results_copy)


# Test 3: Create ModList with object.__new__ (bypassing __init__)
def test_new_empty():
    result = object.__new__(ModList)
    result.allow_dups = True
    result.stackable = False
    result.data = []
    return result

results_new = benchmark(test_new_empty)
print_results("3. object.__new__() + manual setup", results_new)


# Test 4: Copy with optimized implementation (if you've implemented it)
def test_optimized_copy():
    result = object.__new__(ModList)
    result.allow_dups = empty_modlist.allow_dups
    result.stackable = empty_modlist.stackable
    result.data = empty_modlist.data.copy()
    return result

results_optimized = benchmark(test_optimized_copy)
print_results("4. Optimized copy (bypass __init__)", results_optimized)


# Test 5: Create ModList with one item
modlist_with_item = ModList(allow_dups=True, stackable=False)
modlist_with_item.append(Mod(15.99, 1))

def test_copy_with_item():
    return modlist_with_item.copy()

results_copy_item = benchmark(test_copy_with_item)
print_results("5. Copy ModList with 1 item", results_copy_item)


# Test 6: Create ModList with multiple items
modlist_with_items = ModList(allow_dups=True, stackable=False)
for i in range(5):
    modlist_with_items.append(Mod(15.99 + i, 1))

def test_copy_with_items():
    return modlist_with_items.copy()

results_copy_items = benchmark(test_copy_with_items)
print_results("6. Copy ModList with 5 items", results_copy_items)


# Comparison summary
print("\n" + "=" * 70)
print("COMPARISON SUMMARY")
print("=" * 70)

baseline = results_init['mean']
comparisons = [
    ("Init empty ModList()", results_init['mean'], 1.0),
    ("Copy empty ModList", results_copy['mean'], results_copy['mean'] / baseline),
    ("object.__new__() + setup", results_new['mean'], results_new['mean'] / baseline),
    ("Optimized copy", results_optimized['mean'], results_optimized['mean'] / baseline),
]

for name, time_val, ratio in comparisons:
    speedup = baseline / time_val
    print(f"\n{name}:")
    print(f"  Time:    {format_time(time_val)}")
    print(f"  Ratio:   {ratio:.2f}x baseline")
    print(f"  Speedup: {speedup:.2f}x faster" if speedup > 1 else f"  Slowdown: {1/speedup:.2f}x slower")

# Recommendation
print("\n" + "=" * 70)
print("RECOMMENDATION")
print("=" * 70)

fastest_name = min(comparisons, key=lambda x: x[1])[0]
fastest_time = min(comparisons, key=lambda x: x[1])[1]

print(f"\nFastest method: {fastest_name}")
print(f"Time: {format_time(fastest_time)}")

# Calculate potential savings in your use case
print("\n" + "=" * 70)
print("PROJECTED IMPACT ON YOUR WORKLOAD")
print("=" * 70)

# Assume 821 variants per peptide, 2.2M peptides total
num_copies_per_peptide = 821
num_peptides = 1000  # Sample size for projection

total_init = baseline * num_copies_per_peptide * num_peptides
total_optimized = fastest_time * num_copies_per_peptide * num_peptides
savings = total_init - total_optimized

print(f"\nFor {num_copies_per_peptide} variants × {num_peptides} peptides:")
print(f"  Baseline (init):    {format_time(total_init)}")
print(f"  Optimized:          {format_time(total_optimized)}")
print(f"  Time saved:         {format_time(savings)}")
print(f"  Improvement:        {((total_init - total_optimized) / total_init * 100):.1f}%")

print("\n" + "=" * 70)