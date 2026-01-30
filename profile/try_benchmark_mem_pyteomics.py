"""
Memory Footprint: Peptacular vs Pyteomics Annotation Objects
"""
import gc
import sys
import tracemalloc

from pympler import asizeof
from pyteomics import proforma

import peptacular as pt


def measure_memory_scaling():
    """Measure memory footprint at scale to see true per-object cost"""

    test_sequences = [
        ("Simple", "PEPTIDE"),
        ("Single Mod", "PEPTIDEM[+15.995]"),
        ("Multi Mod", "M[+15.995]C[+57.021]PEPTIDE"),
        ("Complex", "{Oxidation}[+15.995]-PEPTIDEM[+15.995]C[+57.021]KPEM[+15.995]PEPTIDEM[+15.995]-[+15.995]"),
        ("Long", "PEPTIDE" * 20 + "M[+15.995]" + "PEPTIDE" * 20),
    ]

    print(f"\n{'='*90}")
    print("Memory Footprint: ProForma Annotation Objects (Scaling)")
    print(f"{'='*90}\n")

    for name, seq in test_sequences:
        print(f"\n{name}: {seq[:60]}{'...' if len(seq) > 60 else ''}")
        print(f"{'-'*90}")

        # Test with multiple counts to see scaling
        for count in [1, 100, 1000]:
            gc.collect()  # Clean up before measurement

            # Peptacular measurement with tracemalloc
            tracemalloc.start()
            pept_objects = [pt.ProFormaAnnotation.parse(seq) for _ in range(count)]
            pept_current, pept_peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            # Pyteomics measurement with tracemalloc
            gc.collect()
            tracemalloc.start()
            pyto_objects = [proforma.ProForma.parse(seq) for _ in range(count)]
            pyto_current, pyto_peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            # Calculate per-object cost
            pept_per_obj = pept_current / count
            pyto_per_obj = pyto_current / count

            print(f"\n{count:>4} objects:")
            print(f"  Peptacular: {pept_current:>10,} bytes total ({pept_per_obj:>8,.1f} bytes/object)")
            print(f"  Pyteomics:  {pyto_current:>10,} bytes total ({pyto_per_obj:>8,.1f} bytes/object)")

            ratio = pyto_per_obj / pept_per_obj if pept_per_obj > 0 else float('inf')
            if ratio > 1:
                print(f"  ✅ Peptacular {ratio:.2f}x smaller ({pyto_per_obj - pept_per_obj:,.1f} bytes/object saved)")
            else:
                print(f"  ⚠️  Pyteomics {1/ratio:.2f}x smaller ({pept_per_obj - pyto_per_obj:,.1f} bytes/object saved)")

            # Clean up
            del pept_objects
            del pyto_objects
            gc.collect()

    print(f"\n{'='*90}\n")


def measure_deep_size_single():
    """Original approach with pympler - useful for understanding object structure"""

    test_sequences = [
        ("Simple", "PEPTIDE"),
        ("Multi Mod", "M[+15.995]C[+57.021]PEPTIDE"),
        ("Complex", "{Oxidation}[+15.995]-PEPTIDEM[+15.995]C[+57.021]KPEM[+15.995]PEPTIDEM[+15.995]-[+15.995]"),
    ]

    print(f"\n{'='*90}")
    print("Deep Object Size Analysis (includes all referenced objects)")
    print(f"{'='*90}\n")

    for name, seq in test_sequences:
        print(f"\n{name}: {seq[:60]}{'...' if len(seq) > 60 else ''}")
        print(f"{'-'*90}")

        pept_annot = pt.ProFormaAnnotation.parse(seq)
        pyto_annot = proforma.ProForma.parse(seq)

        pept_deep = asizeof.asizeof(pept_annot)
        pyto_deep = asizeof.asizeof(pyto_annot)

        # Shallow size for comparison
        pept_shallow = sys.getsizeof(pept_annot)
        pyto_shallow = sys.getsizeof(pyto_annot)

        print(f"\n{'Package':<15} {'Deep Size':<15} {'Shallow Size':<15} {'Overhead':<15}")
        print(f"{'-'*90}")
        print(f"{'Peptacular':<15} {pept_deep:>13,}  {pept_shallow:>13,}  {pept_deep - pept_shallow:>13,}")
        print(f"{'Pyteomics':<15} {pyto_deep:>13,}  {pyto_shallow:>13,}  {pyto_deep - pyto_shallow:>13,}")

        ratio = pyto_deep / pept_deep if pept_deep > 0 else float('inf')
        if ratio > 1:
            print(f"\n✅ Peptacular {ratio:.2f}x smaller deep size ({pyto_deep - pept_deep:,} bytes)")
        else:
            print(f"\n⚠️  Pyteomics {1/ratio:.2f}x smaller deep size ({pept_deep - pyto_deep:,} bytes)")


if __name__ == "__main__":
    # Scaling test shows real per-object cost
    measure_memory_scaling()
    
    # Deep size shows object structure overhead
    measure_deep_size_single()