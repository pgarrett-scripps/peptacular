"""
Memory Footprint: Peptacular vs Pyteomics Annotation Objects
"""
import sys

from pympler import asizeof
from pyteomics import proforma

import peptacular as pt


def profile_annotation_objects():
    """Compare memory footprint of parsed annotation objects"""
    
    # Test sequences with varying complexity
    test_sequences = [
        ("Simple", "PEPTIDE"),
        ("Single Mod", "PEPTIDEM[+15.995]"),
        ("Multi Mod", "M[+15.995]C[+57.021]PEPTIDE"),
        ("Complex", "{Oxidation}[+15.995]-PEPTIDEM[+15.995]C[+57.021]KPEM[+15.995]PEPTIDEM[+15.995]-[+15.995]"),
        ("Long", "PEPTIDE" * 20 + "M[+15.995]" + "PEPTIDE" * 20),
    ]
    
    print(f"\n{'='*90}")
    print(f"Memory Footprint: ProForma Annotation Objects")
    print(f"{'='*90}\n")
    
    for name, seq in test_sequences:
        print(f"\n{name}: {seq[:60]}{'...' if len(seq) > 60 else ''}")
        print(f"{'-'*90}")
        
        # Parse annotations
        pept_annot = pt.ProFormaAnnotation.parse(seq)
        pyto_annot = proforma.ProForma.parse(seq)
        
        # Get memory size using pympler (more accurate than sys.getsizeof)
        pept_size = asizeof.asizeof(pept_annot)
        pyto_size = asizeof.asizeof(pyto_annot)
        
        # Also show sys.getsizeof for comparison
        pept_size_sys = sys.getsizeof(pept_annot)
        pyto_size_sys = sys.getsizeof(pyto_annot)
        
        print(f"\n{'Package':<15} {'Deep Size (bytes)':<20} {'Shallow Size (bytes)':<25}")
        print(f"{'-'*90}")
        print(f"{'Peptacular':<15} {pept_size:>18,}    {pept_size_sys:>23,}")
        print(f"{'Pyteomics':<15} {pyto_size:>18,}    {pyto_size_sys:>23,}")
        print(f"{'-'*90}")
        
        # Memory efficiency
        ratio = pyto_size / pept_size if pept_size > 0 else float('inf')
        
        if ratio > 1:
            print(f"✅ Peptacular object is {ratio:.2f}x SMALLER ({pyto_size - pept_size:,} bytes saved)")
        else:
            print(f"⚠️  Pyteomics object is {1/ratio:.2f}x SMALLER ({pept_size - pyto_size:,} bytes saved)")
    
    print(f"\n{'='*90}\n")


if __name__ == "__main__":
    profile_annotation_objects()