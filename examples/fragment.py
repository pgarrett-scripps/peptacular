"""
Fragment Generation Examples
=============================
Examples of generating fragment ions from ProForma annotations.
All fragment methods return Fragment objects with mass, m/z, and composition.
"""

import peptacular as pt


def run():
    # ============================================================================
    # BASIC FRAGMENTATION
    # ============================================================================

    peptide = pt.parse("PEPT[Phospho]IDE-[Acetyl]")

    print("=" * 60)
    print("BASIC FRAGMENTATION")
    print("=" * 60)
    print(f"Peptide: {peptide}\n")

    # --- b-ions (N-terminal fragments) ---
    print("b-ions (N-terminal):")
    for frag in peptide.fragment(ion_types=["b"]):
        print(f"  {frag}")

    # --- y-ions (C-terminal fragments) ---
    print("\ny-ions (C-terminal):")
    for frag in peptide.fragment(ion_types=["y"]):
        print(f"  {frag}")

    # ============================================================================
    # FRAGMENT ION TYPES
    # ============================================================================

    print("\n" + "=" * 60)
    print("DIFFERENT ION TYPES")
    print("=" * 60)

    # Generate multiple ion types at once
    print("\na, b, c ions:")
    for frag in peptide.fragment(ion_types=["a", "b", "c"]):
        print(f"  {frag}")

    print("\nx, y, z ions:")
    for frag in peptide.fragment(ion_types=["x", "y", "z"]):
        print(f"  {frag}")

    # ============================================================================
    # CHARGED FRAGMENTS
    # ============================================================================

    print("\n" + "=" * 60)
    print("CHARGED FRAGMENTS")
    print("=" * 60)

    # Charge state
    print("\nb-ions at +2 charge:")
    for frag in peptide.fragment(ion_types=["b"], charges=[2]):
        print(f"  {frag}")

    # Adduct charges
    print("\ny-ions with Na+ adduct:")
    for frag in peptide.fragment(ion_types=["y"], charges=["Na:z+1"]):
        print(f"  {frag}")

    # ============================================================================
    # NEUTRAL LOSSES
    # ============================================================================

    print("\n" + "=" * 60)
    print("NEUTRAL LOSSES")
    print("=" * 60)

    # Water loss (Selectively applied to fragments that can lose H2O (containing ["S", "T", "D", "E"]))
    print("\ny-ions with H2O loss:")
    for frag in peptide.fragment(ion_types=["y"], deltas=[pt.NeutralDelta.WATER]):
        print(f"  {frag}")

    # Multiple losses
    print("\nb-ions with H2O and NH3 loss:")
    for frag in peptide.fragment(
        ion_types=["b"], deltas=[pt.NeutralDelta.WATER, pt.NeutralDelta.AMMONIA]
    ):
        print(f"  {frag}")

    # ============================================================================
    # ISOTOPES
    # ============================================================================

    print("\n" + "=" * 60)
    print("ISOTOPIC FRAGMENTS")
    print("=" * 60)

    # C13 isotopes
    print("\ny-ions with 1x 13C:")
    for frag in peptide.fragment(ion_types=["y"], isotope=[1]):
        print(f"  {frag}")

    # Custom isotopes
    print("\nb-ions with 2x 17O:")
    for frag in peptide.fragment(ion_types=["b"], isotope=[{"17O": 2}]):
        print(f"  {frag}")

    # ============================================================================
    # MODIFIED PEPTIDES
    # ============================================================================

    print("\n" + "=" * 60)
    print("MODIFIED PEPTIDES")
    print("=" * 60)

    modified = pt.parse("[Acetyl]-PEM[Oxidation]TIDES[Phospho]")
    print(f"Modified peptide: {modified}\n")

    print("y-ions (modifications preserved in fragments):")
    for frag in modified.fragment(ion_types=["y"]):
        print(f"  {frag}")

    # ============================================================================
    # INTERNAL FRAGMENTS
    # ============================================================================

    print("\n" + "=" * 60)
    print("INTERNAL FRAGMENTS")
    print("=" * 60)

    print("\nInternal fragments (min_len=3, max_len=5):")
    for frag in peptide.fragment(ion_types=["ax"]):
        if frag.position and isinstance(frag.position, tuple):
            start, end = frag.position
            if 3 <= (end - start) <= 5:
                print(f"  {frag}")

    # ============================================================================
    # IMMONIUM IONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("IMMONIUM IONS")
    print("=" * 60)

    print("\nImmonium ions:")
    for frag in peptide.fragment(ion_types=["i"]):
        print(f"  {frag}")

    # ============================================================================
    # PRECURSOR ION
    # ============================================================================

    print("\n" + "=" * 60)
    print("PRECURSOR ION")
    print("=" * 60)

    print("\nPrecursor ion at +2 charge:")
    for frag in peptide.fragment(ion_types=["p"], charges=[2]):
        print(f"  {frag}")

    # ============================================================================
    # COMBINING OPTIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("COMBINING OPTIONS")
    print("=" * 60)

    print("\ny-ions: +2 charge, H2O loss, 1x 13C:")
    for frag in peptide.fragment(
        ion_types=["y"], charges=[2], deltas=[pt.NeutralDelta.WATER], isotope=[1]
    ):
        print(f"  {frag}")

    # ============================================================================
    # ACCESSING FRAGMENT PROPERTIES
    # ============================================================================

    print("\n" + "=" * 60)
    print("FRAGMENT PROPERTIES")
    print("=" * 60)

    """
    Unless otherwise specified, fragments do not include sequence or composition data.
    This can be enabled with the `include_sequence` and `calculate_composition` flags.
    """
    b_ions = peptide.fragment(
        ion_types=["b"], charges=[2], include_sequence=True, calculate_composition=True
    )
    if b_ions:
        frag = b_ions[0]
        print(f"\nExample fragment: {frag}")
        print(f"  Ion type: {frag.ion_type}")
        print(f"  Position: {frag.position}")
        print(f"  Mass: {frag.mass:.4f} Da")
        print(f"  m/z: {frag.mz:.4f}")
        print(f"  Charge: {frag.charge_state}")
        print(f"  Neutral mass: {frag.neutral_mass:.4f} Da")
        if frag.composition:
            comp_str = {str(elem): count for elem, count in frag.composition.items()}
            print(f"  Composition: {comp_str}")
        if frag.sequence:
            print(f"  Sequence: {frag.sequence}")

    # ============================================================================
    # MZPAF OUTPUT
    # ============================================================================

    print("\n" + "=" * 60)
    print("MZPAF OUTPUT")
    print("=" * 60)

    print("\nFragment annotations in mzPAF format:")
    fragments = peptide.fragment(ion_types=["b", "y"], charges=[2])
    for frag in fragments[:8]:  # Show first 8 fragments
        mzpaf = frag.to_mzpaf()
        print(f"  {mzpaf}")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    run()
