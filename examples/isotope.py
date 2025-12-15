"""
Isotopic Distribution Calculations
===================================
Examples of calculating isotopic distributions from ProForma annotations.
"""

import peptacular as pt

def run():

    # Parse a simple peptide sequence
    annot = pt.parse('PEPTIDE')

    # ============================================================================
    # BASIC ISOTOPIC DISTRIBUTION
    # ============================================================================

    print("=" * 60)
    print("BASIC ISOTOPIC DISTRIBUTION")
    print("=" * 60)

    # --- Default Distribution ---
    # Returns list of IsotopicData with mass, neutron_count, and abundance
    # Abundances normalized so max peak = 1.0
    dist = annot.isotopic_distribution()
    print(f"\nPeptide: {annot.serialize()}")
    print(f"Monoisotopic mass: {annot.mass():.3f} Da")
    print(f"Default isotopic distribution:")
    for iso in dist:
        print(f"  mass: {iso.mass:>8.3f} Da, abundance: {iso.abundance:>6.3f}, neutrons: {iso.neutron_count}")

    # --- Control Number of Isotopes ---
    dist_limited = annot.isotopic_distribution(max_isotopes=3)
    print(f"\nLimited to 3 most abundant isotopes:")
    for iso in dist_limited:
        print(f"  mass: {iso.mass:>8.3f} Da, abundance: {iso.abundance:>6.3f}")

    # --- Abundance Threshold ---
    # Only keep isotopes with abundance >= threshold (relative to max peak)
    dist_filtered = annot.isotopic_distribution(min_abundance_threshold=0.05)
    print(f"\nFiltered (≥5% of max peak):")
    for iso in dist_filtered:
        print(f"  mass: {iso.mass:>8.3f} Da, abundance: {iso.abundance:>6.3f}")

    # --- Neutron Offset Mode ---
    # Use neutron count instead of absolute mass (useful for matching patterns)
    dist_neutron = annot.isotopic_distribution(use_neutron_count=True)
    print(f"\nNeutron offset mode:")
    for iso in dist_neutron:
        print(f"  neutron offset: {iso.mass:>3.0f}, abundance: {iso.abundance:>6.3f}")

    # ============================================================================
    # DISTRIBUTION RESOLUTION
    # ============================================================================

    print("\n" + "=" * 60)
    print("DISTRIBUTION RESOLUTION")
    print("=" * 60)

    # --- High Resolution ---
    # More decimal places for precise mass calculations
    dist_high_res = annot.isotopic_distribution(distribution_resolution=5)
    print(f"\nHigh resolution (5 decimals):")
    for iso in dist_high_res[:3]:
        print(f"  mass: {iso.mass:.5f} Da, abundance: {iso.abundance:>6.3f}")

    # --- Low Resolution ---
    # Simulates lower instrument precision, combines nearby masses
    dist_low_res = annot.isotopic_distribution(distribution_resolution=2)
    print(f"\nLow resolution (2 decimals):")
    for iso in dist_low_res[:3]:
        print(f"  mass: {iso.mass:.2f} Da, abundance: {iso.abundance:>6.3f}")

    # ============================================================================
    # COMBINING WITH COMP PARAMETERS
    # ============================================================================

    print("\n" + "=" * 60)
    print("COMBINING WITH COMP PARAMETERS")
    print("=" * 60)

    # isotopic_distribution() accepts same parameters as comp()
    # Combine charge, isotopes, losses, and ion type
    dist_combined = annot.isotopic_distribution(
        ion_type=pt.IonType.Y,
        charge=2,
        isotopes=1,
        losses={pt.NeutralDelta.WATER: 1}
    )
    print(f"\ny-ion, +2 charge, +1 13C, -H2O:")
    for iso in dist_combined[:4]:
        print(f"  m/z: {iso.mass:>8.3f}, abundance: {iso.abundance:>6.3f}")


if __name__ == "__main__":
    run()
