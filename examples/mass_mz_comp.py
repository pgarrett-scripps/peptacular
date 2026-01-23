"""
Mass and Composition Calculations
==================================
Examples of calculating mass, m/z, and elemental composition from ProForma annotations.
"""

import peptacular as pt


def run():
    # Parse a simple peptide sequence
    annot = pt.parse("PEPTIDE")

    # ============================================================================
    # MASS CALCULATIONS
    # ============================================================================

    print("=" * 60)
    print("MASS CALCULATIONS")
    print("=" * 60)

    # --- Basic Mass Calculation ---
    # Default is monoisotopic precursor mass (includes terminal groups H and OH)
    mass = annot.mass()
    print(f"Default mass: {mass:.4f} Da")  # 799.3600

    # Explicitly specify precursor ion type
    mass_p = annot.mass(ion_type="p")
    print(f"Precursor mass: {mass_p:.4f} Da")

    # --- Neutral Mass (no terminal groups) ---
    neutral = annot.mass(ion_type=pt.IonType.NEUTRAL)
    print(f"Neutral mass: {neutral:.4f} Da")

    # --- m/z Calculation ---
    # mz() divides mass by charge
    mz_2plus = annot.mz(charge=2)
    assert mz_2plus == annot.mass(charge=2) / 2
    print(f"m/z at charge +2: {mz_2plus:.4f}")

    # --- Monoisotopic vs Average Mass ---
    mono_mass = annot.mass(monoisotopic=True)  # default
    avg_mass = annot.mass(monoisotopic=False)
    print(f"Monoisotopic: {mono_mass:.4f} Da")
    print(f"Average: {avg_mass:.4f} Da")

    # --- Charge States ---
    # Integer charge assumes protonation/deprotonation
    mass_2plus = annot.mass(charge=2)
    mass_2minus = annot.mass(charge=-2)
    print(f"Mass at +2 charge: {mass_2plus:.4f} Da")
    print(f"Mass at -2 charge: {mass_2minus:.4f} Da")

    # Adduct charges (overrides annotation charge)
    mass_na = annot.mass(charge="Na:z+1")
    mass_multi_adduct = annot.mass(charge=("Na:z+1^2", "H:z+1"))
    print(f"Mass with Na+ adduct: {mass_na:.4f} Da")
    print(f"Mass with multiple adducts: {mass_multi_adduct:.4f} Da")

    # --- Isotopes ---
    # Integer assumes C13 isotopes
    mass_c13 = annot.mass(isotopes=1)
    print(f"Mass with 1x 13C: {mass_c13:.4f} Da")

    # Custom isotope specification
    mass_custom_iso = annot.mass(isotopes={"17O": 2, "13C": 1})
    print(f"Mass with 2x 17O and 1x 13C: {mass_custom_iso:.4f} Da")

    # --- Neutral Losses ---
    # Single loss
    mass_water_loss = annot.mass(deltas={"H2O": 1})
    print(f"Mass with H2O loss: {mass_water_loss:.4f} Da")

    # Multiple losses
    mass_multi_loss = annot.mass(
        deltas={pt.NeutralDelta.WATER: 1, pt.NeutralDelta.AMMONIA: 2}
    )
    print(f"Mass with H2O + 2×NH3 loss: {mass_multi_loss:.4f} Da")

    # ============================================================================
    # COMPOSITION CALCULATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("COMPOSITION CALCULATIONS")
    print("=" * 60)

    # --- Basic Composition ---
    # Returns a Counter of ElementInfo objects
    comp = annot.comp()
    print("\nFull composition (ElementInfo objects):")
    for elem, count in comp.items():
        print(f"  {elem.symbol}: {count}")

    # Convert to simple string representation
    comp_str = {str(elem): count for elem, count in comp.items()}
    print(f"\nSimple composition: {comp_str}")

    # --- Composition with Modifications ---
    # Apply charge and isotopes
    comp_modified = annot.comp(charge="Na:z+1", isotopes={"17O": 2, "13C": 1})
    comp_modified_str = {str(elem): count for elem, count in comp_modified.items()}
    print(f"\nModified composition: {comp_modified_str}")

    # ============================================================================
    # WITH GLOBAL ISOTOPE MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("GLOBAL ISOTOPE MODIFICATIONS")
    print("=" * 60)

    # Applies the global isotope to all residues
    iso_annot = pt.parse("<13C>PEPTIDE")
    iso_comp = iso_annot.comp()
    iso_comp_str = {str(elem): count for elem, count in iso_comp.items()}
    print(f"Composition with global 13C: {iso_comp_str}")

    # will also apply isotopes to modifications where available
    # Charge is applied after isotopes so will reflect in composition
    mod_iso_annot = pt.parse("<2H>PEPT[Phospho]IDE/2")
    mod_iso_comp = mod_iso_annot.comp()
    mod_iso_comp_str = {str(elem): count for elem, count in mod_iso_comp.items()}
    print(f"Composition with global 13C and Phospho mod: {mod_iso_comp_str}")


if __name__ == "__main__":
    run()
