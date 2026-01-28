"""
ProForma Mod Objects Examples
==============================
Demonstrates working with Mods and Mod objects returned from ProForma annotations.
"""

import peptacular as pt


def run():
    # ============================================================================
    # ACCESSING MODS OBJECTS
    # ============================================================================

    # Parse a ProForma annotation with various modifications
    annot = pt.parse("[Acetyl][Formyl]-PEM[Oxidation][Phospho]TIS[Phospho]DE/2")
    print(f"Annotation: {annot.serialize()}\n")

    # Access different mod collections - these return Mods objects
    print("N-terminal mods:", annot.nterm_mods)
    print("Internal mods at pos 2:", annot.get_internal_mods_at_index(2))
    print("Internal mods at pos 5:", annot.get_internal_mods_at_index(5))

    # ============================================================================
    # ITERATING OVER MODS
    # ============================================================================

    print("\n" + "=" * 60)
    print("ITERATING OVER MODS")
    print("=" * 60)

    annot = pt.parse("[Acetyl][Acetyl][Formyl]-PEPTIDE")
    nterm = annot.nterm_mods

    print(f"N-terminal mods: {nterm}\n")

    # Iterate through Mod objects
    print("Individual Mod objects:")
    for mod in nterm:
        print(f"  {mod.value} (count: {mod.count})")

    # ============================================================================
    # WORKING WITH MOD OBJECTS
    # ============================================================================

    print("\n" + "=" * 60)
    print("MOD OBJECT PROPERTIES")
    print("=" * 60)

    annot = pt.parse("PEM[Oxidation][Oxidation][Phospho]TIDE")
    internal_mods = annot.get_internal_mods_at_index(2)

    print(f"Internal mods at position 2: {internal_mods}\n")

    for mod in internal_mods:
        print(f"Modification: {mod.value}")
        print(f"  Count: {mod.count}")
        print(f"  Mass (mono): {mod.get_mass(monoisotopic=True):.4f}")
        print(f"  Mass (avg): {mod.get_mass(monoisotopic=False):.4f}")
        print(f"  Composition: {mod.get_composition()}")
        print(f"  Charge: {mod.get_charge()}")
        print()

    # ============================================================================
    # ACCESSING PARSED MOD VALUES
    # ============================================================================

    print("\n" + "=" * 60)
    print("PARSED MOD VALUES")
    print("=" * 60)

    annot = pt.parse("PEM[Oxidation]TIS[Phospho]DE")

    # Get mods at position 2 (M with Oxidation)
    mods_at_2 = annot.get_internal_mods_at_index(2)
    print(f"Mods at position 2: {mods_at_2}\n")

    # Access parsed items as (modification, count) tuples
    print("Parsed items:")
    for mod_value, count in mods_at_2.parse_items():
        print(f"  {mod_value} × {count}")
        print(f"    Type: {type(mod_value)}")

    # ============================================================================
    # CHECKING MOD PRESENCE
    # ============================================================================

    print("\n" + "=" * 60)
    print("CHECKING MOD PRESENCE")
    print("=" * 60)

    annot = pt.parse("[Acetyl][Formyl]-PEM[Oxidation]TIDE")
    nterm = annot.nterm_mods

    print(f"N-terminal mods: {nterm}\n")
    print(f"Contains 'Acetyl': {'Acetyl' in nterm}")
    print(f"Contains 'Phospho': {'Phospho' in nterm}")
    print(f"Contains 'Formyl': {'Formyl' in nterm}")

    # ============================================================================
    # WORKING WITH DIFFERENT MOD TYPES
    # ============================================================================

    print("\n" + "=" * 60)
    print("DIFFERENT MOD TYPES")
    print("=" * 60)

    # Isotope modifications
    annot = pt.parse("<15N>PEPTIDE")
    isotope = annot.isotope_mods
    print(f"Isotope mods: {isotope}")
    print(f"  Type: {isotope.mod_type}")
    print(f"  Serialized: {isotope.serialize()}\n")

    # Static modifications
    annot = pt.parse("<[Carbamidomethyl]@C>PEPTCDE")
    static = annot.static_mods
    print(f"Static mods: {static}")
    print(f"  Type: {static.mod_type}")
    print(f"  Serialized: {static.serialize()}\n")

    # Labile modifications
    annot = pt.parse("{Glycan:Hex}PEPTIDE")
    labile = annot.labile_mods
    print(f"Labile mods: {labile}")
    print(f"  Type: {labile.mod_type}")
    print(f"  Serialized: {labile.serialize()}\n")

    # Charge adducts
    annot = pt.parse("PEPTIDE/[Na:z+1]")
    charge = annot.charge_adducts
    print(f"Charge adducts: {charge}")
    print(f"  Type: {charge.mod_type}")
    print(f"  Serialized: {charge.serialize()}")

    # ============================================================================
    # MOD COPYING
    # ============================================================================

    print("\n" + "=" * 60)
    print("COPYING MOD OBJECTS")
    print("=" * 60)

    annot = pt.parse("[Acetyl]-PEPTIDE")
    nterm = annot.nterm_mods

    print(f"Original: {nterm}")

    # Copy a Mods collection
    nterm_copy = nterm.copy()
    print(f"Copy: {nterm_copy}")
    print(f"Are they equal: {nterm._mods == nterm_copy._mods}")
    print(f"Are they the same object: {nterm is nterm_copy}")

    # ============================================================================
    # WORKING WITH MULTIPLE MODS AT SAME POSITION
    # ============================================================================

    print("\n" + "=" * 60)
    print("MULTIPLE MODS AT SAME POSITION")
    print("=" * 60)

    annot = pt.parse("PEM[Oxidation][Oxidation][Phospho][Acetyl]TIDE")
    mods = annot.get_internal_mods_at_index(2)

    print(f"Mods at position 2: {mods}\n")

    print("Individual modifications:")
    for mod in mods:
        print(f"  {mod.value} × {mod.count}")

    print(f"\nTotal mass contribution: {mods.get_mass():.4f}")
    print(f"Total composition: {mods.get_composition()}")

    # ============================================================================
    # ACCESSING UNDERLYING PROFORMA COMPONENTS
    # ============================================================================

    print("\n" + "=" * 60)
    print("UNDERLYING PROFORMA COMPONENTS")
    print("=" * 60)

    # Different modification formats
    annot1 = pt.parse("PEM[Oxidation]TIDE")
    annot2 = pt.parse("PEM[UNIMOD:35]TIDE")
    annot3 = pt.parse("PEM[+15.995]TIDE")
    annot4 = pt.parse("PEM[Formula:O]TIDE")

    print("Different representations of Oxidation:\n")

    for i, annot in enumerate([annot1, annot2, annot3, annot4], 1):
        mods = annot.get_internal_mods_at_index(2)
        for mod in mods:
            print(f"{i}. {annot.serialize()}")
            print(f"   Parsed value: {mod.value}")
            print(f"   Type: {type(mod.value).__name__}")
            print(f"   Mass: {mod.get_mass():.4f}\n")

    # ============================================================================
    # VALIDATING MODS
    # ============================================================================

    print("=" * 60)
    print("VALIDATING MODS")
    print("=" * 60)

    annot = pt.parse("[Acetyl]-PEM[Oxidation]TIDE")
    nterm = annot.nterm_mods
    internal = annot.get_internal_mods_at_index(2)

    print(f"N-terminal mods valid: {nterm.is_valid}")
    print(f"Validation result: {nterm.validate()}")

    print(f"\nInternal mods valid: {internal.is_valid}")
    print(f"Validation result: {internal.validate()}")

    # ============================================================================
    # WORKING WITH INTERVALS
    # ============================================================================

    print("\n" + "=" * 60)
    print("INTERVAL MODIFICATIONS")
    print("=" * 60)

    annot = pt.parse("PEP(TIS)[Phospho]DE")

    print(f"Annotation: {annot.serialize()}\n")

    if annot.has_intervals:
        print("Intervals:")
        for interval in annot.intervals:
            print(f"  Range: {interval.start}-{interval.end}")
            print(f"  Ambiguous: {interval.ambiguous}")
            print(f"  Has mods: {interval.has_mods}")
            if interval.has_mods:
                print(f"  Mods: {interval.mods}")
                for mod in interval.mods:
                    print(f"    {mod.value} × {mod.count}")

    # ============================================================================
    # SERIALIZATION
    # ============================================================================

    print("\n" + "=" * 60)
    print("SERIALIZATION")
    print("=" * 60)

    annot = pt.parse("[Acetyl][Formyl]-PEM[Oxidation][Phospho]TIDE/2")

    print("Full annotation:", annot.serialize())
    print("\nIndividual mod serialization:")
    print(f"  N-term: {annot.nterm_mods.serialize()}")
    print(f"  Position 2: {annot.get_internal_mods_at_index(2).serialize()}")

    # Show how mods serialize differently based on type
    print("\nMod type serialization patterns:")
    examples = [
        ("[Acetyl]-PEPTIDE", "nterm_mods"),
        ("PEPTIDE-[Amidated]", "cterm_mods"),
        ("{Glycan:Hex}PEPTIDE", "labile_mods"),
        ("[Phospho]?PEPTIDE", "unknown_mods"),
        ("<15N>PEPTIDE", "isotope_mods"),
        ("PEPTIDE/[Na:z+1]", "charge_adducts"),
    ]

    for proforma, attr in examples:
        annot = pt.parse(proforma)
        mods = getattr(annot, attr)
        print(f"  {proforma:30s} -> {mods.serialize()}")


if __name__ == "__main__":
    run()
