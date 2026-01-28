"""
ProForma Annotation Examples
=============================
Basic examples of parsing, serializing, and manipulating ProForma annotations.
"""

import peptacular as pt


def run():
    # ============================================================================
    # PARSING ANNOTATIONS
    # ============================================================================

    # Simple sequence
    simple: pt.ProFormaAnnotation = pt.parse("PEPTIDE")
    print(f"Simple: {simple.serialize()}")

    # Chimeric sequence
    chimeric: list[pt.ProFormaAnnotation] = pt.parse_chimeric("PEPTIDE+PEPTIDE")
    print(f"Chimeric: {pt.serialize_chimeric(chimeric)}")

    # ============================================================================
    # CREATING ANNOTATIONS PROGRAMMATICALLY
    # ============================================================================

    # Create from scratch
    annot = pt.ProFormaAnnotation(sequence="PEPTIDE", charge=2)
    print(f"New annotation: {annot.serialize()}")

    # Set internal Mods... it takes a dict of position -> {mod: count}
    annot = pt.ProFormaAnnotation(
        sequence="PEPTIDE", charge=2, internal_mods={2: {"Oxidation": 1}}
    )
    print(f"New annotation: {annot.serialize()}")

    # Other modications are just {mod: count}
    annot = pt.ProFormaAnnotation(
        sequence="PEPTIDE",
        nterm_mods={"Acetyl": 1},
        internal_mods={2: {"Oxidation": 1, "Phospho": 1}},
        charge=2,
    )
    print(f"New annotation: {annot.serialize()}")

    # ============================================================================
    # ACCESSING PROPERTIES
    # ============================================================================

    print("\n" + "=" * 60)
    print("ACCESSING PROPERTIES")
    print("=" * 60)

    annot = pt.parse("[Acetyl]-PEM[Oxidation]TIS[Phospho]DE/2")
    print(f"Annotation: {annot.serialize()}\n")

    print(f"Sequence: {annot.sequence}")
    print(f"Length: {len(annot)}")
    print(f"Charge state: {annot.charge_state}")
    print(f"Has N-term mods: {annot.has_nterm_mods}")
    print(f"Has internal mods: {annot.has_internal_mods}")
    print(f"Has charge: {annot.has_charge}")

    # ============================================================================
    # SETTING MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("SETTING MODIFICATIONS")
    print("=" * 60)

    # Start fresh
    annot = pt.ProFormaAnnotation(sequence="PEPTIDE")

    # Set N-terminal modification
    annot.set_nterm_mods({"Acetyl": 1})
    print(f"After N-term: {annot.serialize()}")

    # Set internal modification at specific position
    annot.set_internal_mods_at_index(2, {"Oxidation": 1})
    print(f"After internal: {annot.serialize()}")

    # Set charge
    annot.set_charge(2)
    print(f"After charge: {annot.serialize()}")

    # ============================================================================
    # APPENDING MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("APPENDING MODIFICATIONS")
    print("=" * 60)

    annot = pt.ProFormaAnnotation(sequence="PEPTIDE")

    # Append N-terminal mod
    annot.append_nterm_mod("Acetyl")
    print(f"Append N-term: {annot.serialize()}")

    # Append internal mod
    annot.append_internal_mod_at_index(2, "Oxidation")
    print(f"Append internal: {annot.serialize()}")

    # Append another internal mod at same position
    annot.append_internal_mod_at_index(2, "Phospho")
    print(f"Append another: {annot.serialize()}")

    # ============================================================================
    # EXTENDING MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("EXTENDING MODIFICATIONS")
    print("=" * 60)

    annot = pt.ProFormaAnnotation(sequence="PEPTIDE")

    # Extend with multiple mods
    annot.extend_nterm_mods(["Acetyl", "Formyl"])
    print(f"Extend N-term: {annot.serialize()}")

    # Extend internal mods at position
    annot.extend_internal_mods_at_index(2, ["Oxidation", "Phospho"])
    print(f"Extend internal: {annot.serialize()}")

    # ============================================================================
    # REMOVING MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("REMOVING MODIFICATIONS")
    print("=" * 60)

    annot = pt.parse("[Acetyl]-PEM[Oxidation]TIS[Phospho]DE/2")
    print(f"Original: {annot.serialize()}")

    # Clear specific mod type (removes all)
    annot.clear_nterm_mods()
    print(f"Clear N-term: {annot.serialize()}")

    # Clear internal mod at position (removes all at that position)
    annot.clear_internal_mod_at_index(2)
    print(f"Clear position 2: {annot.serialize()}")

    # Clear all mods
    annot.clear_mods()
    print(f"Clear all: {annot.serialize()}")

    # ============================================================================
    # DECREMENTING MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("DECREMENTING MODIFICATIONS")
    print("=" * 60)

    # When you have multiple copies of a mod, remove() decrements the count
    annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
    annot.extend_nterm_mods(["Acetyl", "Acetyl", "Formyl"])
    print(f"Original: {annot.serialize()}")

    # Remove one Acetyl (decrements count)
    annot.remove_nterm_mod("Acetyl")
    print(f"After removing 1 Acetyl: {annot.serialize()}")

    # Remove another Acetyl
    annot.remove_nterm_mod("Acetyl")
    print(f"After removing another Acetyl: {annot.serialize()}")

    # Remove Formyl
    annot.remove_nterm_mod("Formyl")
    print(f"After removing Formyl: {annot.serialize()}")

    # Works with internal mods too
    annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
    annot.extend_internal_mods_at_index(2, ["Oxidation", "Oxidation", "Phospho"])
    print(f"\nWith internal mods: {annot.serialize()}")

    annot.remove_internal_mod_at_index(2, "Oxidation")
    print(f"Remove 1 Oxidation: {annot.serialize()}")

    annot.remove_internal_mod_at_index(2, "Oxidation")
    print(f"Remove another Oxidation: {annot.serialize()}")

    # ============================================================================
    # POPPING MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("POPPING MODIFICATIONS")
    print("=" * 60)

    annot = pt.parse("[Acetyl]-PEM[Oxidation]TIDE/2")
    print(f"Original: {annot.serialize()}")

    # Pop N-term mods (returns the mods)
    nterm = annot.pop_nterm_mods()
    print(f"Popped N-term: {nterm}")
    print(f"After pop: {annot.serialize()}")

    # Pop charge
    charge = annot.pop_charge()
    print(f"Popped charge: {charge}")
    print(f"After pop charge: {annot.serialize()}")

    # ============================================================================
    # WORKING WITH STATIC MODS
    # ============================================================================

    print("\n" + "=" * 60)
    print("STATIC MODIFICATIONS")
    print("=" * 60)

    # Add static mod by residue
    annot = pt.ProFormaAnnotation(sequence="PEPTIDE")
    annot.add_static_mod_by_residue("E", "Oxidation")
    print(f"Static mod on E: {annot.serialize()}")

    # Condense to internal mods
    annot.condense_static_mods()
    print(f"Condensed: {annot.serialize()}")

    # ============================================================================
    # SLICING
    # ============================================================================

    print("\n" + "=" * 60)
    print("SLICING")
    print("=" * 60)

    annot = pt.parse("[Acetyl]-PEM[Oxidation]TIDE")
    print(f"Original: {annot.serialize()}")

    # Slice using indices
    sub = annot[2:5]
    print(f"Slice [2:5]: {sub.serialize()}")

    # Slice preserves modifications
    sub_with_mod = annot[1:4]
    print(f"Slice [1:4]: {sub_with_mod.serialize()}")

    # ============================================================================
    # COPYING
    # ============================================================================

    print("\n" + "=" * 60)
    print("COPYING")
    print("=" * 60)

    annot = pt.parse("PEM[Oxidation]TIDE")
    print(f"Original: {annot.serialize()}")

    # Make a copy
    copy = annot.copy()
    copy.append_nterm_mod("Acetyl")
    print(f"Copy modified: {copy.serialize()}")
    print(f"Original unchanged: {annot.serialize()}")

    # ============================================================================
    # CHECKING MODIFICATIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("CHECKING FOR MODIFICATIONS")
    print("=" * 60)

    annot = pt.parse("[Acetyl]-PEM[Oxidation]TIDE/2")
    print(f"Annotation: {annot.serialize()}\n")

    # Check for specific mod types
    print(f"Has N-term mods: {annot.has_nterm_mods}")
    print(f"Has C-term mods: {annot.has_cterm_mods}")
    print(f"Has internal mods: {annot.has_internal_mods}")
    print(f"Has charge: {annot.has_charge}")

    # Check if has any mods
    print(f"Has any mods: {annot.has_mods()}")
    print(
        f"Has internal/charge: {annot.has_mods([pt.ModType.INTERNAL, pt.ModType.CHARGE])}"
    )
    print(f"Has internal/charge: {annot.has_mods(['internal', 'charge'])}")

    # ============================================================================
    # SERIALIZATION OPTIONS
    # ============================================================================

    print("\n" + "=" * 60)
    print("SERIALIZATION")
    print("=" * 60)

    annot = pt.parse("[Acetyl]-PEPTIDE/2")

    # Full serialization
    print(f"Full: {annot.serialize()}")

    # Strip mods before serializing
    stripped = annot.copy().strip_mods()
    print(f"Stripped: {stripped.serialize()}")

    # ============================================================================
    # VALIDATION
    # ============================================================================

    print("\n" + "=" * 60)
    print("VALIDATION")
    print("=" * 60)

    # By default, validation is OFF for performance
    annot_no_val = pt.ProFormaAnnotation(sequence="PEPTIDE", validate=False)
    print(f"No validation (default): {annot_no_val.serialize()}")

    # Enable validation when creating (and for methods that modify the annotation)
    annot_with_val = pt.ProFormaAnnotation(sequence="PEPTIDE", validate=True)
    print(f"With validation: {annot_with_val.serialize()}")

    # Validation checks modification syntax (can be disabled per method)
    print("\nAttempting to add invalid modification with validation ON:")
    try:
        annot_with_val.append_internal_mod_at_index(0, "InvalidMod123")
        print(f"  Error: Modification added unexpectedly: {annot_with_val.serialize()}")
    except Exception as e:
        # successfully raises error
        print(f"Successfully caught error: {e}")

    # Validation checks modification syntax (can be disabled per method)
    print("\nAttempting to add invalid modification with validation OFF:")
    try:
        annot_with_val.append_internal_mod_at_index(0, "InvalidMod123", validate=False)
        print(f"  Success (no validation): {annot_no_val.serialize()}")
    except Exception as e:
        print(f"  Error: {e}")


if __name__ == "__main__":
    run()
