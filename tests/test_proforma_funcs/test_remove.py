# filepath: /workspaces/peptacular/tests/test_proforma_funcs/test_remove.py
import unittest

import peptacular as pt


class TestRemoveMods(unittest.TestCase):
    """
    TESTS FOR: remove_nterm_mods
    """

    def test_remove_nterm_mods_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_nterm_mods(inplace=True)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE-[Amide]")
        self.assertEqual(
            annotation.serialize(), "PEPTIDE-[Amide]"
        )  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_nterm_mods_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_nterm_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE-[Amide]")
        self.assertEqual(
            annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]"
        )  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    def test_remove_nterm_mods_no_nterm_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_nterm_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE-[Amide]")

    def test_remove_nterm_mods_multiple_nterm_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl][Phospho]-PEPTIDE")
        removed_annotation = annotation.remove_nterm_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: remove_cterm_mods
    """

    def test_remove_cterm_mods_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_cterm_mods(inplace=True)
        self.assertEqual(removed_annotation.serialize(), "[Acetyl]-PEPTIDE")
        self.assertEqual(
            annotation.serialize(), "[Acetyl]-PEPTIDE"
        )  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_cterm_mods_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_cterm_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "[Acetyl]-PEPTIDE")
        self.assertEqual(
            annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]"
        )  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    def test_remove_cterm_mods_no_cterm_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        removed_annotation = annotation.remove_cterm_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_remove_cterm_mods_multiple_cterm_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide][Phospho]")
        removed_annotation = annotation.remove_cterm_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: remove_internal_mods
    """

    def test_remove_internal_mods_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        removed_annotation = annotation.remove_internal_mods(inplace=True)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(annotation.serialize(), "PEPTIDE")  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_internal_mods_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        removed_annotation = annotation.remove_internal_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(
            annotation.serialize(), "PE[Phospho]PTI[Methyl]DE"
        )  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    def test_remove_internal_mods_no_internal_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        removed_annotation = annotation.remove_internal_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_internal_mods_multiple_mods_same_position(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho][Methyl]PTIDE")
        removed_annotation = annotation.remove_internal_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: remove_intervals
    """

    def test_remove_intervals_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        removed_annotation = annotation.remove_intervals(inplace=True)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(annotation.serialize(), "PEPTIDE")  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_intervals_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        removed_annotation = annotation.remove_intervals(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(
            annotation.serialize(), "P(EP)[Phospho]TIDE"
        )  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    def test_remove_intervals_no_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        removed_annotation = annotation.remove_intervals(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_intervals_ambiguous(self):
        annotation = pt.ProFormaAnnotation.parse("P(?EP)TIDE")
        removed_annotation = annotation.remove_intervals(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_intervals_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]T(ID)[Methyl]E")
        removed_annotation = annotation.remove_intervals(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: remove_charge
    """

    def test_remove_charge_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        removed_annotation = annotation.remove_charge(inplace=True)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(annotation.serialize(), "PEPTIDE")  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_charge_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        removed_annotation = annotation.remove_charge(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(annotation.serialize(), "PEPTIDE/2")  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    def test_remove_charge_no_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        removed_annotation = annotation.remove_charge(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_charge_negative_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/-2")
        removed_annotation = annotation.remove_charge(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: remove_charge_adducts
    """

    def test_remove_charge_adducts_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2[Na+]")
        removed_annotation = annotation.remove_charge_adducts(inplace=True)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE/2")
        self.assertEqual(annotation.serialize(), "PEPTIDE/2")  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_charge_adducts_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2[Na+]")
        removed_annotation = annotation.remove_charge_adducts(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE/2")
        self.assertEqual(annotation.serialize(), "PEPTIDE/2[Na+]")  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    def test_remove_charge_adducts_no_charge_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        removed_annotation = annotation.remove_charge_adducts(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE/2")

    def test_remove_charge_adducts_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2[Na+][K+]")
        removed_annotation = annotation.remove_charge_adducts(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE/2")

    """
    TESTS FOR: remove_isotope_mods
    """

    def test_remove_isotope_mods_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        removed_annotation = annotation.remove_isotope_mods(inplace=True)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(annotation.serialize(), "PEPTIDE")  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_isotope_mods_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        removed_annotation = annotation.remove_isotope_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(annotation.serialize(), "<15N>PEPTIDE")  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    def test_remove_isotope_mods_no_isotope_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        removed_annotation = annotation.remove_isotope_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_isotope_mods_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("<15N><13C>PEPTIDE")
        removed_annotation = annotation.remove_isotope_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: remove_static_mods
    """

    def test_remove_static_mods_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        removed_annotation = annotation.remove_static_mods(inplace=True)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(annotation.serialize(), "PEPTIDE")  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_static_mods_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        removed_annotation = annotation.remove_static_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")
        self.assertEqual(annotation.serialize(), "<57@C>PEPTIDE")  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    def test_remove_static_mods_no_static_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        removed_annotation = annotation.remove_static_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_static_mods_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("<57@C><15@M>PEPTIDE")
        removed_annotation = annotation.remove_static_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    """
    TESTS FOR: remove_mods (general method)
    """

    def test_remove_mods_all_default(self):
        annotation = pt.ProFormaAnnotation.parse(
            "{Glycan}<57@C>[Acetyl]-PE[Phospho]PTIDE-[Amide]/2"
        )
        removed_annotation = annotation.remove_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_mods_specific_types_single_string(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        removed_annotation = annotation.remove_mods(mods="nterm", inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PE[Phospho]PTIDE-[Amide]")

    def test_remove_mods_specific_types_list(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        removed_annotation = annotation.remove_mods(
            mods=["nterm", "cterm"], inplace=False
        )
        self.assertEqual(removed_annotation.serialize(), "PE[Phospho]PTIDE")

    def test_remove_mods_mod_type_enum(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        removed_annotation = annotation.remove_mods(
            mods=pt.ModType.INTERNAL, inplace=False
        )
        self.assertEqual(removed_annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")

    def test_remove_mods_mixed_list(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")
        removed_annotation = annotation.remove_mods(
            mods=["nterm", pt.ModType.CHARGE], inplace=False
        )
        self.assertEqual(removed_annotation.serialize(), "PE[Phospho]PTIDE-[Amide]")

    def test_remove_mods_inplace_true(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_mods(mods=["nterm"], inplace=True)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE-[Amide]")
        self.assertEqual(
            annotation.serialize(), "PEPTIDE-[Amide]"
        )  # Original changed too
        self.assertIs(removed_annotation, annotation)

    def test_remove_mods_inplace_false(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_mods(mods=["nterm"], inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE-[Amide]")
        self.assertEqual(
            annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]"
        )  # Original unchanged
        self.assertIsNot(removed_annotation, annotation)

    """
    TESTS FOR: edge cases
    """

    def test_remove_mods_no_mods_to_remove(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        removed_annotation = annotation.remove_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_mods_empty_list(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_mods(mods=[], inplace=False)
        self.assertEqual(removed_annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")

    def test_remove_mods_none_parameter(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_mods(mods=None, inplace=False)
        self.assertEqual(removed_annotation.serialize(), "PEPTIDE")

    def test_remove_methods_with_labile_mods(self):
        # Test that labile mods are preserved when removing other mod types
        annotation = pt.ProFormaAnnotation.parse("{Glycan}[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_nterm_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "{Glycan}PEPTIDE-[Amide]")
        self.assertTrue(removed_annotation.has_labile_mods)

    def test_remove_methods_with_unknown_mods(self):
        # Test that unknown mods are preserved when removing other mod types
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?[Acetyl]-PEPTIDE-[Amide]")
        removed_annotation = annotation.remove_nterm_mods(inplace=False)
        self.assertEqual(removed_annotation.serialize(), "[Unknown]?PEPTIDE-[Amide]")
        self.assertTrue(removed_annotation.has_unknown_mods)


if __name__ == "__main__":
    unittest.main()
