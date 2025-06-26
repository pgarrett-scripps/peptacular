# filepath: /workspaces/peptacular/tests/test_proforma_funcs/test_filter.py
import unittest

import peptacular as pt
from peptacular.constants import ModType


class TestFilterMods(unittest.TestCase):
    
    def test_filter_mods_basic_inplace_true(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['nterm'], inplace=True)
        # Should keep only nterm mods
        self.assertEqual(filtered_annotation.serialize(), "[Acetyl]-PEPTIDE")
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")  # Original changed too

    def test_filter_mods_basic_inplace_false(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['nterm'], inplace=False)
        # Should keep only nterm mods
        self.assertEqual(filtered_annotation.serialize(), "[Acetyl]-PEPTIDE")
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")  # Original unchanged

    def test_filter_mods_single_string(self):
        annotation = pt.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods='internal', inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "PE[Phospho]PTIDE")

    def test_filter_mods_single_mod_type(self):
        annotation = pt.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=ModType.INTERNAL, inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "PE[Phospho]PTIDE")

    def test_filter_mods_list_of_strings(self):
        annotation = pt.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['nterm', 'cterm'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")

    def test_filter_mods_list_of_mod_types(self):
        annotation = pt.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=[ModType.NTERM, ModType.CTERM], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")

    def test_filter_mods_mixed_list(self):
        annotation = pt.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['nterm', ModType.INTERNAL], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "[Acetyl]-PE[Phospho]PTIDE")

    """
    TESTS FOR: filtering different modification types
    """
    def test_filter_mods_keep_labile(self):
        annotation = pt.parse("{Glycan}[Acetyl]-PEPTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['labile'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "{Glycan}PEPTIDE")

    def test_filter_mods_keep_static(self):
        annotation = pt.parse("<57@C>[Acetyl]-PEPTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['static'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "<57@C>PEPTIDE")

    def test_filter_mods_keep_isotope(self):
        annotation = pt.parse("<15N>[Acetyl]-PEPTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['isotope'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "<15N>PEPTIDE")

    def test_filter_mods_keep_unknown(self):
        annotation = pt.parse("[Unknown]?[Acetyl]-PEPTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['unknown'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "[Unknown]?PEPTIDE")

    def test_filter_mods_keep_charge(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE-[Amide]/2")
        filtered_annotation = annotation.filter_mods(mods=['charge'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "PEPTIDE/2")

    def test_filter_mods_keep_charge_adducts(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE-[Amide]/2[Na+]")
        filtered_annotation = annotation.filter_mods(mods=['charge_adducts'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "PEPTIDE[Na+]")

    def test_filter_mods_keep_intervals(self):
        annotation = pt.parse("[Acetyl]-P(EP)[Phospho]TIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['interval'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "P(EP)[Phospho]TIDE")

    """
    TESTS FOR: filtering multiple modification types
    """
    def test_filter_mods_keep_multiple_types(self):
        annotation = pt.parse("{Glycan}<57@C>[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")
        filtered_annotation = annotation.filter_mods(mods=['nterm', 'internal', 'charge'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "[Acetyl]-PE[Phospho]PTIDE/2")

    def test_filter_mods_keep_all_terminal_mods(self):
        annotation = pt.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['nterm', 'cterm'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")

    """
    TESTS FOR: edge cases
    """
    def test_filter_mods_no_mods_to_filter(self):
        annotation = pt.parse("PEPTIDE")
        filtered_annotation = annotation.filter_mods(mods=['nterm'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "PEPTIDE")

    def test_filter_mods_none_parameter(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=None, inplace=False)
        # Should remove all mods when mods=None (empty list)
        self.assertEqual(filtered_annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")

    def test_filter_mods_empty_list(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=[], inplace=False)
        # Should remove all mods when mods=[] (empty list)
        self.assertEqual(filtered_annotation.serialize(), "PEPTIDE")

    def test_filter_mods_all_mod_types(self):
        # Test keeping all modification types (should not change anything)
        annotation = pt.parse("{Glycan}<57@C>[Unknown]?[Acetyl]-PE[Phospho]PTIDE-[Amide]/2[Na+]")
        all_mod_types = ['labile', 'static', 'unknown', 'nterm', 'cterm', 'internal', 'charge', 'charge_adducts']
        filtered_annotation = annotation.filter_mods(mods=all_mod_types, inplace=False)
        self.assertEqual(filtered_annotation.serialize(), annotation.serialize())

    def test_filter_mods_complex_annotation(self):
        annotation = pt.parse("{Glycan}<57@C>[Acetyl]-P(EP)[Phospho]TI[Methyl]DE-[Amide]/2")
        # Keep only internal modifications
        filtered_annotation = annotation.filter_mods(mods=['internal'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "PEPTI[Methyl]DE")

    def test_filter_mods_with_intervals_ambiguous(self):
        annotation = pt.parse("[Acetyl]-P(?EP)TIDE-[Amide]")
        filtered_annotation = annotation.filter_mods(mods=['interval'], inplace=False)
        self.assertEqual(filtered_annotation.serialize(), "P(?EP)TIDE")


if __name__ == '__main__':
    unittest.main()