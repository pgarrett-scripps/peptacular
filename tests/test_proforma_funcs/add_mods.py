# filepath: /workspaces/peptacular/tests/test_proforma_funcs/add_mods.py
import unittest

import peptacular as pt
from peptacular.constants import ModType
from peptacular.proforma_dataclasses import Mod, Interval


class TestAddMods(unittest.TestCase):
    
    """
    TESTS FOR: add_labile_mods
    """
    def test_add_labile_mods_replace_inplace_true(self):
        annotation = pt.parse("{Glycan}PEPTIDE")
        result = annotation.add_labile_mods("HexNAc", append=False, inplace=True)
        
        self.assertIs(result, annotation)  # Should return same object
        self.assertEqual(annotation.serialize(), "{HexNAc}PEPTIDE")

    def test_add_labile_mods_replace_inplace_false(self):
        annotation = pt.parse("{Glycan}PEPTIDE")
        result = annotation.add_labile_mods("HexNAc", append=False, inplace=False)
        
        self.assertIsNot(result, annotation)  # Should return different object
        self.assertEqual(result.serialize(), "{HexNAc}PEPTIDE")
        # Original should be unchanged
        self.assertEqual(annotation.serialize(), "{Glycan}PEPTIDE")

    def test_add_labile_mods_append_to_existing(self):
        annotation = pt.parse("{Glycan}PEPTIDE")
        annotation.add_labile_mods("HexNAc", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "{Glycan}{HexNAc}PEPTIDE")

    def test_add_labile_mods_append_to_empty(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_labile_mods("Glycan", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "{Glycan}PEPTIDE")

    def test_add_labile_mods_list_input(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_labile_mods(["Glycan", "HexNAc"], append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "{Glycan}{HexNAc}PEPTIDE")

    def test_add_labile_mods_mod_object_input(self):
        annotation = pt.parse("PEPTIDE")
        mod = Mod("Glycan", 1)
        annotation.add_labile_mods(mod, append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "{Glycan}PEPTIDE")

    def test_add_labile_mods_none_input_no_append(self):
        annotation = pt.parse("{Glycan}PEPTIDE")
        annotation.add_labile_mods(None, append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_add_labile_mods_none_input_append(self):
        annotation = pt.parse("{Glycan}PEPTIDE")
        annotation.add_labile_mods(None, append=True, inplace=True)
        
        # Should remain unchanged when appending None
        self.assertEqual(annotation.serialize(), "{Glycan}PEPTIDE")

    """
    TESTS FOR: add_unknown_mods
    """
    def test_add_unknown_mods_replace_inplace_true(self):
        annotation = pt.parse("[Unknown]?PEPTIDE")
        result = annotation.add_unknown_mods("Mystery", append=False, inplace=True)
        
        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "[Mystery]?PEPTIDE")

    def test_add_unknown_mods_replace_inplace_false(self):
        annotation = pt.parse("[Unknown]?PEPTIDE")
        result = annotation.add_unknown_mods("Mystery", append=False, inplace=False)
        
        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "[Mystery]?PEPTIDE")
        # Original should be unchanged
        self.assertEqual(annotation.serialize(), "[Unknown]?PEPTIDE")

    def test_add_unknown_mods_append_to_existing(self):
        annotation = pt.parse("[Unknown]?PEPTIDE")
        annotation.add_unknown_mods("Mystery", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Unknown][Mystery]?PEPTIDE")

    def test_add_unknown_mods_append_to_empty(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_unknown_mods("Unknown", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Unknown]?PEPTIDE")

    def test_add_unknown_mods_list_input(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_unknown_mods(["Unknown1", "Unknown2"], append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Unknown1][Unknown2]?PEPTIDE")

    def test_add_unknown_mods_none_input_no_append(self):
        annotation = pt.parse("[Unknown]?PEPTIDE")
        annotation.add_unknown_mods(None, append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_add_unknown_mods_none_input_append(self):
        annotation = pt.parse("[Unknown]?PEPTIDE")
        annotation.add_unknown_mods(None, append=True, inplace=True)
        
        # Should remain unchanged when appending None
        self.assertEqual(annotation.serialize(), "[Unknown]?PEPTIDE")

    """
    TESTS FOR: add_nterm_mods
    """
    def test_add_nterm_mods_replace_inplace_true(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        result = annotation.add_nterm_mods("Formyl", append=False, inplace=True)
        
        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "[Formyl]-PEPTIDE")

    def test_add_nterm_mods_replace_inplace_false(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        result = annotation.add_nterm_mods("Formyl", append=False, inplace=False)
        
        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "[Formyl]-PEPTIDE")
        # Original should be unchanged
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_add_nterm_mods_append_to_existing(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        annotation.add_nterm_mods("Formyl", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Acetyl][Formyl]-PEPTIDE")

    def test_add_nterm_mods_append_to_empty(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_nterm_mods("Acetyl", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_add_nterm_mods_list_input(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_nterm_mods(["Acetyl", "Formyl"], append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Acetyl][Formyl]-PEPTIDE")

    def test_add_nterm_mods_none_input_no_append(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        annotation.add_nterm_mods(None, append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_add_nterm_mods_none_input_append(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        annotation.add_nterm_mods(None, append=True, inplace=True)
        
        # Should remain unchanged when appending None
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    """
    TESTS FOR: add_cterm_mods
    """
    def test_add_cterm_mods_replace_inplace_true(self):
        annotation = pt.parse("PEPTIDE-[Amide]")
        result = annotation.add_cterm_mods("Hydroxyl", append=False, inplace=True)
        
        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "PEPTIDE-[Hydroxyl]")

    def test_add_cterm_mods_replace_inplace_false(self):
        annotation = pt.parse("PEPTIDE-[Amide]")
        result = annotation.add_cterm_mods("Hydroxyl", append=False, inplace=False)
        
        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "PEPTIDE-[Hydroxyl]")
        # Original should be unchanged
        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide]")

    def test_add_cterm_mods_append_to_existing(self):
        annotation = pt.parse("PEPTIDE-[Amide]")
        annotation.add_cterm_mods("Hydroxyl", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide][Hydroxyl]")

    def test_add_cterm_mods_append_to_empty(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_cterm_mods("Amide", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide]")

    def test_add_cterm_mods_list_input(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_cterm_mods(["Amide", "Hydroxyl"], append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide][Hydroxyl]")

    def test_add_cterm_mods_none_input_no_append(self):
        annotation = pt.parse("PEPTIDE-[Amide]")
        annotation.add_cterm_mods(None, append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_add_cterm_mods_none_input_append(self):
        annotation = pt.parse("PEPTIDE-[Amide]")
        annotation.add_cterm_mods(None, append=True, inplace=True)
        
        # Should remain unchanged when appending None
        self.assertEqual(annotation.serialize(), "PEPTIDE-[Amide]")

    """
    TESTS FOR: add_mod_by_type
    """
    def test_add_mod_by_type_string_parameter(self):
        annotation = pt.parse("PEPTIDE")
        result = annotation.add_mod_by_type("Acetyl", "nterm", append=False, inplace=True)
        
        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_add_mod_by_type_mod_type_parameter(self):
        annotation = pt.parse("PEPTIDE")
        result = annotation.add_mod_by_type("Acetyl", ModType.NTERM, append=False, inplace=True)
        
        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE")

    def test_add_mod_by_type_inplace_false(self):
        annotation = pt.parse("PEPTIDE")
        result = annotation.add_mod_by_type("Acetyl", "nterm", append=False, inplace=False)
        
        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "[Acetyl]-PEPTIDE")
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_add_mod_by_type_all_types(self):
        # Test each mod type
        test_cases = [
            ("Glycan", "labile", "{Glycan}PEPTIDE"),
            ("Unknown", "unknown", "[Unknown]?PEPTIDE"),
            ("Acetyl", "nterm", "[Acetyl]-PEPTIDE"),
            ("Amide", "cterm", "PEPTIDE-[Amide]"),
            ("15N", "isotope", "<15N>PEPTIDE"),
            ("57@C", "static", "<57@C>PEPTIDE"),
        ]
        
        for mod_val, mod_type, expected_serialization in test_cases:
            with self.subTest(mod_type=mod_type):
                annotation = pt.parse("PEPTIDE")
                annotation.add_mod_by_type(mod_val, mod_type, append=False, inplace=True)
                self.assertEqual(annotation.serialize(), expected_serialization)

    def test_add_mod_by_type_charge(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_mod_by_type(2, "charge", append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE/2")

    def test_add_mod_by_type_charge_adducts(self):
        annotation = pt.parse("PEPTIDE/2")
        annotation.add_mod_by_type("Na+", "charge_adducts", append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE/2[Na+]")

    def test_add_mod_by_type_internal(self):
        annotation = pt.parse("PEPTIDE")
        annotation.add_mod_by_type("Phospho", "internal", append=False, inplace=True)
        
        # Internal mods are added at position 0 by default in add_mod_by_type
        self.assertEqual(annotation.serialize(), "P[Phospho]EPTIDE")

    def test_add_mod_by_type_intervals(self):
        annotation = pt.parse("PEPTIDE")
        interval = Interval(1, 3, False, [Mod("Phospho", 1)])
        annotation.add_mod_by_type(interval, "interval", append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "P(EP)[Phospho]TIDE")

    def test_add_mod_by_type_append_behavior(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        annotation.add_mod_by_type("Formyl", "nterm", append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Acetyl][Formyl]-PEPTIDE")

    """
    TESTS FOR: add_mod_dict
    """
    def test_add_mod_dict_inplace_true(self):
        annotation = pt.parse("PEPTIDE")
        mod_dict = {
            "nterm": "Acetyl",
            "cterm": "Amide",
            "labile": "Glycan"
        }
        result = annotation.add_mod_dict(mod_dict, append=False, inplace=True)
        
        self.assertIs(result, annotation)
        self.assertEqual(annotation.serialize(), "{Glycan}[Acetyl]-PEPTIDE-[Amide]")

    def test_add_mod_dict_inplace_false(self):
        annotation = pt.parse("PEPTIDE")
        mod_dict = {
            "nterm": "Acetyl",
            "cterm": "Amide"
        }
        result = annotation.add_mod_dict(mod_dict, append=False, inplace=False)
        
        self.assertIsNot(result, annotation)
        self.assertEqual(result.serialize(), "[Acetyl]-PEPTIDE-[Amide]")
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_add_mod_dict_with_internal_mods(self):
        annotation = pt.parse("PEPTIDE")
        mod_dict = {
            "nterm": "Acetyl",
            1: "Phospho",  # Internal mod at position 1
            3: "Methyl"    # Internal mod at position 3
        }
        annotation.add_mod_dict(mod_dict, append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Acetyl]-PE[Phospho]PT[Methyl]IDE")

    def test_add_mod_dict_append_behavior(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        mod_dict = {
            "nterm": "Formyl",
            "cterm": "Amide"
        }
        annotation.add_mod_dict(mod_dict, append=True, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Acetyl][Formyl]-PEPTIDE-[Amide]")

    def test_add_mod_dict_replace_behavior(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE-[Amide]")
        mod_dict = {
            "nterm": "Formyl",
            "cterm": "Hydroxyl"
        }
        annotation.add_mod_dict(mod_dict, append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Formyl]-PEPTIDE-[Hydroxyl]")

    """
    TESTS FOR: edge cases and validation
    """
    def test_add_methods_preserve_sequence(self):
        original_sequence = "PEPTIDE"
        
        # Test each add method preserves the sequence
        test_cases = [
            ("add_labile_mods", "Glycan"),
            ("add_unknown_mods", "Unknown"),
            ("add_nterm_mods", "Acetyl"),
            ("add_cterm_mods", "Amide"),
        ]
        
        for method_name, mod_val in test_cases:
            with self.subTest(method=method_name):
                annotation = pt.parse("PEPTIDE")
                method = getattr(annotation, method_name)
                method(mod_val, append=False, inplace=True)
                self.assertEqual(annotation.sequence, original_sequence)

    def test_add_methods_return_self_when_inplace(self):
        annotation = pt.parse("PEPTIDE")
        
        methods = [
            ("add_labile_mods", "Glycan"),
            ("add_unknown_mods", "Unknown"),
            ("add_nterm_mods", "Acetyl"),
            ("add_cterm_mods", "Amide"),
        ]
        
        for method_name, mod_val in methods:
            with self.subTest(method=method_name):
                method = getattr(annotation, method_name)
                result = method(mod_val, append=False, inplace=True)
                self.assertIs(result, annotation)

    def test_add_methods_return_copy_when_not_inplace(self):
        annotation = pt.parse("PEPTIDE")
        
        methods = [
            ("add_labile_mods", "Glycan"),
            ("add_unknown_mods", "Unknown"),
            ("add_nterm_mods", "Acetyl"),
            ("add_cterm_mods", "Amide"),
        ]
        
        for method_name, mod_val in methods:
            with self.subTest(method=method_name):
                method = getattr(annotation, method_name)
                result = method(mod_val, append=False, inplace=False)
                self.assertIsNot(result, annotation)

    def test_add_methods_with_multipliers(self):
        annotation = pt.parse("PEPTIDE")
        mod_with_multiplier = Mod("Phospho", 2)
        
        annotation.add_nterm_mods(mod_with_multiplier, append=False, inplace=True)
        
        self.assertEqual(annotation.serialize(), "[Phospho]^2-PEPTIDE")

    def test_add_empty_list_to_existing_mods(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        original_serialization = annotation.serialize()
        
        annotation.add_nterm_mods([], append=True, inplace=True)
        
        # Adding empty list should not change the serialization
        self.assertEqual(annotation.serialize(), original_serialization)

    def test_complex_mod_dict_scenario(self):
        annotation = pt.parse("PEPTIDE")
        complex_mod_dict = {
            "labile": ["Glycan", "HexNAc"],
            "nterm": "Acetyl",
            "cterm": ["Amide", "Hydroxyl"],
            0: "Phospho",  # Internal mod at position 0
            2: ["Methyl", "Oxidation"],  # Multiple internal mods at position 2
            "charge": 2,
            "isotope": "15N"
        }
        
        annotation.add_mod_dict(complex_mod_dict, append=False, inplace=True)
        
        # Verify the complex modification produces expected serialization
        self.assertEqual(annotation.serialize(), "{Glycan}{HexNAc}<15N>[Acetyl]-P[Phospho]EP[Methyl][Oxidation]TIDE-[Amide][Hydroxyl]/2")


if __name__ == '__main__':
    unittest.main()