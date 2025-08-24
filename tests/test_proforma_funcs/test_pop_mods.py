# filepath: /workspaces/peptacular/tests/test_proforma_funcs/test_pop_mods.py
import unittest

import peptacular as pt


class TestPopMods(unittest.TestCase):
    
    """
    TESTS FOR: individual pop methods
    """
    def test_pop_labile_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        popped_mods = annotation.pop_labile_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0].val, "Glycan")
        self.assertFalse(annotation.has_labile_mods())
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_labile_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_labile_mods()
        
        self.assertEqual(len(popped_mods), 0)
        self.assertFalse(annotation.has_labile_mods())

    def test_pop_labile_mods_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}{HexNAc}PEPTIDE")
        popped_mods = annotation.pop_labile_mods()
        
        self.assertEqual(len(popped_mods), 2)
        self.assertEqual([mod.val for mod in popped_mods], ["Glycan", "HexNAc"])
        self.assertFalse(annotation.has_labile_mods())

    def test_pop_unknown_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        popped_mods = annotation.pop_unknown_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0].val, "Unknown")
        self.assertFalse(annotation.has_unknown_mods())

    def test_pop_unknown_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_unknown_mods()
        
        self.assertEqual(len(popped_mods), 0)
        self.assertFalse(annotation.has_unknown_mods())

    def test_pop_nterm_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_nterm_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0].val, "Acetyl")
        self.assertFalse(annotation.has_nterm_mods())
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_nterm_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_nterm_mods()
        
        self.assertEqual(len(popped_mods), 0)
        self.assertFalse(annotation.has_nterm_mods())

    def test_pop_nterm_mods_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl][Phospho]-PEPTIDE")
        popped_mods = annotation.pop_nterm_mods()
        
        self.assertEqual(len(popped_mods), 2)
        self.assertEqual([mod.val for mod in popped_mods], ["Acetyl", "Phospho"])
        self.assertFalse(annotation.has_nterm_mods())

    def test_pop_cterm_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        popped_mods = annotation.pop_cterm_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0].val, "Amide")
        self.assertFalse(annotation.has_cterm_mods())
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_cterm_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_cterm_mods()
        
        self.assertEqual(len(popped_mods), 0)
        self.assertFalse(annotation.has_cterm_mods())

    def test_pop_internal_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        popped_mods = annotation.pop_internal_mods()
        
        self.assertEqual(len(popped_mods), 2)
        self.assertIn(1, popped_mods)  # Position 1 (E)
        self.assertIn(4, popped_mods)  # Position 4 (I)
        self.assertEqual(popped_mods[1][0].val, "Phospho")
        self.assertEqual(popped_mods[4][0].val, "Methyl")
        self.assertFalse(annotation.has_internal_mods())
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_internal_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_internal_mods()
        
        self.assertEqual(len(popped_mods), 0)
        self.assertFalse(annotation.has_internal_mods())

    def test_pop_intervals_with_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        popped_intervals = annotation.pop_intervals()
        
        self.assertEqual(len(popped_intervals), 1)
        self.assertEqual(popped_intervals[0].start, 1)
        self.assertEqual(popped_intervals[0].end, 3)
        self.assertFalse(annotation.has_intervals())
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_intervals_no_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_intervals = annotation.pop_intervals()
        
        self.assertEqual(len(popped_intervals), 0)
        self.assertFalse(annotation.has_intervals())

    def test_pop_charge_with_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        popped_charge = annotation.pop_charge()
        
        self.assertEqual(popped_charge, 2)
        self.assertFalse(annotation.has_charge())
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_charge_no_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_charge = annotation.pop_charge()
        
        self.assertIsNone(popped_charge)
        self.assertFalse(annotation.has_charge())

    def test_pop_charge_negative_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/-2")
        popped_charge = annotation.pop_charge()
        
        self.assertEqual(popped_charge, -2)
        self.assertFalse(annotation.has_charge())

    def test_pop_charge_adducts_with_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2[Na+]")
        popped_adducts = annotation.pop_charge_adducts()
        
        self.assertEqual(len(popped_adducts), 1)
        self.assertEqual(popped_adducts[0].val, "Na+")
        self.assertFalse(annotation.has_charge_adducts())
        self.assertEqual(annotation.serialize(), "PEPTIDE/2")

    def test_pop_charge_adducts_no_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        popped_adducts = annotation.pop_charge_adducts()
        
        self.assertEqual(len(popped_adducts), 0)
        self.assertFalse(annotation.has_charge_adducts())

    def test_pop_charge_adducts_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2[Na+][K+]")
        popped_adducts = annotation.pop_charge_adducts()
        
        self.assertEqual(len(popped_adducts), 2)
        self.assertEqual([mod.val for mod in popped_adducts], ["Na+", "K+"])
        self.assertFalse(annotation.has_charge_adducts())

    def test_pop_isotope_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        popped_mods = annotation.pop_isotope_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0].val, "15N")
        self.assertFalse(annotation.has_isotope_mods())
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_isotope_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_isotope_mods()
        
        self.assertEqual(len(popped_mods), 0)
        self.assertFalse(annotation.has_isotope_mods())

    def test_pop_static_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        popped_mods = annotation.pop_static_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0].val, "57@C")
        self.assertFalse(annotation.has_static_mods())
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_static_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_static_mods()
        
        self.assertEqual(len(popped_mods), 0)
        self.assertFalse(annotation.has_static_mods())

    """
    TESTS FOR: pop_mod_by_type method
    """
    def test_pop_mod_by_type_string_parameter(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_mod_by_type('nterm')
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0].val, "Acetyl")
        self.assertFalse(annotation.has_nterm_mods())

    def test_pop_mod_by_type_mod_type_parameter(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_mod_by_type(ModType.NTERM)
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0].val, "Acetyl")
        self.assertFalse(annotation.has_nterm_mods())

    def test_pop_mod_by_type_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_mod_by_type('nterm')
        
        self.assertEqual(len(popped_mods), 0)

    def test_pop_mod_by_type_all_types(self):
        # Test each mod type
        test_cases = [
            ("{Glycan}PEPTIDE", "labile", "Glycan"),
            ("[Unknown]?PEPTIDE", "unknown", "Unknown"),
            ("[Acetyl]-PEPTIDE", "nterm", "Acetyl"),
            ("PEPTIDE-[Amide]", "cterm", "Amide"),
            ("<15N>PEPTIDE", "isotope", "15N"),
            ("<57@C>PEPTIDE", "static", "57@C"),
        ]
        
        for proforma_string, mod_type, expected_val in test_cases:
            with self.subTest(mod_type=mod_type):
                annotation = pt.ProFormaAnnotation.parse(proforma_string)
                popped_mods = annotation.pop_mod_by_type(mod_type)
                self.assertEqual(len(popped_mods), 1)
                self.assertEqual(popped_mods[0].val, expected_val)

    def test_pop_mod_by_type_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        popped_charge = annotation.pop_mod_by_type('charge')
        
        self.assertEqual(popped_charge, 2)
        self.assertFalse(annotation.has_charge())

    def test_pop_mod_by_type_charge_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2[Na+]")
        popped_adducts = annotation.pop_mod_by_type('charge_adducts')
        
        self.assertEqual(len(popped_adducts), 1)
        self.assertEqual(popped_adducts[0].val, "Na+")

    def test_pop_mod_by_type_internal(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTIDE")
        popped_mods = annotation.pop_mod_by_type('internal')
        
        self.assertEqual(len(popped_mods), 1)
        self.assertIn(1, popped_mods)
        self.assertEqual(popped_mods[1][0].val, "Phospho")

    def test_pop_mod_by_type_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        popped_intervals = annotation.pop_mod_by_type('interval')
        
        self.assertEqual(len(popped_intervals), 1)
        self.assertEqual(popped_intervals[0].start, 1)

    def test_pop_mod_by_type_invalid_type(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        with self.assertRaises(ValueError):
            annotation.pop_mod_by_type('invalid_type')

    def test_pop_mod_by_type_invalid_object(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        with self.assertRaises(ValueError):
            annotation.pop_mod_by_type(123)

    """
    TESTS FOR: pop_mods method (general)
    """
    def test_pop_mods_all_default(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}<57@C>[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")
        popped_mods = annotation.pop_mods()
        
        # Should pop all modification types
        expected_keys = ['labile', 'static', 'nterm', 'internal', 'cterm', 'charge']
        for key in expected_keys:
            self.assertIn(key, popped_mods)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_mods_specific_types_single_string(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods='nterm')
        
        self.assertIn('nterm', popped_mods)
        self.assertEqual(len(popped_mods['nterm']), 1)
        self.assertEqual(popped_mods['nterm'][0].val, "Acetyl")
        
        # Other mods should still be present
        self.assertTrue(annotation.has_internal_mods())
        self.assertTrue(annotation.has_cterm_mods())

    def test_pop_mods_specific_types_list(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods=['nterm', 'cterm'])
        
        self.assertIn('nterm', popped_mods)
        self.assertIn('cterm', popped_mods)
        self.assertEqual(len(popped_mods['nterm']), 1)
        self.assertEqual(len(popped_mods['cterm']), 1)
        
        # Internal mods should still be present
        self.assertTrue(annotation.has_internal_mods())

    def test_pop_mods_mod_type_enum(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods=ModType.INTERNAL)
        
        self.assertIn('internal', popped_mods)
        self.assertEqual(len(popped_mods['internal']), 1)
        
        # Other mods should still be present
        self.assertTrue(annotation.has_nterm_mods())
        self.assertTrue(annotation.has_cterm_mods())

    def test_pop_mods_mixed_list(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")
        popped_mods = annotation.pop_mods(mods=['nterm', ModType.CHARGE])
        
        self.assertIn('nterm', popped_mods)
        self.assertIn('charge', popped_mods)
        self.assertEqual(len(popped_mods['nterm']), 1)
        self.assertEqual(popped_mods['charge'], 2)

    def test_pop_mods_condense_true(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_mods(condense=True)
        
        # Should only contain keys with actual values
        self.assertIn('nterm', popped_mods)
        self.assertNotIn('cterm', popped_mods)  # No cterm mods
        self.assertNotIn('internal', popped_mods)  # No internal mods

    def test_pop_mods_condense_false(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_mods(condense=False)
        
        # Should contain all keys, even empty ones
        all_mod_types = ['labile', 'static', 'unknown', 'nterm', 'cterm', 
                        'internal', 'interval', 'charge', 'charge_adducts', 'isotope']
        for mod_type in all_mod_types:
            self.assertIn(mod_type, popped_mods)

    def test_pop_mods_none_parameter(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods=None)
        
        # Should pop all mods
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_mods_empty_list(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods=[])
        
        # Should pop nothing
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")
        self.assertEqual(len(popped_mods), 0)

    """
    TESTS FOR: complex scenarios and edge cases
    """
    def test_pop_methods_preserve_sequence(self):
        original_sequence = "PEPTIDE"
        
        # Test each pop method preserves the sequence
        test_cases = [
            ("{Glycan}PEPTIDE", "pop_labile_mods"),
            ("[Unknown]?PEPTIDE", "pop_unknown_mods"),
            ("[Acetyl]-PEPTIDE", "pop_nterm_mods"),
            ("PEPTIDE-[Amide]", "pop_cterm_mods"),
            ("PE[Phospho]PTIDE", "pop_internal_mods"),
            ("P(EP)TIDE", "pop_intervals"),
            ("PEPTIDE/2", "pop_charge"),
            ("PEPTIDE/2[Na+]", "pop_charge_adducts"),
            ("<15N>PEPTIDE", "pop_isotope_mods"),
            ("<57@C>PEPTIDE", "pop_static_mods"),
        ]
        
        for proforma_string, method_name in test_cases:
            with self.subTest(method=method_name):
                annotation = pt.ProFormaAnnotation.parse(proforma_string)
                method = getattr(annotation, method_name)
                popped_data = method()
                self.assertEqual(annotation.sequence, original_sequence)

    def test_pop_methods_return_types(self):
        # Test return types are correct
        annotation = pt.ProFormaAnnotation.parse("{Glycan}<15N><57@C>[Unknown]?[Acetyl]-PE[Phospho]PTIDE-[Amide]/2[Na+]")
        
        # List return types
        list_methods = [
            'pop_labile_mods', 'pop_unknown_mods', 'pop_nterm_mods', 
            'pop_cterm_mods', 'pop_charge_adducts', 'pop_isotope_mods', 
            'pop_static_mods', 'pop_intervals'
        ]
        
        for method_name in list_methods:
            with self.subTest(method=method_name):
                annotation_copy = annotation.copy()
                method = getattr(annotation_copy, method_name)
                result = method()
                self.assertIsInstance(result, list)

        # Dict return type
        annotation_copy = annotation.copy()
        result = annotation_copy.pop_internal_mods()
        self.assertIsInstance(result, dict)

        # Optional int return type
        annotation_copy = annotation.copy()
        result = annotation_copy.pop_charge()
        self.assertIsInstance(result, int)

    def test_pop_mods_complex_annotation(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}<15N><57@C>[Unknown]?[Acetyl]-P(EP)[Phospho]TI[Methyl]DE-[Amide]/2[Na+]")
        popped_mods = annotation.pop_mods(mods=['internal', 'interval'])
        
        # Should only pop internal and interval modifications
        self.assertIn('internal', popped_mods)
        self.assertIn('interval', popped_mods)
        
        # Other modifications should still be present
        self.assertTrue(annotation.has_labile_mods())
        self.assertTrue(annotation.has_nterm_mods())
        self.assertTrue(annotation.has_cterm_mods())
        self.assertTrue(annotation.has_charge())

    def test_pop_mods_empty_annotation(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_mods()
        
        # Should return empty dict when condensed
        self.assertEqual(len(popped_mods), 0)

    def test_pop_mod_by_type_returns_empty_list_for_none(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        # Test that pop_mod_by_type returns empty list when no mods exist
        for mod_type in ['labile', 'unknown', 'nterm', 'cterm', 'isotope', 'static']:
            with self.subTest(mod_type=mod_type):
                result = annotation.pop_mod_by_type(mod_type)
                self.assertEqual(result, [])

    def test_multiple_pops_same_type(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        
        # First pop should return the mods
        first_pop = annotation.pop_nterm_mods()
        self.assertEqual(len(first_pop), 1)
        
        # Second pop should return empty list
        second_pop = annotation.pop_nterm_mods()
        self.assertEqual(len(second_pop), 0)

    def test_pop_mods_with_multipliers(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]^2PTIDE")
        popped_mods = annotation.pop_internal_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertIn(1, popped_mods)
        self.assertEqual(popped_mods[1][0].val, "Phospho")
        self.assertEqual(popped_mods[1][0].mult, 2)


if __name__ == '__main__':
    unittest.main()