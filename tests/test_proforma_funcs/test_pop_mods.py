import unittest
import peptacular as pt


class TestPopMods(unittest.TestCase):
    
    def test_pop_labile_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        popped_mods = annotation.pop_labile_mods()

        self.assertIsNotNone(popped_mods)
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0], "Glycan")
        self.assertFalse(annotation.has_labile_mods)
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_labile_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_labile_mods()
        self.assertIsNone(popped_mods)

    def test_pop_labile_mods_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}{HexNAc}PEPTIDE")
        popped_mods = annotation.pop_labile_mods()

        self.assertEqual(len(popped_mods), 2)
        self.assertEqual([mod for mod in popped_mods], ["Glycan", "HexNAc"])
        self.assertFalse(annotation.has_labile_mods)

    def test_pop_unknown_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        popped_mods = annotation.pop_unknown_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0], "Unknown")
        self.assertFalse(annotation.has_unknown_mods)

    def test_pop_unknown_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_unknown_mods()
        self.assertIsNone(popped_mods)

    def test_pop_nterm_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_nterm_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0], "Acetyl")
        self.assertFalse(annotation.has_nterm_mods)

    def test_pop_nterm_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_nterm_mods()
        self.assertIsNone(popped_mods)

    def test_pop_nterm_mods_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl][Phospho]-PEPTIDE")
        popped_mods = annotation.pop_nterm_mods()
        
        self.assertEqual(len(popped_mods), 2)
        self.assertEqual([mod for mod in popped_mods], ["Acetyl", "Phospho"])

    def test_pop_cterm_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        popped_mods = annotation.pop_cterm_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0], "Amide")
        self.assertFalse(annotation.has_cterm_mods)

    def test_pop_cterm_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_cterm_mods()
        self.assertIsNone(popped_mods)

    def test_pop_internal_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        popped_mods = annotation.pop_internal_mods()
        
        self.assertEqual(len(popped_mods), 2)
        self.assertIn(1, popped_mods)
        self.assertIn(4, popped_mods)
        self.assertEqual(popped_mods[1][0], "Phospho")
        self.assertEqual(popped_mods[4][0], "Methyl")

    def test_pop_internal_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_internal_mods()
        self.assertIsNone(popped_mods)

    def test_pop_intervals_with_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        popped_intervals = annotation.pop_intervals()
        
        self.assertEqual(len(popped_intervals), 1)
        self.assertEqual(popped_intervals[0].start, 1)
        self.assertEqual(popped_intervals[0].end, 3)

    def test_pop_intervals_no_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_intervals = annotation.pop_intervals()
        self.assertIsNone(popped_intervals)

    def test_pop_charge_with_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        popped_charge = annotation.pop_charge()
        
        self.assertEqual(popped_charge, 2)
        self.assertFalse(annotation.has_charge)

    def test_pop_charge_no_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_charge = annotation.pop_charge()
        self.assertIsNone(popped_charge)

    def test_pop_charge_negative(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/-2")
        popped_charge = annotation.pop_charge()
        self.assertEqual(popped_charge, -2)

    def test_pop_charge_adducts_with_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2[Na+]")
        popped_adducts = annotation.pop_charge_adducts()
        
        self.assertEqual(len(popped_adducts), 1)
        self.assertEqual(popped_adducts[0], "Na+")

    def test_pop_charge_adducts_no_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        popped_adducts = annotation.pop_charge_adducts()
        self.assertIsNone(popped_adducts)

    def test_pop_charge_adducts_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2[Na+][K+]")
        popped_adducts = annotation.pop_charge_adducts()
        
        self.assertEqual(len(popped_adducts), 2)
        self.assertEqual([mod for mod in popped_adducts], ["Na+", "K+"])

    def test_pop_isotope_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        popped_mods = annotation.pop_isotope_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0], "15N")

    def test_pop_isotope_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_isotope_mods()
        self.assertIsNone(popped_mods)

    def test_pop_static_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        popped_mods = annotation.pop_static_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0], "57@C")

    def test_pop_static_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_static_mods()
        self.assertIsNone(popped_mods)

    def test_pop_mod_by_type_string(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_mod_by_type('nterm')
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0], "Acetyl")

    def test_pop_mod_by_type_enum(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_mod_by_type(pt.ModType.NTERM)
        
        self.assertEqual(len(popped_mods), 1)
        self.assertEqual(popped_mods[0], "Acetyl")

    def test_pop_mod_by_type_invalid(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        with self.assertRaises(ValueError):
            annotation.pop_mod_by_type('invalid_type')

    def test_pop_mods_all_default(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}<57@C>[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")
        popped_mods = annotation.pop_mods()
        
        expected_keys = ['labile', 'static', 'nterm', 'internal', 'cterm', 'charge']
        for key in expected_keys:
            self.assertIn(key, popped_mods)
        
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_mods_specific_type(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods='nterm')
        
        self.assertIn('nterm', popped_mods)
        self.assertEqual(len(popped_mods['nterm']), 1)
        self.assertTrue(annotation.has_internal_mods)
        self.assertTrue(annotation.has_cterm_mods)

    def test_pop_mods_multiple_types(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods=['nterm', 'cterm'])
        
        self.assertIn('nterm', popped_mods)
        self.assertIn('cterm', popped_mods)
        self.assertTrue(annotation.has_internal_mods)

    def test_pop_mods_enum_type(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods=pt.ModType.INTERNAL)
        
        self.assertIn('internal', popped_mods)
        self.assertTrue(annotation.has_nterm_mods)
        self.assertTrue(annotation.has_cterm_mods)

    def test_pop_mods_mixed_types(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")
        popped_mods = annotation.pop_mods(mods=['nterm', pt.ModType.CHARGE])
        
        self.assertIn('nterm', popped_mods)
        self.assertIn('charge', popped_mods)
        self.assertEqual(popped_mods['charge'], 2)

    def test_pop_mods_empty_annotation(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_mods()
        self.assertIsNone(popped_mods)

    def test_pop_mods_empty_list(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mods=[])
        
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")
        self.assertIsNone(popped_mods)

    def test_multiple_pops_return_none(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        
        first_pop = annotation.pop_nterm_mods()
        self.assertEqual(len(first_pop), 1)
        
        second_pop = annotation.pop_nterm_mods()
        self.assertIsNone(second_pop)

    def test_pop_mods_with_multipliers(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]^2PTIDE")
        popped_mods = annotation.pop_internal_mods()
        
        self.assertEqual(len(popped_mods), 1)
        self.assertIn(1, popped_mods)
        self.assertEqual(popped_mods[1][0], "Phospho")


if __name__ == '__main__':
    unittest.main()