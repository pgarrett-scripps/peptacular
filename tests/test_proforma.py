import unittest

import peptacular as pt
from peptacular.proforma.proforma_dataclasses import Mod


class TestProForma(unittest.TestCase):

    def setUp(self):
        self.basic_mod = pt.Mod('HelloWorld', 1)

    def test_parse(self):
        # Test the case where there are no modifications in the peptide _sequence
        proforma = "PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE")
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_parse_empty(self):
        # Test the case where there are no modifications in the peptide _sequence
        proforma = ""
        expected_output = pt.ProFormaAnnotation(_sequence="")
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_single_mod(self):
        # Test the case where there is a single modification in the peptide _sequence
        proforma = "PEPTIDE[HelloWorld]"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={6: [self.basic_mod]})
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_single_double_mod(self):
        # Test the case where there is a single modification in the peptide _sequence
        proforma = "PEPTIDE[Hello][World]"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE",
                                             _internal_mods={6: [pt.Mod('Hello', 1), pt.Mod('World', 1)]})
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_single_double_multi_mod(self):
        # Test the case where there is a single modification in the peptide _sequence
        proforma = "PEPTIDE[Hello]^2[World]^3"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE",
                                             _internal_mods={6: [pt.Mod('Hello', 2), pt.Mod('World', 3)]})
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_nterm_mod(self):
        # Test the case where there is a single N-Term mod
        proforma = "[HelloWorld]-PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _nterm_mods=[self.basic_mod])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_nterm_double_mod(self):
        # Test the case where there is a single N-Term mod
        proforma = "[Hello][World]-PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _nterm_mods=[pt.Mod('Hello', 1), pt.Mod('World', 1)])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_nterm_double_multi_mod(self):
        # Test the case where there is a single N-Term mod
        proforma = "[Hello]^2[World]^3-PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _nterm_mods=[pt.Mod('Hello', 2), pt.Mod('World', 3)])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_cterm_mod(self):
        # Test the case where there is a single C-Term mod
        proforma = "PEPTIDE-[HelloWorld]"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _cterm_mods=[self.basic_mod])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_cterm_double_mod(self):
        # Test the case where there is a single C-Term mod
        proforma = "PEPTIDE-[Hello][World]"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _cterm_mods=[pt.Mod('Hello', 1), pt.Mod('World', 1)])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_cterm_double_multi_mod(self):
        # Test the case where there is a single C-Term mod
        proforma = "PEPTIDE-[Hello]^2[World]^3"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _cterm_mods=[pt.Mod('Hello', 2), pt.Mod('World', 3)])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_multiple_mods(self):
        # Test the case with multiple modifications at different positions
        proforma = "P[HelloWorld]EPT[HelloWorld]IDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE",
                                             _internal_mods={0: [self.basic_mod], 3: [self.basic_mod]})
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_interval_mods(self):
        # Test the case where there is an interval with modifications
        proforma = "PEP(TI)[HelloWorld]DE"
        expected_interval = pt.Interval(start=3, end=5, ambiguous=False, mods=[self.basic_mod])
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _intervals=[expected_interval])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_interval_mods_double(self):
        proforma = "PEP(TI)[Hello][World]DE"
        expected_interval = pt.Interval(start=3, end=5, ambiguous=False, mods=[pt.Mod('Hello', 1), pt.Mod('World', 1)])
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _intervals=[expected_interval])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_ambiguous_interval(self):
        # Test an ambiguous interval
        proforma = "PEP(?TI)DE"
        expected_interval = pt.Interval(start=3, end=5, ambiguous=True, mods=None)
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _intervals=[expected_interval])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_ambiguous_interval_with_mod(self):
        # Test an ambiguous interval with a modification
        proforma = "PEP(?TI)[HelloWorld]DE"
        expected_interval = pt.Interval(start=3, end=5, ambiguous=True, mods=[self.basic_mod])
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _intervals=[expected_interval])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_charge_state(self):
        # Test the case with a charge state at the end of the _sequence
        proforma = "PEPTIDE/2"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _charge=2)
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_negative_charge_state(self):
        # Test the case with a negative charge state at the end of the _sequence
        proforma = "PEPTIDE/-2"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _charge=-2)
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_charge_state_with_mod(self):
        # Test the case with a charge state at the end of the _sequence
        proforma = "PEPTIDE[HelloWorld]/2"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={6: [self.basic_mod]}, _charge=2)
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_charge_state_with_adduct(self):
        # Test the case with a charge state at the end of the _sequence
        proforma = "PEPTIDE/2[+H]"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _charge_adducts=[pt.Mod('+H', 1)], _charge=2)
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_labile_mod(self):
        # Test the case where there is a labile modification
        proforma = "{HelloWorld}PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _labile_mods=[self.basic_mod])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_labile_mod_double(self):
        # Test the case where there is a labile modification
        proforma = "{Hello}{World}PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _labile_mods=[pt.Mod('Hello', 1), pt.Mod('World', 1)])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_labile_mod_double_multi(self):
        # Test the case where there is a labile modification
        proforma = "{Hello}^2{World}^3PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _labile_mods=[pt.Mod('Hello', 2), pt.Mod('World', 3)])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_unknown_mod(self):
        # Test the case with an unknown modification
        proforma = "[HelloWorld]?PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _unknown_mods=[self.basic_mod])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_unknown_mod_double(self):
        # Test the case with an unknown modification
        proforma = "[Hello][World]?PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _unknown_mods=[pt.Mod('Hello', 1), pt.Mod('World', 1)])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_unknown_mod_double_multi(self):
        # Test the case with an unknown modification
        proforma = "[Hello]^2[World]^3?PEPTIDE"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _unknown_mods=[pt.Mod('Hello', 2), pt.Mod('World', 3)])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_isotope_mod(self):
        # Test the case with an isotope modification
        proforma = "<13C>PEPTIDE"
        expected_mod = pt.Mod('13C', 1)
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _isotope_mods=[expected_mod])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_isotope_mod_double(self):
        # Test the case with an isotope modification
        proforma = "<13C><15N>PEPTIDE"
        expected_mod1 = pt.Mod('13C', 1)
        expected_mod2 = pt.Mod('15N', 1)
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _isotope_mods=[expected_mod1, expected_mod2])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_isotope_mod_double_multi(self):
        # Test the case with an isotope modification
        proforma = "<13C>^2<15N>^3PEPTIDE"
        # expect ValueError
        self.assertRaises(pt.ProFormaFormatError, pt.parse, proforma)

    def test_static_mod_with_target(self):
        # Test static modification with a specified target
        proforma = "<Phospho@S>PEPTIDE"
        expected_mod = pt.Mod('Phospho@S', 1)
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _static_mods=[expected_mod])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_static_mod_with_target_double(self):
        # Test static modification with a specified target
        proforma = "<Phospho@S><Acetyl@K>PEPTIDE"
        expected_mod1 = pt.Mod('Phospho@S', 1)
        expected_mod2 = pt.Mod('Acetyl@K', 1)
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _static_mods=[expected_mod1, expected_mod2])
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_static_mod_with_target_double_multi(self):
        # Test static modification with a specified target
        proforma = "<Phospho@S>^2<Acetyl@K>^3PEPTIDE"
        self.assertRaises(ValueError, pt.parse, proforma)

    def test_multiple_types_of_mods(self):
        # Test a _sequence with multiple types of modifications
        nterm_mod = pt.Mod('Acetyl', 1)
        cterm_mod = pt.Mod('Amide', 1)
        internal_mod = pt.Mod('T[Phospho]', 1)
        proforma = "[Acetyl]-PEP[T[Phospho]]TIDE-[Amide]"
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _nterm_mods=[nterm_mod], _cterm_mods=[cterm_mod],
                                             _internal_mods={2: [internal_mod]})
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_mod_with_multiplier(self):
        # Test modification with a multiplier
        proforma = "PEPTIDE[HelloWorld]^3"
        expected_mod = pt.Mod('HelloWorld', 3)
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={6: [expected_mod]})
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_sequence_with_multiple_intervals(self):
        # Test a _sequence with multiple intervals
        interval_mod1 = pt.Mod('Mod1', 1)
        interval_mod2 = pt.Mod('Mod2', 1)
        proforma = "PEP(TI)[Mod1](DE)[Mod2]"
        expected_intervals = [pt.Interval(start=3, end=5, ambiguous=False, mods=[interval_mod1]),
                              pt.Interval(start=5, end=7, ambiguous=False, mods=[interval_mod2])]
        expected_output = pt.ProFormaAnnotation(_sequence="PEPTIDE", _intervals=expected_intervals)
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)

    def test_full_complex_sequence(self):
        # Test a complex _sequence with various modifications and features
        proforma = "{Oxidation}{Oxidation}<13C>[Acetyl]^2[ABC]-PEP[Phospho]T[Oxidation]^2IDE-[Amide]/3+PEPTIDE"
        expected_nterm_mod1 = pt.Mod('Acetyl', 2)
        expected_nterm_mod2 = pt.Mod('ABC', 1)
        expected_isotope_mod = pt.Mod('13C', 1)
        expected_internal_mod1 = pt.Mod('Phospho', 1)
        expected_internal_mod2 = pt.Mod('Oxidation', 2)
        expected_cterm_mod = pt.Mod('Amide', 1)
        expected_output = pt.MultiProFormaAnnotation(
            [
                pt.ProFormaAnnotation(
                    _sequence="PEPTIDE",
                    _labile_mods=[pt.Mod('Oxidation', 1), pt.Mod('Oxidation', 1)],
                    _nterm_mods=[expected_nterm_mod1, expected_nterm_mod2],
                    _isotope_mods=[expected_isotope_mod],
                    _internal_mods={2: [expected_internal_mod1], 3: [expected_internal_mod2]},
                    _cterm_mods=[expected_cterm_mod],
                    _charge=3),
                pt.ProFormaAnnotation(_sequence="PEPTIDE")
            ],
            connections=[False]
        )
        self.assertEqual(pt.parse(proforma), expected_output)
        self.assertEqual(pt.serialize(expected_output), proforma)


    def test_multi_annotation_crosslink(self):
        # Test creating and serializing crosslinked peptides
        pfa1 = pt.ProFormaAnnotation(_sequence='PEPTIDE')
        pfa2 = pt.ProFormaAnnotation(_sequence='ANOTHER')
        multi = pt.MultiProFormaAnnotation([pfa1, pfa2], [True])
        self.assertEqual(pt.serialize(multi), r'PEPTIDE\\ANOTHER')
        
    def test_multi_annotation_chain(self):
        # Test a chain of multiple peptides with different connection types
        pfa1 = pt.ProFormaAnnotation(_sequence='PEPTIDE')
        pfa2 = pt.ProFormaAnnotation(_sequence='ANOTHER')
        pfa3 = pt.ProFormaAnnotation(_sequence='THIRD')
        multi = pt.MultiProFormaAnnotation([pfa1, pfa2, pfa3], [True, False])
        self.assertEqual(pt.serialize(multi), r'PEPTIDE\\ANOTHER+THIRD')
        
    def test_multi_annotation_with_modifications(self):
        # Test with modifications on multiple peptides
        pfa1 = pt.ProFormaAnnotation(
            _sequence='PEPTIDE', 
            _internal_mods={2: [pt.Mod('Phospho', 1)]}
        )
        pfa2 = pt.ProFormaAnnotation(
            _sequence='ANOTHER', 
            _nterm_mods=[pt.Mod('Acetyl', 1)]
        )
        multi = pt.MultiProFormaAnnotation([pfa1, pfa2], [False])
        self.assertEqual(pt.parse(pt.serialize(multi)), multi)

    def test_slice_annotation(self):
        # Test slicing a ProForma annotation
        proforma = "P[Phospho]EPTIDE"
        annotation = pt.parse(proforma)
        sliced = annotation.slice(0, 4)
        self.assertEqual(sliced.sequence, "PEPT")
        self.assertEqual(sliced.internal_mods, {0: [Mod('Phospho', 1)]})
        
    def test_reverse_annotation(self):
        # Test reversing a ProForma annotation
        proforma = "P[Phospho]EPTIDE"
        annotation = pt.parse(proforma)
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.sequence, "EDITPEP")
        self.assertEqual(reversed_annotation.internal_mods, {6: [Mod('Phospho', 1)]})
        
    def test_shuffle_annotation(self):
        # Test shuffling a ProForma annotation with a fixed seed
        proforma = "P[Phospho]EPTIDE"
        annotation = pt.parse(proforma)
        shuffled = annotation.shuffle(seed=42)
        self.assertEqual(len(shuffled.sequence), len(annotation.sequence))
        self.assertEqual(shuffled.serialize(), "ETIPEP[Phospho]D")
        
    def test_shift_annotation(self):
        # Test shifting a ProForma annotation
        proforma = "P[Phospho]EPTIDE"
        annotation = pt.parse(proforma)
        shifted = annotation.shift(2)
        self.assertEqual(shifted.serialize(), "PTIDEP[Phospho]E")

    def test_is_subsequence(self):
        # Test subsequence detection
        parent = pt.parse("PEPTIDEPEPTIDE")
        child = pt.parse("TIDE")
        self.assertTrue(child.is_subsequence(parent))
        
    def test_find_indices(self):
        # Test finding indices of a subsequence
        parent = pt.parse("PEPTIDEPEPTIDE")
        child = pt.parse("PEP")
        indices = child.find_indices(parent)
        self.assertEqual(indices, [0, 7])
        
    def test_count_residues(self):
        # Test counting residues
        annotation = pt.parse("PEPTIDE")
        counts = annotation.count_residues()
        self.assertEqual(counts['P'], 2)
        self.assertEqual(counts['E'], 2)
        self.assertEqual(counts['T'], 1)
        self.assertEqual(counts['I'], 1)
        self.assertEqual(counts['D'], 1)
        
    def test_count_residues_with_modifications(self):
        # Test counting residues with modifications
        proforma = "P[Phospho]EPTIDE"
        annotation = pt.parse(proforma)
        counts = annotation.count_residues()
        self.assertEqual(counts['P'], 1)
        self.assertEqual(counts['P[Phospho]'], 1)
        self.assertEqual(counts['E'], 2)
        self.assertEqual(counts['T'], 1)
        self.assertEqual(counts['I'], 1)
        self.assertEqual(counts['D'], 1)