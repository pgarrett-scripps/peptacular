import unittest

from peptacular.sequence.proforma import parse, ProFormaAnnotation, Mod, Interval, serialize


class TestProForma(unittest.TestCase):

    def setUp(self):
        self.basic_mod = Mod('HelloWorld', 1)

    def test_parse(self):
        # Test the case where there are no modifications in the peptide _sequence
        proforma = "PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE")]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_parse_empty(self):
        # Test the case where there are no modifications in the peptide _sequence
        proforma = ""
        expected_output = [ProFormaAnnotation(_sequence="")]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_single_mod(self):
        # Test the case where there is a single modification in the peptide _sequence
        proforma = "PEPTIDE[HelloWorld]"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={6: [self.basic_mod]})]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_single_double_mod(self):
        # Test the case where there is a single modification in the peptide _sequence
        proforma = "PEPTIDE[Hello][World]"
        expected_output = [
            ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={6: [Mod('Hello', 1), Mod('World', 1)]})]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_single_double_multi_mod(self):
        # Test the case where there is a single modification in the peptide _sequence
        proforma = "PEPTIDE[Hello]^2[World]^3"
        expected_output = [
            ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={6: [Mod('Hello', 2), Mod('World', 3)]})]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_nterm_mod(self):
        # Test the case where there is a single N-Term mod
        proforma = "[HelloWorld]-PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _nterm_mods=[self.basic_mod])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_nterm_double_mod(self):
        # Test the case where there is a single N-Term mod
        proforma = "[Hello][World]-PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _nterm_mods=[Mod('Hello', 1), Mod('World', 1)])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_nterm_double_multi_mod(self):
        # Test the case where there is a single N-Term mod
        proforma = "[Hello]^2[World]^3-PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _nterm_mods=[Mod('Hello', 2), Mod('World', 3)])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_cterm_mod(self):
        # Test the case where there is a single C-Term mod
        proforma = "PEPTIDE-[HelloWorld]"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _cterm_mods=[self.basic_mod])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_cterm_double_mod(self):
        # Test the case where there is a single C-Term mod
        proforma = "PEPTIDE-[Hello][World]"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _cterm_mods=[Mod('Hello', 1), Mod('World', 1)])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_cterm_double_multi_mod(self):
        # Test the case where there is a single C-Term mod
        proforma = "PEPTIDE-[Hello]^2[World]^3"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _cterm_mods=[Mod('Hello', 2), Mod('World', 3)])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_multiple_mods(self):
        # Test the case with multiple modifications at different positions
        proforma = "P[HelloWorld]EPT[HelloWorld]IDE"
        expected_output = [
            ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={0: [self.basic_mod], 3: [self.basic_mod]})]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_interval_mods(self):
        # Test the case where there is an interval with modifications
        proforma = "PEP(TI)[HelloWorld]DE"
        expected_interval = Interval(start=3, end=5, ambiguous=False, mods=[self.basic_mod])
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _intervals=[expected_interval])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_interval_mods_double(self):
        proforma = "PEP(TI)[Hello][World]DE"
        expected_interval = Interval(start=3, end=5, ambiguous=False, mods=[Mod('Hello', 1), Mod('World', 1)])
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _intervals=[expected_interval])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_ambiguous_interval(self):
        # Test an ambiguous interval
        proforma = "PEP(?TI)DE"
        expected_interval = Interval(start=3, end=5, ambiguous=True, mods=None)
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _intervals=[expected_interval])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_ambiguous_interval_with_mod(self):
        # Test an ambiguous interval with a modification
        proforma = "PEP(?TI)[HelloWorld]DE"
        expected_interval = Interval(start=3, end=5, ambiguous=True, mods=[self.basic_mod])
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _intervals=[expected_interval])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_charge_state(self):
        # Test the case with a charge state at the end of the _sequence
        proforma = "PEPTIDE/2"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _charge=2)]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_negative_charge_state(self):
        # Test the case with a negative charge state at the end of the _sequence
        proforma = "PEPTIDE/-2"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _charge=-2)]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_charge_state_with_mod(self):
        # Test the case with a charge state at the end of the _sequence
        proforma = "PEPTIDE[HelloWorld]/2"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={6: [self.basic_mod]}, _charge=2)]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_charge_state_with_adduct(self):
        # Test the case with a charge state at the end of the _sequence
        proforma = "PEPTIDE/2[+H]"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _charge_adducts=[Mod('+H', 1)], _charge=2)]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_labile_mod(self):
        # Test the case where there is a labile modification
        proforma = "{HelloWorld}PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _labile_mods=[self.basic_mod])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_labile_mod_double(self):
        # Test the case where there is a labile modification
        proforma = "{Hello}{World}PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _labile_mods=[Mod('Hello', 1), Mod('World', 1)])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_labile_mod_double_multi(self):
        # Test the case where there is a labile modification
        proforma = "{Hello}^2{World}^3PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _labile_mods=[Mod('Hello', 2), Mod('World', 3)])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_unknown_mod(self):
        # Test the case with an unknown modification
        proforma = "[HelloWorld]?PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _unknown_mods=[self.basic_mod])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_unknown_mod_double(self):
        # Test the case with an unknown modification
        proforma = "[Hello][World]?PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _unknown_mods=[Mod('Hello', 1), Mod('World', 1)])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_unknown_mod_double_multi(self):
        # Test the case with an unknown modification
        proforma = "[Hello]^2[World]^3?PEPTIDE"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _unknown_mods=[Mod('Hello', 2), Mod('World', 3)])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_isotope_mod(self):
        # Test the case with an isotope modification
        proforma = "<13C>PEPTIDE"
        expected_mod = Mod('13C', 1)
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _isotope_mods=[expected_mod])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_isotope_mod_double(self):
        # Test the case with an isotope modification
        proforma = "<13C><15N>PEPTIDE"
        expected_mod1 = Mod('13C', 1)
        expected_mod2 = Mod('15N', 1)
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _isotope_mods=[expected_mod1, expected_mod2])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_isotope_mod_double_multi(self):
        # Test the case with an isotope modification
        proforma = "<13C>^2<15N>^3PEPTIDE"
        # expect ValueError
        self.assertRaises(ValueError, parse, proforma)

    def test_static_mod_with_target(self):
        # Test static modification with a specified target
        proforma = "<Phospho@S>PEPTIDE"
        expected_mod = Mod('Phospho@S', 1)
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _static_mods=[expected_mod])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_static_mod_with_target_double(self):
        # Test static modification with a specified target
        proforma = "<Phospho@S><Acetyl@K>PEPTIDE"
        expected_mod1 = Mod('Phospho@S', 1)
        expected_mod2 = Mod('Acetyl@K', 1)
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _static_mods=[expected_mod1, expected_mod2])]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_static_mod_with_target_double_multi(self):
        # Test static modification with a specified target
        proforma = "<Phospho@S>^2<Acetyl@K>^3PEPTIDE"
        self.assertRaises(ValueError, parse, proforma)

    def test_multiple_types_of_mods(self):
        # Test a _sequence with multiple types of modifications
        nterm_mod = Mod('Acetyl', 1)
        cterm_mod = Mod('Amide', 1)
        internal_mod = Mod('T[Phospho]', 1)
        proforma = "[Acetyl]-PEP[T[Phospho]]TIDE-[Amide]"
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _nterm_mods=[nterm_mod], _cterm_mods=[cterm_mod],
                                              _internal_mods={2: [internal_mod]})]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_mod_with_multiplier(self):
        # Test modification with a multiplier
        proforma = "PEPTIDE[HelloWorld]^3"
        expected_mod = Mod('HelloWorld', 3)
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _internal_mods={6: [expected_mod]})]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_sequence_with_multiple_intervals(self):
        # Test a _sequence with multiple intervals
        interval_mod1 = Mod('Mod1', 1)
        interval_mod2 = Mod('Mod2', 1)
        proforma = "PEP(TI)[Mod1](DE)[Mod2]"
        expected_intervals = [Interval(start=3, end=5, ambiguous=False, mods=[interval_mod1]),
                              Interval(start=5, end=7, ambiguous=False, mods=[interval_mod2])]
        expected_output = [ProFormaAnnotation(_sequence="PEPTIDE", _intervals=expected_intervals)]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)

    def test_full_complex_sequence(self):
        # Test a complex _sequence with various modifications and features
        proforma = "<13C>{Oxidation}{Oxidation}[Acetyl]^2[ABC]-PEP[Phospho]T[Oxidation]^2IDE-[Amide]/3+PEPTIDE"
        expected_nterm_mod1 = Mod('Acetyl', 2)
        expected_nterm_mod2 = Mod('ABC', 1)
        expected_isotope_mod = Mod('13C', 1)
        expected_internal_mod1 = Mod('Phospho', 1)
        expected_internal_mod2 = Mod('Oxidation', 2)
        expected_cterm_mod = Mod('Amide', 1)
        expected_output = [
            ProFormaAnnotation(
                _sequence="PEPTIDE",
                _labile_mods=[Mod('Oxidation', 1), Mod('Oxidation', 1)],
                _nterm_mods=[expected_nterm_mod1, expected_nterm_mod2],
                _isotope_mods=[expected_isotope_mod],
                _internal_mods={2: [expected_internal_mod1], 3: [expected_internal_mod2]},
                _cterm_mods=[expected_cterm_mod],
                _charge=3
            ),
            ProFormaAnnotation(_sequence="PEPTIDE")
        ]
        self.assertEqual(parse(proforma), expected_output)
        self.assertEqual(serialize(expected_output), proforma)
