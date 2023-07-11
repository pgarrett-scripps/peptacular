import numpy as np

from peptacular.sequence import parse_modified_sequence, create_modified_sequence, strip_modifications, \
    get_semi_sequences, get_non_enzymatic_sequences, convert_to_ms2_pip_style, convert_ip2_mod_to_uniprot_mod, \
    get_left_semi_sequences, get_right_semi_sequences, add_static_mods, add_variable_mods, sequence_generator, \
    _sequence_generator_back, _sequence_generator_front, fragment_sequence, calculate_mass_array

import unittest


class TestPeptide(unittest.TestCase):

    def test_read_modified_peptide(self):
        # Test the case where there are no modifications in the peptide sequence
        peptide_sequence = "ACDEFG"
        expected_output = {}
        self.assertEqual(parse_modified_sequence(peptide_sequence), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        peptide_sequence = "A(2)CDEFG"
        expected_output = {0: '2'}
        self.assertEqual(parse_modified_sequence(peptide_sequence), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        peptide_sequence = "A(2)C(3)DE(1)FG"
        expected_output = {0: '2', 1: '3', 3: '1'}
        self.assertEqual(parse_modified_sequence(peptide_sequence), expected_output)

        # Test the case where the peptide sequence contains incorrect modification notation
        peptide_sequence = "A(2)CDEFG("
        with self.assertRaises(ValueError):
            parse_modified_sequence(peptide_sequence)

    def test_convert_to_ms2_pip_style(self):
        # Test the case where there are no modifications in the peptide sequence
        peptide_sequence = "ACDEFG"
        expected_output = ''
        self.assertEqual(convert_to_ms2_pip_style(peptide_sequence), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        peptide_sequence = "A(Acetyl)CDEFG"
        expected_output = '1|Acetyl'
        self.assertEqual(convert_to_ms2_pip_style(peptide_sequence), expected_output)

        # Test the case where there is a single N-Term mod
        peptide_sequence = "(Acetyl)CDEFG"
        expected_output = '0|Acetyl'
        self.assertEqual(convert_to_ms2_pip_style(peptide_sequence), expected_output)

        # Test the case where there is a single C-Term mod
        peptide_sequence = "CDEFG(Acetyl)"
        expected_output = '-1|Acetyl'
        self.assertEqual(convert_to_ms2_pip_style(peptide_sequence), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        peptide_sequence = "A(2)C(Acetyl)DE(1)FG"
        expected_output = '1|2|2|Acetyl|4|1'
        self.assertEqual(convert_to_ms2_pip_style(peptide_sequence), expected_output)

        # Test the case where the peptide sequence contains incorrect modification notation
        peptide_sequence = "A(2)CDEFG("
        with self.assertRaises(ValueError):
            convert_to_ms2_pip_style(peptide_sequence)

    def test_convert_ip2_mod_to_uniprot_mod(self):
        # Test the case where there are multiple modifications in the peptide sequence
        peptide_sequence = "A(2)C(2)DE(2)FG"
        mod_map = {'2': 'Acetyl'}
        expected_output = "A(Acetyl)C(Acetyl)DE(Acetyl)FG"
        self.assertEqual(convert_ip2_mod_to_uniprot_mod(peptide_sequence, mod_map), expected_output)

    def test_create_modified_peptide(self):
        # Test the case where there are no modifications in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {}
        expected_output = "ACDEFG"
        self.assertEqual(create_modified_sequence(unmodified_sequence, modifications), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {0: 2}
        expected_output = "A(2)CDEFG"
        self.assertEqual(create_modified_sequence(unmodified_sequence, modifications), expected_output)

        # Test the case where there is a single N-Term mod
        unmodified_sequence = "ACDEFG"
        modifications = {-1: 2}
        expected_output = "(2)ACDEFG"
        self.assertEqual(create_modified_sequence(unmodified_sequence, modifications), expected_output)

        # Test the case where there is a single C-Term mod
        unmodified_sequence = "ACDEFG"
        modifications = {5: 2}
        expected_output = "ACDEFG(2)"
        self.assertEqual(create_modified_sequence(unmodified_sequence, modifications), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {0: 2, 1: 3, 3: 1}
        expected_output = "A(2)C(3)DE(1)FG"
        self.assertEqual(create_modified_sequence(unmodified_sequence, modifications), expected_output)

        # Test the case where the index of a modification is invalid for the given peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {10: 2}
        with self.assertRaises(ValueError):
            create_modified_sequence(unmodified_sequence, modifications)

    def test_get_unmodified_peptide(self):
        self.assertEqual(strip_modifications('ACDEFG'), 'ACDEFG')
        self.assertEqual(strip_modifications('A(2)CDE(2)FG'), 'ACDEFG')
        self.assertEqual(strip_modifications('A(2)CDE(Acetyl)FG'), 'ACDEFG')

    def test_get_semi_sequences(self):
        assert get_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=5) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP'}
        assert get_semi_sequences(sequence='PEPTIDE', min_len=2, max_len=5) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE'}
        assert get_semi_sequences(sequence='PEPTIDE', min_len=2, max_len=100) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'}
        assert get_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None) == \
               {'P', 'E', 'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'}

    def test_get_left_semi_sequences(self):
        self.assertEqual(get_left_semi_sequences('PEPTIDE', None, None, None),
                         {'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', 3, None, None), {'PEP', 'PEPT', 'PEPTI', 'PEPTID'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', None, 5, None), {'P', 'PE', 'PEP', 'PEPT', 'PEPTI'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', None, None, 1), {'PEPTID'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', None, None, 2), {'PEPTID', 'PEPTI'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', 3, 5, 1), {'PEPTI'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', 3, 5, 2), {'PEPTI', 'PEPT'})

        assert get_left_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None) == \
               {'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'}
        assert get_left_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=None) == \
               {'PEP', 'PEPT', 'PEPTI', 'PEPTID'}
        assert get_left_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=5) == \
               {'P', 'PE', 'PEP', 'PEPT', 'PEPTI'}
        assert get_left_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=4) == \
               {'PEP', 'PEPT'}

    def test_get_right_semi_sequences(self):
        self.assertEqual(get_right_semi_sequences('PEPTIDE', None, None, None),
                         {'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', 3, None, None), {'EPTIDE', 'PTIDE', 'TIDE', 'IDE'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', None, 5, None), {'PTIDE', 'TIDE', 'IDE', 'DE', 'E'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', None, None, 1), {'EPTIDE'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', None, None, 2), {'EPTIDE', 'PTIDE'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', 3, 5, 1), {'PTIDE'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', 3, 5, 2), {'PTIDE', 'TIDE'})

        assert get_right_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None) == \
               {'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E'}
        assert get_right_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=None) == \
               {'EPTIDE', 'PTIDE', 'TIDE', 'IDE'}
        assert get_right_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=5) == \
               {'PTIDE', 'TIDE', 'IDE', 'DE', 'E'}
        assert get_right_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=4) == \
               {'TIDE', 'IDE'}

    def test_get_all_non_enzymatic_peptides(self):
        self.assertEqual(get_non_enzymatic_sequences('PEPT'),
                         {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'})
        self.assertEqual(get_non_enzymatic_sequences('PEPT', min_len=1, max_len=2),
                         {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT'})
        self.assertEqual(get_non_enzymatic_sequences('PEPT', min_len=2, max_len=4),
                         {'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'})

    def test_add_static_mod(self):
        self.assertEqual(add_static_mods('PEPCTIDE', {'C': 57.021464}), 'PEPC(57.021464)TIDE')
        self.assertEqual(add_static_mods('PEPC(57.021464)TIDE', {'C': 57.021464}), 'PEPC(57.021464)TIDE')
        self.assertEqual(add_static_mods('CPEPTIDEC', {'C': 57.021464}), 'C(57.021464)PEPTIDEC(57.021464)')

    def test_add_variable_mod(self):
        self.assertEqual(set(add_variable_mods('PEPCTIDCE', {'C': 57.021464}, 2)), {'PEPC(57.021464)TIDCE',
                                                                                    'PEPC(57.021464)TIDC(57.021464)E',
                                                                                    'PEPCTIDC(57.021464)E',
                                                                                    'PEPCTIDCE'})

        self.assertEqual(set(add_variable_mods('PEPC(20)TIDCE', {'C': 57.021464}, 2)), {'PEPC(20)TIDCE',
                                                                                        'PEPC(20)TIDC(57.021464)E'})

    def test_peptide_generator(self):
        # Test forward generator
        generator = sequence_generator("A(P)C(P)G", True)
        self.assertEqual(next(generator), "A(P)C(P)G")
        self.assertEqual(next(generator), "C(P)G")
        self.assertEqual(next(generator), "G")
        with self.assertRaises(StopIteration):
            next(generator)

        # Test backward generator
        generator = sequence_generator("A(P)C(P)G", False)
        self.assertEqual(next(generator), "A(P)C(P)G")
        self.assertEqual(next(generator), "A(P)C(P)")
        self.assertEqual(next(generator), "A(P)")
        with self.assertRaises(StopIteration):
            next(generator)

    def test_peptide_generator_back(self):
        generator = _sequence_generator_back("A(P)C(P)G")
        self.assertEqual(next(generator), "A(P)C(P)G")
        self.assertEqual(next(generator), "A(P)C(P)")
        self.assertEqual(next(generator), "A(P)")
        with self.assertRaises(StopIteration):
            next(generator)

    def test_peptide_generator_front(self):
        generator = _sequence_generator_front("(P)A(P)C(P)G")
        self.assertEqual(next(generator), "(P)A(P)C(P)G")
        self.assertEqual(next(generator), "C(P)G")
        self.assertEqual(next(generator), "G")
        with self.assertRaises(StopIteration):
            next(generator)

    def test_peptide_fragments(self):
        pyteomics_y_frags = [346.1608765557, 173.584076511235, 249.10811270685, 125.05769458680999, 120.06551961887999,
                             60.536398042825]
        pyteomics_b_frags = [98.06004031562, 49.533658391195004, 227.10263340359, 114.05495493517999, 328.150311872,
                             164.578794169385]
        y_frags = list(fragment_sequence(sequence='PET', types=('y'), max_charge=2))
        b_frags = list(fragment_sequence(sequence='PET', types=('b'), max_charge=2))

        self.assertEqual(len(pyteomics_y_frags), len(y_frags))
        for f in pyteomics_y_frags:
            self.assertTrue(f in y_frags)

        self.assertEqual(len(pyteomics_b_frags), len(b_frags))
        for f in pyteomics_b_frags:
            self.assertTrue(f in b_frags)

    def test_peptide_fragments_mod(self):
        pyteomics_y_frags = [348.1608765557, 174.584076511235, 249.10811270685, 125.05769458680999, 120.06551961887999,
                             60.536398042825]
        pyteomics_b_frags = [100.06004031562, 50.533658391195004, 229.10263340359, 115.05495493517999, 330.150311872,
                             165.578794169385]

        y_frags = list(fragment_sequence(sequence='P(2)ET', types=('y'), max_charge=2))
        b_frags = list(fragment_sequence(sequence='P(2)ET', types=('b'), max_charge=2))

        self.assertEqual(len(pyteomics_y_frags), len(y_frags))
        for f in pyteomics_y_frags:
            self.assertTrue(f in y_frags)

        self.assertEqual(len(pyteomics_b_frags), len(b_frags))
        for f in pyteomics_b_frags:
            self.assertTrue(f in b_frags)

    def test_convert_to_mass_array(self):

        FRAGS = np.array([97.05276384885, 129.04259308796998, 97.05276384885, 101.04767846841,
                          113.08406397713001, 115.02694302383001, 129.04259308796998], dtype=np.float32)
        mass_arr = calculate_mass_array('PEPTIDE')
        for m in FRAGS:
            self.assertTrue(m in mass_arr)

        FRAGS = np.array([197.05276384885, 129.04259308796998, 197.05276384885, 101.04767846841,
                          113.08406397713001, 115.02694302383001, 229.04259308796998], dtype=np.float32)
        mass_arr = calculate_mass_array('(100)PEP(100)TIDE(100)')
        for m in FRAGS:
            self.assertTrue(m in mass_arr)

    def test_fragment2(self):
        pyteomics = {'a': [70.06512569606001, 199.10771878403003, 300.15539725243997],
                     'b': [98.06004031562, 227.10263340359, 328.150311872],
                     'c': [115.08658941662999, 244.1291825046, 345.17686097301],
                     'x': [372.14014111112, 275.08737726226997, 146.0447841743],
                     'y': [346.1608765557, 249.10811270685, 120.06551961887999],
                     'z': [329.13432745469, 232.08156360584, 103.03897051786998],
                     }

        for ion_type in 'abcxyz':
            frags = sorted(list(fragment_sequence(sequence='PET', types=(ion_type), max_charge=1)))
            pyteomics_frags = sorted(pyteomics[ion_type])
            for f, pf in zip(frags, pyteomics_frags):
                self.assertAlmostEqual(f, pf, 6)


if __name__ == "__main__":
    unittest.main()
