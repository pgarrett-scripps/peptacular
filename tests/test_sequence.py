from peptacular.sequence import parse_modified_sequence, create_modified_sequence, strip_modifications, \
    get_semi_sequences, get_non_enzymatic_sequences, convert_to_ms2_pip_style, \
    get_left_semi_sequences, get_right_semi_sequences, add_static_mods, add_variable_mods, sequence_generator, \
    _sequence_generator_back, _sequence_generator_front, get_fragment_sequences, get_internal_fragment_sequences

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
        assert set(get_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=5)) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP'}
        assert set(get_semi_sequences(sequence='PEPTIDE', min_len=2, max_len=5)) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE'}
        assert set(get_semi_sequences(sequence='PEPTIDE', min_len=2, max_len=100)) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'}
        assert set(get_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None)) == \
               {'P', 'E', 'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'}

    def test_get_left_semi_sequences(self):
        self.assertEqual(set(get_left_semi_sequences('PEPTIDE', None, None)),
                         {'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'})
        self.assertEqual(set(get_left_semi_sequences('PEPTIDE', 3, None)), {'PEP', 'PEPT', 'PEPTI', 'PEPTID'})
        self.assertEqual(set(get_left_semi_sequences('PEPTIDE', None, 5)), {'P', 'PE', 'PEP', 'PEPT', 'PEPTI'})


        assert set(get_left_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None)) == \
               {'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'}
        assert set(get_left_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=None)) == \
               {'PEP', 'PEPT', 'PEPTI', 'PEPTID'}
        assert set(get_left_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=5)) == \
               {'P', 'PE', 'PEP', 'PEPT', 'PEPTI'}
        assert set(get_left_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=4)) == \
               {'PEP', 'PEPT'}

    def test_get_right_semi_sequences(self):
        self.assertEqual(set(get_right_semi_sequences('PEPTIDE', None, None)),
                         {'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E'})
        self.assertEqual(set(get_right_semi_sequences('PEPTIDE', 3, None)), {'EPTIDE', 'PTIDE', 'TIDE', 'IDE'})
        self.assertEqual(set(get_right_semi_sequences('PEPTIDE', None, 5)), {'PTIDE', 'TIDE', 'IDE', 'DE', 'E'})
        assert set(get_right_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None)) == \
               {'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E'}
        assert set(get_right_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=None)) == \
               {'EPTIDE', 'PTIDE', 'TIDE', 'IDE'}
        assert set(get_right_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=5)) == \
               {'PTIDE', 'TIDE', 'IDE', 'DE', 'E'}
        assert set(get_right_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=4)) == \
               {'TIDE', 'IDE'}

    def test_get_all_non_enzymatic_peptides(self):
        self.assertEqual(set(get_non_enzymatic_sequences('PEPT')),
                         {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'})
        self.assertEqual(set(get_non_enzymatic_sequences('PEPT', min_len=1, max_len=2)),
                         {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT'})
        self.assertEqual(set(get_non_enzymatic_sequences('PEPT', min_len=2, max_len=4)),
                         {'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'})

    def test_add_static_mod(self):
        self.assertEqual(add_static_mods('PEPCTIDE', {'C': 57.021464}), 'PEPC(57.021464)TIDE')
        self.assertEqual(add_static_mods('PEPC(57.021464)TIDE', {'C': 57.021464}), 'PEPC(57.021464)TIDE')
        self.assertEqual(add_static_mods('CPEPTIDEC', {'C': 57.021464}), 'C(57.021464)PEPTIDEC(57.021464)')
        self.assertEqual(add_static_mods('P(1)EPCTIDE(1)', {'C': 57.021464}), 'P(1)EPC(57.021464)TIDE(1)')

    def test_add_variable_mod(self):
        self.assertEqual(set(add_variable_mods('PEPCTIDCE', {'C': 57.021464}, 2)), {'PEPC(57.021464)TIDCE',
                                                                                    'PEPC(57.021464)TIDC(57.021464)E',
                                                                                    'PEPCTIDC(57.021464)E',
                                                                                    'PEPCTIDCE'})

        self.assertEqual(set(add_variable_mods('P(1)EPC(20)TIDCE', {'C': 57.021464}, 2)), {'P(1)EPC(20)TIDCE',
                                                                                        'P(1)EPC(20)TIDC(57.021464)E'})

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

    def test_get_fragment_sequences(self):
        sequence = "PEPTIDE"

        for ion_type in 'abc':
            fragments = get_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID', 'PEPTIDE'][::-1])

        for ion_type in 'xyz':
            fragments = get_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['E', 'DE', 'IDE', 'TIDE', 'PTIDE', 'EPTIDE', 'PEPTIDE'][::-1])

    def test_get_internal_fragment_sequences(self):
        sequence = "PEPTIDE"

        for ion_type in 'xyz':
            fragments = get_internal_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID', 'PEPTIDE'][::-1])

        for ion_type in 'abc':
            fragments = get_internal_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['E', 'DE', 'IDE', 'TIDE', 'PTIDE', 'EPTIDE', 'PEPTIDE'][::-1])

    def test_get_fragment_sequences_modified(self):
        sequence = "(-10)PEP(2)TIDE(100)"

        for ion_type in 'abc':
            fragments = get_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['(-10)P', '(-10)PE', '(-10)PEP(2)', '(-10)PEP(2)T', '(-10)PEP(2)TI', '(-10)PEP(2)TID', '(-10)PEP(2)TIDE(100)'][::-1])

        for ion_type in 'xyz':
            fragments = get_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['E(100)', 'DE(100)', 'IDE(100)', 'TIDE(100)', 'P(2)TIDE(100)', 'EP(2)TIDE(100)', '(-10)PEP(2)TIDE(100)'][::-1])

    def test_get_internal_fragment_sequences_modified(self):
        sequence = "(-10)PEP(2)TIDE(100)"

        for ion_type in 'xyz':
            fragments = get_internal_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['(-10)P', '(-10)PE', '(-10)PEP(2)', '(-10)PEP(2)T', '(-10)PEP(2)TI', '(-10)PEP(2)TID', '(-10)PEP(2)TIDE(100)'][::-1])

        for ion_type in 'abc':
            fragments = get_internal_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['E(100)', 'DE(100)', 'IDE(100)', 'TIDE(100)', 'P(2)TIDE(100)', 'EP(2)TIDE(100)', '(-10)PEP(2)TIDE(100)'][::-1])


if __name__ == "__main__":
    unittest.main()
