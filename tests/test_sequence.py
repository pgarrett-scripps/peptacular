from peptacular.sequence import get_modifications, add_modifications, strip_modifications, \
    apply_static_modifications, apply_variable_modifications, pop_modifications, calculate_sequence_length, \
    shift_sequence, reverse_sequence, is_sequence_valid, span_to_sequence, split_sequence

import unittest


class TestSequence(unittest.TestCase):

    def test_read_modified_peptide(self):
        # Test the case where there are no modifications in the peptide sequence
        peptide_sequence = "ACDEFG"
        expected_output = {}
        self.assertEqual(get_modifications(peptide_sequence), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        peptide_sequence = "A(2)CDEFG"
        expected_output = {0: 2}
        self.assertEqual(get_modifications(peptide_sequence), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        peptide_sequence = "A(2)C(3)DE(1)FG"
        expected_output = {0: 2, 1: 3, 3: 1}
        self.assertEqual(get_modifications(peptide_sequence), expected_output)

    def test_create_modified_peptide(self):
        # Test the case where there are no modifications in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {}
        expected_output = "ACDEFG"
        self.assertEqual(add_modifications(unmodified_sequence, modifications), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {0: 2}
        expected_output = "A(2)CDEFG"
        self.assertEqual(add_modifications(unmodified_sequence, modifications), expected_output)

        # Test the case where there is a single N-Term mod
        unmodified_sequence = "ACDEFG"
        modifications = {-1: 2}
        expected_output = "[2]ACDEFG"
        self.assertEqual(add_modifications(unmodified_sequence, modifications), expected_output)

        # Test the case where there is a single C-Term mod
        unmodified_sequence = "ACDEFG"
        modifications = {5: 2}
        expected_output = "ACDEFG(2)"
        self.assertEqual(add_modifications(unmodified_sequence, modifications), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {0: 2, 1: 3, 3: 1}
        expected_output = "A(2)C(3)DE(1)FG"
        self.assertEqual(add_modifications(unmodified_sequence, modifications), expected_output)

        # Test the case where the index of a modification is invalid for the given peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {10: 2}
        with self.assertRaises(ValueError):
            add_modifications(unmodified_sequence, modifications)

    def test_get_unmodified_peptide(self):
        self.assertEqual(strip_modifications('ACDEFG'), 'ACDEFG')
        self.assertEqual(strip_modifications('A(2)CDE(2)FG'), 'ACDEFG')
        self.assertEqual(strip_modifications('A(2)CDE(Acetyl)FG'), 'ACDEFG')

    def test_add_static_mod(self):
        self.assertEqual(apply_static_modifications('PEPCTIDE', {'C': 57.021464}), 'PEPC(57.021464)TIDE')
        self.assertEqual(apply_static_modifications('PEPC(57.021464)TIDE', {'C': 57.021464}), 'PEPC(57.021464)TIDE')
        self.assertEqual(apply_static_modifications('CPEPTIDEC', {'C': 57.021464}), 'C(57.021464)PEPTIDEC(57.021464)')
        self.assertEqual(apply_static_modifications('P(1)EPCTIDE(1)', {'C': 57.021464}), 'P(1)EPC(57.021464)TIDE(1)')

    def test_add_variable_mod(self):
        self.assertEqual(set(apply_variable_modifications('PEPCTIDCE', {'C': 57.021464}, 2)), {'PEPC(57.021464)TIDCE',
                                                                                    'PEPC(57.021464)TIDC(57.021464)E',
                                                                                    'PEPCTIDC(57.021464)E',
                                                                                    'PEPCTIDCE'})

        self.assertEqual(set(apply_variable_modifications('P(1)EPC(20)TIDCE', {'C': 57.021464}, 2)), {'P(1)EPC(20)TIDCE',
                                                                                           'P(1)EPC(20)TIDC(57.021464)E'})

    def test_calculate_sequence_length(self):
        self.assertEqual(calculate_sequence_length("[Acetyl]PEP(1.2345)TID(3.14)E[Amide]"), 7)
        self.assertEqual(calculate_sequence_length("PEPTIDE"), 7)

    def test_get_modifications(self):
        self.assertEqual(get_modifications('PEP(Phospho)T(1)IDE(3.14)'), {2: 'Phospho', 3: 1, 6: 3.14})
        self.assertEqual(get_modifications('[Acetyl]PEPTIDE[Amide]'), {-1: 'Acetyl', 7: 'Amide'})
        self.assertEqual(get_modifications('PEPTIDE'), {})

    def test_add_modifications(self):
        self.assertEqual(add_modifications('PEPTIDE', {2: 'phospho'}), 'PEP(phospho)TIDE')
        self.assertEqual(add_modifications('PEPTIDE', {-1: 'Acetyl', 6: '1.234', 7: 'Amide'}),
                         '[Acetyl]PEPTIDE(1.234)[Amide]')
        self.assertEqual(add_modifications('PEPTIDE', {}), 'PEPTIDE')

    def test_pop_modifications(self):
        self.assertEqual(pop_modifications('PEP(phospho)TIDE'), ('PEPTIDE', {2: 'phospho'}))

    def test_strip_modifications(self):
        self.assertEqual(strip_modifications('PEP(phospho)TIDE'), 'PEPTIDE')

    def test_apply_static_modifications(self):
        self.assertEqual(apply_static_modifications('PEPTIDE', {'P': 'phospho'}), 'P(phospho)EP(phospho)TIDE')

    def test_apply_variable_modifications(self):
        self.assertEqual(apply_variable_modifications('PEPTIDE', {'P': 'phospho'}, 2),
                         ['P(phospho)EP(phospho)TIDE', 'P(phospho)EPTIDE', 'PEP(phospho)TIDE', 'PEPTIDE'])

    def test_reverse_sequence(self):
        self.assertEqual(reverse_sequence('PEPTIDE'), 'EDITPEP')
        self.assertEqual(reverse_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]'),
                         '[Acetyl]EDITP(phospho)EP(phospho)[Amide]')

    def test_shift_sequence(self):
        self.assertEqual(shift_sequence('PEPTIDE', 2), 'PTIDEPE')

    def test_is_sequence_valid(self):
        self.assertTrue(is_sequence_valid('PEPTIDE'))
        self.assertTrue(is_sequence_valid('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]'))
        self.assertFalse(is_sequence_valid('[Acetyl]P(phospho)EP((phospho))TIDE[Amide]'))

    def test_span_to_sequence(self):
        self.assertEqual(span_to_sequence('PEPTIDE', (0, 4, 0)), 'PEPT')
        self.assertEqual(span_to_sequence('[Acetyl]PEP(1.2345)TID(3.14)E[Amide]', (0, 4, 0)), '[Acetyl]PEP(1.2345)T')

    def test_split_sequence(self):
        self.assertEqual(split_sequence('PEPTIDE'), ['P', 'E', 'P', 'T', 'I', 'D', 'E'])
        self.assertEqual(split_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]'),
                         ['[Acetyl]P(phospho)', 'E', 'P(phospho)', 'T', 'I', 'D', 'E[Amide]'])


if __name__ == "__main__":
    unittest.main()
