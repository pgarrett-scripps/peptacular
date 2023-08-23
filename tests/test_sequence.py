from peptacular.sequence import parse_modifications, add_modifications, strip_modifications, \
    apply_static_modifications, apply_variable_modifications

import unittest


class TestSequence(unittest.TestCase):

    def test_read_modified_peptide(self):
        # Test the case where there are no modifications in the peptide sequence
        peptide_sequence = "ACDEFG"
        expected_output = {}
        self.assertEqual(parse_modifications(peptide_sequence), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        peptide_sequence = "A(2)CDEFG"
        expected_output = {0: 2}
        self.assertEqual(parse_modifications(peptide_sequence), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        peptide_sequence = "A(2)C(3)DE(1)FG"
        expected_output = {0: 2, 1: 3, 3: 1}
        self.assertEqual(parse_modifications(peptide_sequence), expected_output)

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


if __name__ == "__main__":
    unittest.main()
