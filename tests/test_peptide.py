from peptacular.peptide import parse_modified_peptide, create_modified_peptide, strip_modifications, \
    get_semi_sequences, get_non_enzymatic_sequences

import unittest


class TestPeptide(unittest.TestCase):

    def test_read_modified_peptide(self):
        # Test the case where there are no modifications in the peptide sequence
        peptide_sequence = "ACDEFG"
        expected_output = {}
        self.assertEqual(parse_modified_peptide(peptide_sequence), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        peptide_sequence = "A(2)CDEFG"
        expected_output = {0: 2}
        self.assertEqual(parse_modified_peptide(peptide_sequence), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        peptide_sequence = "A(2)C(3)DE(1)FG"
        expected_output = {0: 2, 1: 3, 3: 1}
        self.assertEqual(parse_modified_peptide(peptide_sequence), expected_output)

        # Test the case where the peptide sequence contains incorrect modification notation
        peptide_sequence = "A(2)CDEFG("
        with self.assertRaises(ValueError):
            parse_modified_peptide(peptide_sequence)

    def test_create_modified_peptide(self):
        # Test the case where there are no modifications in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {}
        expected_output = "ACDEFG"
        self.assertEqual(create_modified_peptide(unmodified_sequence, modifications), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {0: 2}
        expected_output = "A(2)CDEFG"
        self.assertEqual(create_modified_peptide(unmodified_sequence, modifications), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {0: 2, 1: 3, 3: 1}
        expected_output = "A(2)C(3)DE(1)FG"
        self.assertEqual(create_modified_peptide(unmodified_sequence, modifications), expected_output)

        # Test the case where the index of a modification is invalid for the given peptide sequence
        unmodified_sequence = "ACDEFG"
        modifications = {10: 2}
        with self.assertRaises(ValueError):
            create_modified_peptide(unmodified_sequence, modifications)

    def test_get_unmodified_peptide(self):
        self.assertEqual(strip_modifications('ACDEFG'), 'ACDEFG')
        self.assertEqual(strip_modifications('ACDEFG123'), 'ACDEFG')
        self.assertEqual(strip_modifications('ACDEFG123*'), 'ACDEFG')
        self.assertEqual(strip_modifications('&^*(%$'), '')
        self.assertEqual(strip_modifications('A1C2D3E4F5G6'), 'ACDEFG')
        self.assertEqual(strip_modifications('A-C-D-E-F-G'), 'ACDEFG')
        self.assertEqual(strip_modifications('A(2)CDEFG'), 'ACDEFG')

    def test_get_all_substrings(self):
        self.assertEqual(get_semi_sequences('PEPTIDE', 3, 5), {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP'})
        self.assertEqual(get_semi_sequences('PEPTIDE', 2, 5), {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE'})
        self.assertEqual(get_semi_sequences('PEPTIDE', 2, 100), {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'PEPTIDE', 'EPTIDE'})
        self.assertEqual(get_semi_sequences('PEPTIDE'), {'P', 'E', 'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'PEPTIDE', 'EPTIDE'})

    def test_get_all_non_enzymatic_peptides(self):
        self.assertEqual(get_non_enzymatic_sequences('PEPT'), {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'})
        self.assertEqual(get_non_enzymatic_sequences('PEPT', min_len=1, max_len=2), {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT'})
        self.assertEqual(get_non_enzymatic_sequences('PEPT', min_len=2, max_len=4), {'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'})

if __name__ == "__main__":
    unittest.main()
