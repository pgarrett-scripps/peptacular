from peptacular.peptide import parse_modified_peptide, create_modified_peptide, strip_modifications, \
    get_semi_sequences, get_non_enzymatic_sequences, convert_to_ms2_pip_style, convert_ip2_mod_to_uniprot_mod, \
    get_left_semi_sequences, get_right_semi_sequences

import unittest


class TestPeptide(unittest.TestCase):

    def test_read_modified_peptide(self):
        # Test the case where there are no modifications in the peptide sequence
        peptide_sequence = "ACDEFG"
        expected_output = {}
        self.assertEqual(parse_modified_peptide(peptide_sequence), expected_output)

        # Test the case where there is a single modification in the peptide sequence
        peptide_sequence = "A(2)CDEFG"
        expected_output = {0: '2'}
        self.assertEqual(parse_modified_peptide(peptide_sequence), expected_output)

        # Test the case where there are multiple modifications in the peptide sequence
        peptide_sequence = "A(2)C(3)DE(1)FG"
        expected_output = {0: '2', 1: '3', 3: '1'}
        self.assertEqual(parse_modified_peptide(peptide_sequence), expected_output)

        # Test the case where the peptide sequence contains incorrect modification notation
        peptide_sequence = "A(2)CDEFG("
        with self.assertRaises(ValueError):
            parse_modified_peptide(peptide_sequence)

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
        self.assertEqual(strip_modifications('A(2)CDE(2)FG'), 'ACDEFG')
        self.assertEqual(strip_modifications('A(2)CDE(Acetyl)FG'), 'ACDEFG')

    def test_get_all_substrings(self):
        self.assertEqual(get_semi_sequences('PEPTIDE', 3, 5), {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP'})
        self.assertEqual(get_semi_sequences('PEPTIDE', 2, 5),
                         {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE'})
        self.assertEqual(get_semi_sequences('PEPTIDE', 2, 100),
                         {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'})
        self.assertEqual(get_semi_sequences('PEPTIDE'),
                         {'P', 'E', 'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'})

    def test_get_left_substrings(self):
        self.assertEqual(get_left_semi_sequences('PEPTIDE', None, None, None),
                         {'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', 3, None, None), {'PEP', 'PEPT', 'PEPTI', 'PEPTID'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', None, 5, None), {'P', 'PE', 'PEP', 'PEPT', 'PEPTI'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', None, None, 1), {'PEPTID'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', None, None, 2), {'PEPTID', 'PEPTI'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', 3, 5, 1), {'PEPTI'})
        self.assertEqual(get_left_semi_sequences('PEPTIDE', 3, 5, 2), {'PEPTI', 'PEPT'})

    def test_get_right_substrings(self):
        self.assertEqual(get_right_semi_sequences('PEPTIDE', None, None, None),
                         {'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', 3, None, None), {'EPTIDE', 'PTIDE', 'TIDE', 'IDE'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', None, 5, None), {'PTIDE', 'TIDE', 'IDE', 'DE', 'E'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', None, None, 1), {'EPTIDE'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', None, None, 2), {'EPTIDE', 'PTIDE'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', 3, 5, 1), {'PTIDE'})
        self.assertEqual(get_right_semi_sequences('PEPTIDE', 3, 5, 2), {'PTIDE', 'TIDE'})

    def test_get_all_substrings_with_max_semi(self):
        self.assertEqual(get_semi_sequences('PEPTIDE', None, None, 1), {'PEPTID', 'EPTIDE'})

    def test_get_all_non_enzymatic_peptides(self):
        self.assertEqual(get_non_enzymatic_sequences('PEPT'),
                         {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'})
        self.assertEqual(get_non_enzymatic_sequences('PEPT', min_len=1, max_len=2),
                         {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT'})
        self.assertEqual(get_non_enzymatic_sequences('PEPT', min_len=2, max_len=4),
                         {'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'})

    def test_get_semi_digest(self):
        assert get_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=5) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP'}
        assert get_semi_sequences(sequence='PEPTIDE', min_len=2, max_len=5) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE'}
        assert get_semi_sequences(sequence='PEPTIDE', min_len=2, max_len=100) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'}
        assert get_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None) == \
               {'P', 'E', 'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'}

    def test_get_left_semi_sequences(self):
        assert get_left_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None) == \
               {'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'}
        assert get_left_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=None) == \
               {'PEP', 'PEPT', 'PEPTI', 'PEPTID'}
        assert get_left_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=5) == \
               {'P', 'PE', 'PEP', 'PEPT', 'PEPTI'}
        assert get_left_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=4) == \
               {'PEP', 'PEPT'}

    def test_get_right_semi_sequences(self):
        assert get_right_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=None) == \
               {'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E'}
        assert get_right_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=None) == \
               {'EPTIDE', 'PTIDE', 'TIDE', 'IDE'}
        assert get_right_semi_sequences(sequence='PEPTIDE', min_len=None, max_len=5) == \
               {'PTIDE', 'TIDE', 'IDE', 'DE', 'E'}
        assert get_right_semi_sequences(sequence='PEPTIDE', min_len=3, max_len=4) == \
               {'TIDE', 'IDE'}

    def test_get_all_non_enzymatic_sequences(self):
        assert get_non_enzymatic_sequences(sequence='PEPT', min_len=None, max_len=None) == \
               {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'}
        assert get_non_enzymatic_sequences(sequence='PEPT', min_len=1, max_len=2) == \
               {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT'}
        assert get_non_enzymatic_sequences(sequence='PEPT', min_len=2, max_len=4) == \
               {'PE', 'EP', 'PT', 'PEP', 'EPT', 'PEPT'}



if __name__ == "__main__":
    unittest.main()
