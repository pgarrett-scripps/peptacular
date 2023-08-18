from peptacular.protein import build_coverage_array, find_peptide_indexes, \
    calculate_percent_coverage

import unittest

class TestProtein(unittest.TestCase):

    def test_get_peptide_indexes_in_protein(self):
        # Test basic functionality
        self.assertEqual(find_peptide_indexes("AAPEPTIDEAA", "PEPTIDE"), [2])
        # Test with multiple occurrences
        self.assertEqual(find_peptide_indexes("AAPEPTIDEPEPTIDEAA", "PEPTIDE"), [2, 9])
        # Test peptide not in protein
        self.assertEqual(find_peptide_indexes("AAPEPTIDEAA", "XYZ"), [])

    def test_calculate_protein_coverage(self):
        # Test basic functionality
        protein = "AAPEPTIDEAA"
        peptides = ["PEPTIDE"]
        self.assertEqual(build_coverage_array(protein, peptides), [0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0])

        # Test with multiple peptides
        peptides = ["PEP", "TIDE"]
        self.assertEqual(build_coverage_array(protein, peptides), [0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0])

        # Test with peptide not in protein
        peptides = ["XYZ"]
        self.assertEqual(build_coverage_array(protein, peptides), [0] * 11)

    def test_calculate_protein_coverage_percent(self):
        # Test basic functionality
        protein = "AAPEPTIDEAA"
        peptides = ["PEPTIDE"]
        self.assertAlmostEqual(calculate_percent_coverage(protein, peptides), 7 / 11)

        # Test with multiple peptides
        peptides = ["PEP", "TIDE"]
        self.assertAlmostEqual(calculate_percent_coverage(protein, peptides), 7 / 11)

        # Test with peptide not in protein
        peptides = ["XYZ"]
        self.assertAlmostEqual(calculate_percent_coverage(protein, peptides), 0 / 11)

    def test_overlapping_peptides(self):
        protein = "AAPEPTIDEPEPTIDEAA"
        peptides = ["PEPTIDEPEP"]
        self.assertEqual(build_coverage_array(protein, peptides),
                         [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0])

    def test_empty_inputs(self):
        # Empty protein and peptides
        self.assertEqual(find_peptide_indexes("", ""), [])
        self.assertEqual(build_coverage_array("", []), [])
        self.assertEqual(calculate_percent_coverage("", []), 0)

        # Empty protein
        self.assertEqual(find_peptide_indexes("", "PEP"), [])
        self.assertEqual(build_coverage_array("", ["PEP"]), [])
        self.assertEqual(calculate_percent_coverage("", ["PEP"]), 0)

        # Empty peptide
        self.assertEqual(find_peptide_indexes("AAPEPTIDEAA", ""), [])
        self.assertEqual(build_coverage_array("AAPEPTIDEAA", [""]), [0] * 11)
        self.assertEqual(calculate_percent_coverage("AAPEPTIDEAA", [""]), 0)

    def test_peptide_longer_than_protein(self):
        # Peptide is longer than the protein
        self.assertEqual(find_peptide_indexes("PEP", "PEPTIDE"), [])
        self.assertEqual(build_coverage_array("PEP", ["PEPTIDE"]), [0, 0, 0])
        self.assertEqual(calculate_percent_coverage("PEP", ["PEPTIDE"]), 0)

    def test_full_coverage(self):
        # Peptides fully cover the protein
        peptides = ["AAPEP", "TIDEAA"]
        self.assertEqual(build_coverage_array("AAPEPTIDEAA", peptides), [1] * 11)
        self.assertEqual(calculate_percent_coverage("AAPEPTIDEAA", peptides), 1)

    def test_no_coverage(self):
        # No peptide matches in the protein
        peptides = ["XYZ", "LMN"]
        self.assertEqual(build_coverage_array("AAPEPTIDEAA", peptides), [0] * 11)
        self.assertEqual(calculate_percent_coverage("AAPEPTIDEAA", peptides), 0)

    def test_case_sensitivity(self):
        # Check if the functions are case-sensitive
        self.assertEqual(find_peptide_indexes("AAPEPTIDEAA", "peptide"), [])
        peptides = ["pep", "TIDE"]
        self.assertNotEqual(build_coverage_array("AAPEPTIDEAA", peptides), [1] * 11)
        self.assertNotEqual(calculate_percent_coverage("AAPEPTIDEAA", peptides), 1)


if __name__ == "__main__":
    unittest.main()
