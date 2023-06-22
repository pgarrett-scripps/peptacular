from peptacular import constants
from peptacular.protein import calculate_protein_coverage

import unittest

from peptacular.sequence import identify_cleavage_sites, digest_sequence

PROTEIN = 'MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTE' \



class TestProtein(unittest.TestCase):

    def test_trypsin_sites(self):
        cleavage_sites = identify_cleavage_sites(PROTEIN, constants.PROTEASES['trypsin'])

        sites = [23, 31, 35, 42, 45, 47, 50, 56, 60, 62, 70, 73, 79, 88, 96]

        self.assertEqual(cleavage_sites, sites)

    def test_digest_protein(self):
        peptides = set(digest_sequence(sequence='TIDERTIDEKTIDE',
                               enzyme_regexes=constants.PROTEASES['trypsin'],
                               missed_cleavages=2,
                               min_len=0,
                               max_len=100,
                               semi_enzymatic=False))
        self.assertEqual(peptides,
                         {'TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE','TIDE'})

        peptides = set(digest_sequence(sequence='TIDERTIDEKTIDE',
                               enzyme_regexes=constants.PROTEASES['trypsin'],
                               missed_cleavages=1,
                               min_len=0,
                               max_len=100,
                               semi_enzymatic=False))
        self.assertEqual(peptides, {'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDE', 'TIDE'})

        peptides = set(digest_sequence(sequence='KTIDERTIDEKTIDE',
                               enzyme_regexes=constants.PROTEASES['trypsin'],
                               missed_cleavages=1,
                               min_len=0,
                               max_len=100,
                               semi_enzymatic=False))
        self.assertEqual(peptides,
                         {'K', 'KTIDER', 'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDE', 'TIDE'})

        peptides = set(digest_sequence(sequence='TIDERTIDEKTIDEK',
                               enzyme_regexes=constants.PROTEASES['trypsin'],
                               missed_cleavages=1,
                               min_len=0,
                               max_len=100,
                               semi_enzymatic=False))
        self.assertEqual(peptides, {'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDEK', 'TIDEK'})

        peptides = set(digest_sequence(sequence='TIDERTIDEKKTIDE',
                               enzyme_regexes=constants.PROTEASES['trypsin'],
                               missed_cleavages=1,
                               min_len=0,
                               max_len=100,
                               semi_enzymatic=False))
        self.assertEqual(peptides,
                         {'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKK', 'K', 'KTIDE', 'TIDE'})

        peptides = set(digest_sequence(sequence='TIDERTIDEKTIDE',
                                       enzyme_regexes=constants.PROTEASES['trypsin'],
                                       missed_cleavages=0,
                                       min_len=0,
                                       max_len=100,
                                       semi_enzymatic=False))
        self.assertEqual(peptides, {'TIDER', 'TIDEK', 'TIDE'})

        peptides = set(digest_sequence(sequence='TIDERTIDEKTIDE',
                               enzyme_regexes=constants.PROTEASES['trypsin'],
                               missed_cleavages=10,
                               min_len=0,
                               max_len=100,
                               semi_enzymatic=False))
        self.assertEqual(peptides,
                         {'TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE'})

    def test_calculate_protein_coverage(self):
        protein_sequence = 'TIDERTIDEKTIDE'
        cov_arr = calculate_protein_coverage(protein=protein_sequence, peptides=['TIDER', 'TIDEK', 'TIDE'])
        self.assertEqual(sum(cov_arr), len(protein_sequence))
        self.assertEqual(cov_arr, [1] * len(protein_sequence))

        cov_arr = calculate_protein_coverage(protein=protein_sequence, peptides=['TIDER', 'TIDEK'])
        self.assertEqual(sum(cov_arr), len(protein_sequence) - 4)
        self.assertEqual(cov_arr, [1] * (len(protein_sequence) - 4) + [0] * 4)


if __name__ == "__main__":
    unittest.main()
