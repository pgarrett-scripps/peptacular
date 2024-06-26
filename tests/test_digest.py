import unittest

import peptacular as pt

PROTEIN = 'MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTE'


class TestDigest(unittest.TestCase):

    def test_trypsin_sites(self):
        cleavage_sites = pt.get_cleavage_sites(PROTEIN, pt.PROTEASES['trypsin'])

        sites = [23, 31, 35, 42, 45, 47, 50, 56, 60, 62, 70, 73, 79, 88, 96]

        self.assertEqual(cleavage_sites, sites)

    def test_digest_protein(self):
        peptides = set(pt.digest(sequence='TIDERTIDEKTIDE',
                              enzyme_regex=pt.PROTEASES['trypsin'],
                              missed_cleavages=2,
                              min_len=0,
                              max_len=100,
                              semi=False))

        self.assertEqual(peptides,
                         {'TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE'})

        peptides = set(pt.digest(sequence='TIDERTIDEKTIDE',
                              enzyme_regex=pt.PROTEASES['trypsin'],
                              missed_cleavages=1,
                              min_len=0,
                              max_len=100,
                              semi=False))
        self.assertEqual(peptides, {'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDE', 'TIDE'})

        peptides = set(pt.digest(sequence='KTIDERTIDEKTIDE',
                              enzyme_regex=pt.PROTEASES['trypsin'],
                              missed_cleavages=1,
                              min_len=0,
                              max_len=100,
                              semi=False))
        self.assertEqual(peptides,
                         {'K', 'KTIDER', 'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDE', 'TIDE'})

        peptides = set(pt.digest(sequence='TIDERTIDEKTIDEK',
                              enzyme_regex=pt.PROTEASES['trypsin'],
                              missed_cleavages=1,
                              min_len=0,
                              max_len=100,
                              semi=False))
        self.assertEqual(peptides, {'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDEK', 'TIDEK'})

        peptides = set(pt.digest(sequence='TIDERTIDEKKTIDE',
                              enzyme_regex=pt.PROTEASES['trypsin'],
                              missed_cleavages=1,
                              min_len=0,
                              max_len=100,
                              semi=False))
        self.assertEqual(peptides,
                         {'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKK', 'K', 'KTIDE', 'TIDE'})

        peptides = set(pt.digest(sequence='TIDERTIDEKTIDE',
                              enzyme_regex=pt.PROTEASES['trypsin'],
                              missed_cleavages=0,
                              min_len=0,
                              max_len=100,
                              semi=False))
        self.assertEqual(peptides, {'TIDER', 'TIDEK', 'TIDE'})

        peptides = set(pt.digest(sequence='TIDERTIDEKTIDE',
                              enzyme_regex=pt.PROTEASES['trypsin'],
                              missed_cleavages=10,
                              min_len=0,
                              max_len=100,
                              semi=False))
        self.assertEqual(peptides,
                         {'TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE'})

    def test_get_semi_sequences(self):
        assert set(pt.get_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=3, max_len=5)) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP'}
        assert set(pt.get_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=2, max_len=5)) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE'}
        assert set(pt.get_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=2, max_len=100)) == \
               {'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'}
        assert set(pt.get_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=None, max_len=None)) == \
               {'P', 'E', 'TIDE', 'PTIDE', 'PEPT', 'PEPTI', 'IDE', 'PEP', 'PE', 'DE', 'PEPTID', 'EPTIDE'}

    def test_get_left_semi_sequences(self):
        self.assertEqual(set(pt.get_left_semi_enzymatic_sequences('PEPTIDE', None, None)),
                         {'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'})
        self.assertEqual(set(pt.get_left_semi_enzymatic_sequences('PEPTIDE', 3, None)), {'PEP', 'PEPT', 'PEPTI', 'PEPTID'})
        self.assertEqual(set(pt.get_left_semi_enzymatic_sequences('PEPTIDE', None, 5)), {'P', 'PE', 'PEP', 'PEPT', 'PEPTI'})

        assert set(pt.get_left_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=None, max_len=None)) == \
               {'P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID'}
        assert set(pt.get_left_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=3, max_len=None)) == \
               {'PEP', 'PEPT', 'PEPTI', 'PEPTID'}
        assert set(pt.get_left_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=None, max_len=5)) == \
               {'P', 'PE', 'PEP', 'PEPT', 'PEPTI'}
        assert set(pt.get_left_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=3, max_len=4)) == \
               {'PEP', 'PEPT'}

    def test_get_right_semi_sequences(self):
        self.assertEqual(set(pt.get_right_semi_enzymatic_sequences('PEPTIDE', None, None)),
                         {'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E'})
        self.assertEqual(set(pt.get_right_semi_enzymatic_sequences('PEPTIDE', 3, None)), {'EPTIDE', 'PTIDE', 'TIDE', 'IDE'})
        self.assertEqual(set(pt.get_right_semi_enzymatic_sequences('PEPTIDE', None, 5)), {'PTIDE', 'TIDE', 'IDE', 'DE', 'E'})
        assert set(pt.get_right_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=None, max_len=None)) == \
               {'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E'}
        assert set(pt.get_right_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=3, max_len=None)) == \
               {'EPTIDE', 'PTIDE', 'TIDE', 'IDE'}
        assert set(pt.get_right_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=None, max_len=5)) == \
               {'PTIDE', 'TIDE', 'IDE', 'DE', 'E'}
        assert set(pt.get_right_semi_enzymatic_sequences(sequence='PEPTIDE', min_len=3, max_len=4)) == \
               {'TIDE', 'IDE'}

    def test_get_all_non_enzymatic_peptides(self):
        self.assertEqual(set(pt.get_non_enzymatic_sequences('PEPT')),
                         {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT', 'PEP', 'EPT'})
        self.assertEqual(set(pt.get_non_enzymatic_sequences('PEPT', min_len=1, max_len=2)),
                         {'P', 'E', 'P', 'T', 'PE', 'EP', 'PT'})
        self.assertEqual(set(pt.get_non_enzymatic_sequences('PEPT', min_len=2, max_len=4)),
                         {'PT', 'EP', 'EPT', 'PE', 'PEP'})
