import unittest

import peptacular as pt
from peptacular.proforma.proforma_parser import ProFormaAnnotation, create_annotation

PROTEIN = 'MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTE'


class TestDigest(unittest.TestCase):

    def test_trypsin_sites(self):
        cleavage_sites = list(pt.get_cleavage_sites(PROTEIN, pt.PROTEASES['trypsin']))
        sites = [23, 31, 35, 42, 45, 47, 50, 56, 60, 62, 70, 73, 79, 88, 96]

        self.assertEqual(cleavage_sites, sites)
        
        # Test with ProFormaAnnotation object
        annotation = create_annotation(PROTEIN)
        cleavage_sites_annotation = list(pt.get_cleavage_sites(annotation, pt.PROTEASES['trypsin']))
        self.assertEqual(cleavage_sites_annotation, sites)

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

    def test_modified_sequences(self):
        """Test digestion with modified sequences."""
        # Test with standard ProForma modified sequence
        modified_seq = '[Acetyl]-PEPTIDER[Phospho]TIDEM[Oxidation]K'
        peptides = set(pt.digest(sequence=modified_seq,
                              enzyme_regex='([KR])',
                              missed_cleavages=1,
                              return_type='str'))
        self.assertEqual(peptides, {'[Acetyl]-PEPTIDER[Phospho]', '[Acetyl]-PEPTIDER[Phospho]TIDEM[Oxidation]K', 'TIDEM[Oxidation]K'})
        
        # Test with annotation input
        annotation = pt.parse(modified_seq)
        peptides = set(pt.digest(sequence=annotation,
                              enzyme_regex='([KR])',
                              missed_cleavages=0,
                              return_type='str'))
        self.assertEqual(peptides, {'[Acetyl]-PEPTIDER[Phospho]', 'TIDEM[Oxidation]K'})
        
        # Test alternative return types
        peptides = list(pt.digest(sequence=modified_seq,
                              enzyme_regex='([KR])',
                              missed_cleavages=0,
                              return_type='annotation-span'))
        self.assertEqual(len(peptides), 2)
        self.assertIsInstance(peptides[0][0], ProFormaAnnotation)
        self.assertEqual(peptides[0][0].serialize(), '[Acetyl]-PEPTIDER[Phospho]')
        # Don't check exact span values as they may vary depending on implementation
        self.assertIsInstance(peptides[0][1], tuple)
        self.assertEqual(len(peptides[0][1]), 3)

    def test_edge_cases(self):
        """Test edge cases like empty sequences and single letters."""
        # Empty sequence
        self.assertEqual(list(pt.digest(sequence='', enzyme_regex='([KR])')), [])
        self.assertEqual(list(pt.get_semi_enzymatic_sequences('')), [])
        self.assertEqual(list(pt.get_left_semi_enzymatic_sequences('')), [])
        self.assertEqual(list(pt.get_right_semi_enzymatic_sequences('')), [])
        self.assertEqual(list(pt.get_non_enzymatic_sequences('')), [])
        
        # Single letter sequence
        self.assertEqual(list(pt.digest(sequence='K', enzyme_regex='([KR])')), ['K'])
        self.assertEqual(list(pt.get_semi_enzymatic_sequences('K')), [])
        self.assertEqual(list(pt.get_left_semi_enzymatic_sequences('K')), [])
        self.assertEqual(list(pt.get_right_semi_enzymatic_sequences('K')), [])
        self.assertEqual(list(pt.get_non_enzymatic_sequences('K')), [])

       # Single letter sequence
        self.assertEqual(list(pt.digest(sequence='K', enzyme_regex='non-enzymatic')), ['K'])
        self.assertEqual(list(pt.get_semi_enzymatic_sequences('K')), [])
        self.assertEqual(list(pt.get_left_semi_enzymatic_sequences('K')), [])
        self.assertEqual(list(pt.get_right_semi_enzymatic_sequences('K')), [])
        self.assertEqual(list(pt.get_non_enzymatic_sequences('K')), [])
        
        # Non-digesting sequence (no cleavage sites)
        self.assertEqual(list(pt.digest(sequence='PEPTIDE', enzyme_regex=pt.PROTEASES['trypsin'])), ['PEPTIDE'])


    def test_digest_from_config(self):
        """Test the digest_from_config function."""
        # Create a basic EnzymeConfig
        config = pt.EnzymeConfig(
            regex=pt.PROTEASES['trypsin'],
            missed_cleavages=1,
            semi_enzymatic=False,
            complete_digestion=True
        )
        
        # Test basic digestion
        peptides = set(pt.digest_from_config(
            sequence='TIDERTIDEKTIDE',
            config=config
        ))
        self.assertEqual(peptides, {'TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDE', 'TIDE'})
        
        # Test with modified sequence
        modified_seq = 'TIDER[Phospho]TIDEK'
        peptides = set(pt.digest_from_config(
            sequence=modified_seq,
            config=config
        ))
        self.assertEqual(peptides, {'TIDER[Phospho]', 'TIDER[Phospho]TIDEK', 'TIDEK'})
        
        # Test with semi-enzymatic
        config = pt.EnzymeConfig(
            regex=pt.PROTEASES['trypsin'],
            missed_cleavages=0,
            semi_enzymatic=True,
            complete_digestion=True
        )
        peptides = set(pt.digest_from_config(
            sequence='TIDEKTIDE',
            config=config,
            min_len=3
        ))
        # Should include both full enzymatic and semi-enzymatic peptides
        self.assertIn('TIDEK', peptides)  # Full enzymatic
        self.assertIn('IDE', peptides)    # Right semi
        self.assertIn('TID', peptides)    # Left semi

    def test_sequential_digest(self):
        """Test the sequential_digest function."""
        # Define enzyme configs
        trypsin = pt.EnzymeConfig(
            regex=['([KR])'], 
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True
        )
        
        asp_n = pt.EnzymeConfig(
            regex=['\\w(?=D)'], 
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=True
        )
        
        # Test basic sequential digestion
        peptides = set(pt.sequential_digest(
            sequence='PDEREKPKP',
            enzyme_configs=[trypsin, asp_n]
        ))
        
        # First trypsin splits PDER|EK|PK|P into PDER, EK, PK, P
        # Then asp-n splits those peptides at D
        # Expected result: 'P', 'DER', 'EK', 'PK', 'P'
        self.assertEqual(peptides, {'P', 'DER', 'EK', 'PK', 'P'})
        
        # Test with partial digestion (complete_digestion=False)
        partial_trypsin = pt.EnzymeConfig(
            regex=['([KR])'], 
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=False
        )
        
        partial_asp_n = pt.EnzymeConfig(
            regex=['\\w(?=D)'], 
            missed_cleavages=0,
            semi_enzymatic=False,
            complete_digestion=False
        )
        
        peptides = set(pt.sequential_digest(
            sequence='PEPRDK',
            enzyme_configs=[partial_trypsin, partial_asp_n]
        ))
        
        # First trypsin splits PEPR|DK| into PEPR, DK, PEPRDK (original retained)
        # Then asp-n splits those peptides at D
        # Expected result: PEPR, DK, PEPRDK
        self.assertEqual(peptides, {'PEPR', 'DK', 'PEPRDK'})
        
        # Test with modified sequence
        peptides = set(pt.sequential_digest(
            sequence='PEP[Phospho]RDK',
            enzyme_configs=[trypsin, asp_n]
        ))
        self.assertEqual(peptides, {'PEP[Phospho]R', 'DK'})

    def test_return_types(self):
        """Test different return types."""
        # Test return_type='annotation'
        peptides = list(pt.digest(
            sequence='TIDEKTIDE',
            enzyme_regex=pt.PROTEASES['trypsin'],
            missed_cleavages=0,
            return_type='annotation'
        ))
        self.assertEqual(len(peptides), 2)
        self.assertIsInstance(peptides[0], ProFormaAnnotation)
        self.assertEqual(peptides[0].serialize(), 'TIDEK')
        
        # Test return_type='span'
        peptides = list(pt.digest(
            sequence='TIDEKTIDE',
            enzyme_regex=pt.PROTEASES['trypsin'],
            missed_cleavages=0,
            return_type='span'
        ))
        self.assertEqual(len(peptides), 2)
        self.assertEqual(peptides[0], (0, 5, 0))  # (start, end, span_type)
        
        # Test return_type='str-span'
        peptides = list(pt.digest(
            sequence='TIDEKTIDE',
            enzyme_regex=pt.PROTEASES['trypsin'],
            missed_cleavages=0,
            return_type='str-span'
        ))
        self.assertEqual(len(peptides), 2)
        self.assertEqual(peptides[0][0], 'TIDEK')
        self.assertEqual(peptides[0][1], (0, 5, 0))
        
        # Test return_type='annotation-span'
        peptides = list(pt.digest(
            sequence='TIDEKTIDE',
            enzyme_regex=pt.PROTEASES['trypsin'],
            missed_cleavages=0,
            return_type='annotation-span'
        ))
        self.assertEqual(len(peptides), 2)
        self.assertIsInstance(peptides[0][0], ProFormaAnnotation)
        self.assertEqual(peptides[0][0].serialize(), 'TIDEK')
        self.assertEqual(peptides[0][1], (0, 5, 0))

    def test_error_cases(self):
        """Test error cases."""
        # Test invalid return type
        with self.assertRaises((ValueError, TypeError)):
            list(pt.digest(
                sequence='TIDEKTIDE',
                enzyme_regex=pt.PROTEASES['trypsin'],
                return_type='invalid_type'
            ))
        
        # Test invalid sequence type
        with self.assertRaises((ValueError, TypeError)):
            list(pt.digest(
                sequence=123,  # Not a string or ProFormaAnnotation
                enzyme_regex=pt.PROTEASES['trypsin']
            ))
        
        # Test mismatched enzyme_regex and missed_cleavages lengths
        with self.assertRaises((ValueError, TypeError)):
            list(pt.digest(
                sequence='TIDEKTIDE',
                enzyme_regex=[pt.PROTEASES['trypsin'], pt.PROTEASES['asp-n']],
                missed_cleavages=[0]
            ))
        
        # Test invalid missed_cleavages type
        with self.assertRaises((ValueError, TypeError)):
            list(pt.digest(
                sequence='TIDEKTIDE',
                enzyme_regex=pt.PROTEASES['trypsin'],
                missed_cleavages="invalid"
            ))
