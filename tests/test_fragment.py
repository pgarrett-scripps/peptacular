from peptacular.fragment import fragment, build_fragment_sequences, sequence_trimmer

import unittest


class TestFragment(unittest.TestCase):
    def test_calculate_fragment_mz_values_with_unmodified_peptide(self):
        pyteomics_y1_frags = [346.1608765557, 249.10811270685, 120.06551961887999]
        pyteomics_y2_frags = [173.584076511235, 125.05769458680999, 60.536398042825]
        pyteomics_b1_frags = [328.150311872, 227.10263340359, 98.06004031562]
        pyteomics_b2_frags = [164.578794169385, 114.05495493517999, 49.533658391195004]
        y1_frags = list(fragment(sequence='PET', ion_types='y', charges=1))
        y2_frags = list(fragment(sequence='PET', ion_types='y', charges=2))
        b1_frags = list(fragment(sequence='PET', ion_types='b', charges=1))
        b2_frags = list(fragment(sequence='PET', ion_types='b', charges=2))

        self.assertEqual(len(pyteomics_y1_frags), len(y1_frags))
        for f, pf in zip(y1_frags, pyteomics_y1_frags):
            self.assertAlmostEqual(f, pf, 6)

        self.assertEqual(len(pyteomics_y2_frags), len(y2_frags))
        for f, pf in zip(y2_frags, pyteomics_y2_frags):
            self.assertAlmostEqual(f, pf, 6)

        self.assertEqual(len(pyteomics_b1_frags), len(b1_frags))
        for f, pf in zip(b2_frags, pyteomics_b2_frags):
            self.assertAlmostEqual(f, pf, 6)

        self.assertEqual(len(pyteomics_b2_frags), len(b2_frags))
        for f, pf in zip(b2_frags, pyteomics_b2_frags):
            self.assertAlmostEqual(f, pf, 6)

    def test_calculate_fragment_mz_values_with_modified_peptide(self):
        pyteomics_y_frags = [174.584076511235, 125.05769458680999, 60.536398042825]
        pyteomics_b_frags = [165.578794169385, 115.05495493517999, 50.533658391195004]

        y_frags = list(fragment(sequence='P(2)ET', ion_types='y', charges=2))
        b_frags = list(fragment(sequence='P(2)ET', ion_types='b', charges=2))

        self.assertEqual(len(pyteomics_y_frags), len(y_frags))
        for f, pf in zip(y_frags, pyteomics_y_frags):
            self.assertAlmostEqual(f, pf, 6)

        self.assertEqual(len(pyteomics_b_frags), len(b_frags))
        for f, pf in zip(b_frags, pyteomics_b_frags):
            self.assertAlmostEqual(f, pf, 6)

    def test_calculate_fragment_mz_values_with_all_ion_types(self):
        pyteomics = {'a': [70.06512569606001, 199.10771878403003, 300.15539725243997],
                     'b': [98.06004031562, 227.10263340359, 328.150311872],
                     'c': [115.08658941662999, 244.1291825046, 345.17686097301],
                     'x': [372.14014111112, 275.08737726226997, 146.0447841743],
                     'y': [346.1608765557, 249.10811270685, 120.06551961887999],
                     'z': [329.13432745469, 232.08156360584, 103.03897051786998],
                     }

        for ion_type in 'abcxyz':
            frags = sorted(list(fragment(sequence='PET', ion_types=ion_type, charges=1)))
            pyteomics_frags = sorted(pyteomics[ion_type])
            for f, pf in zip(frags, pyteomics_frags):
                self.assertAlmostEqual(f, pf, 6)

    def test_create_fragment_sequences(self):
        sequence = "PEPTIDE"

        for ion_type in 'abc':
            fragments = list(build_fragment_sequences(sequence, ion_type))
            self.assertEqual(fragments, ['P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID', 'PEPTIDE'][::-1])

        for ion_type in 'xyz':
            fragments = list(build_fragment_sequences(sequence, ion_type))
            self.assertEqual(fragments, ['E', 'DE', 'IDE', 'TIDE', 'PTIDE', 'EPTIDE', 'PEPTIDE'][::-1])

        modified_sequence = "[-10]PEP(2)TIDE(100)"

        for ion_type in 'abc':
            fragments = list(build_fragment_sequences(modified_sequence, ion_type))
            self.assertEqual(fragments,
                             ['[-10]P', '[-10]PE', '[-10]PEP(2)', '[-10]PEP(2)T', '[-10]PEP(2)TI', '[-10]PEP(2)TID',
                              '[-10]PEP(2)TIDE(100)'][::-1])

        for ion_type in 'xyz':
            fragments = list(build_fragment_sequences(modified_sequence, ion_type))
            self.assertEqual(fragments,
                             ['E(100)', 'DE(100)', 'IDE(100)', 'TIDE(100)', 'P(2)TIDE(100)', 'EP(2)TIDE(100)',
                              '[-10]PEP(2)TIDE(100)'][::-1])

    def test_create_internal_sequences(self):
        sequence = "PEPTIDE"

        for ion_type in 'xyz':
            fragments = list(build_fragment_sequences(sequence, ion_type, True))
            self.assertEqual(fragments, ['P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID', 'PEPTIDE'][::-1])

        for ion_type in 'abc':
            fragments = list(build_fragment_sequences(sequence, ion_type, True))
            self.assertEqual(fragments, ['E', 'DE', 'IDE', 'TIDE', 'PTIDE', 'EPTIDE', 'PEPTIDE'][::-1])

        modified_sequence = "[-10]PEP(2)TIDE(100)"

        for ion_type in 'xyz':
            fragments = list(build_fragment_sequences(modified_sequence, ion_type, True))
            self.assertEqual(fragments,
                             ['[-10]P', '[-10]PE', '[-10]PEP(2)', '[-10]PEP(2)T', '[-10]PEP(2)TI', '[-10]PEP(2)TID',
                              '[-10]PEP(2)TIDE(100)'][::-1])

        for ion_type in 'abc':
            fragments = list(build_fragment_sequences(modified_sequence, ion_type, True))
            self.assertEqual(fragments,
                             ['E(100)', 'DE(100)', 'IDE(100)', 'TIDE(100)', 'P(2)TIDE(100)', 'EP(2)TIDE(100)',
                              '[-10]PEP(2)TIDE(100)'][::-1])


if __name__ == "__main__":
    unittest.main()
