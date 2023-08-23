from peptacular.fragment import calculate_fragment_mz_values, sequence_trimmer, _trim_from_end, \
    _trim_from_start, create_fragment_sequences, create_fragment_internal_sequences

import unittest


class TestFragment(unittest.TestCase):
    def test_calculate_fragment_mz_values_with_unmodified_peptide(self):
        pyteomics_y_frags = [346.1608765557, 173.584076511235, 249.10811270685, 125.05769458680999, 120.06551961887999,
                             60.536398042825]
        pyteomics_b_frags = [98.06004031562, 49.533658391195004, 227.10263340359, 114.05495493517999, 328.150311872,
                             164.578794169385]
        y_frags = list(calculate_fragment_mz_values(sequence='PET', types=('y'), max_charge=2))
        b_frags = list(calculate_fragment_mz_values(sequence='PET', types=('b'), max_charge=2))

        self.assertEqual(len(pyteomics_y_frags), len(y_frags))
        for f in pyteomics_y_frags:
            self.assertTrue(f in y_frags)

        self.assertEqual(len(pyteomics_b_frags), len(b_frags))
        for f in pyteomics_b_frags:
            self.assertTrue(f in b_frags)

    def test_calculate_fragment_mz_values_with_modified_peptide(self):
        pyteomics_y_frags = [348.1608765557, 174.584076511235, 249.10811270685, 125.05769458680999, 120.06551961887999,
                             60.536398042825]
        pyteomics_b_frags = [100.06004031562, 50.533658391195004, 229.10263340359, 115.05495493517999, 330.150311872,
                             165.578794169385]

        y_frags = list(calculate_fragment_mz_values(sequence='P(2)ET', types=('y'), max_charge=2))
        b_frags = list(calculate_fragment_mz_values(sequence='P(2)ET', types=('b'), max_charge=2))

        self.assertEqual(len(pyteomics_y_frags), len(y_frags))
        for f in pyteomics_y_frags:
            self.assertTrue(f in y_frags)

        self.assertEqual(len(pyteomics_b_frags), len(b_frags))
        for f in pyteomics_b_frags:
            self.assertTrue(f in b_frags)

    def test_calculate_fragment_mz_values_with_all_ion_types(self):
        pyteomics = {'a': [70.06512569606001, 199.10771878403003, 300.15539725243997],
                     'b': [98.06004031562, 227.10263340359, 328.150311872],
                     'c': [115.08658941662999, 244.1291825046, 345.17686097301],
                     'x': [372.14014111112, 275.08737726226997, 146.0447841743],
                     'y': [346.1608765557, 249.10811270685, 120.06551961887999],
                     'z': [329.13432745469, 232.08156360584, 103.03897051786998],
                     }

        for ion_type in 'abcxyz':
            frags = sorted(list(calculate_fragment_mz_values(sequence='PET', types=(ion_type), max_charge=1)))
            pyteomics_frags = sorted(pyteomics[ion_type])
            for f, pf in zip(frags, pyteomics_frags):
                self.assertAlmostEqual(f, pf, 6)

    def test_sequence_trimmer(self):
        # Test forward generator
        generator = sequence_trimmer("A(P)C(P)G", True)
        self.assertEqual(next(generator), "A(P)C(P)G")
        self.assertEqual(next(generator), "C(P)G")
        self.assertEqual(next(generator), "G")
        with self.assertRaises(StopIteration):
            next(generator)

        # Test backward generator
        generator = sequence_trimmer("A(P)C(P)G", False)
        self.assertEqual(next(generator), "A(P)C(P)G")
        self.assertEqual(next(generator), "A(P)C(P)")
        self.assertEqual(next(generator), "A(P)")
        with self.assertRaises(StopIteration):
            next(generator)

    def test_trim_from_start(self):
        generator = _trim_from_end("A(P)C(P)G")
        self.assertEqual(next(generator), "A(P)C(P)G")
        self.assertEqual(next(generator), "A(P)C(P)")
        self.assertEqual(next(generator), "A(P)")
        with self.assertRaises(StopIteration):
            next(generator)

    def test_trim_from_end(self):
        generator = _trim_from_start("[P]A(P)C(P)G")
        self.assertEqual(next(generator), "[P]A(P)C(P)G")
        self.assertEqual(next(generator), "C(P)G")
        self.assertEqual(next(generator), "G")
        with self.assertRaises(StopIteration):
            next(generator)

    def test_create_fragment_sequences(self):
        sequence = "PEPTIDE"

        for ion_type in 'abc':
            fragments = create_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID', 'PEPTIDE'][::-1])

        for ion_type in 'xyz':
            fragments = create_fragment_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['E', 'DE', 'IDE', 'TIDE', 'PTIDE', 'EPTIDE', 'PEPTIDE'][::-1])

        modified_sequence = "[-10]PEP(2)TIDE(100)"

        for ion_type in 'abc':
            fragments = create_fragment_sequences(modified_sequence, ion_type)
            self.assertEqual(fragments,
                             ['[-10]P', '[-10]PE', '[-10]PEP(2)', '[-10]PEP(2)T', '[-10]PEP(2)TI', '[-10]PEP(2)TID',
                              '[-10]PEP(2)TIDE(100)'][::-1])

        for ion_type in 'xyz':
            fragments = create_fragment_sequences(modified_sequence, ion_type)
            self.assertEqual(fragments,
                             ['E(100)', 'DE(100)', 'IDE(100)', 'TIDE(100)', 'P(2)TIDE(100)', 'EP(2)TIDE(100)',
                              '[-10]PEP(2)TIDE(100)'][::-1])

    def test_create_internal_sequences(self):
        sequence = "PEPTIDE"

        for ion_type in 'xyz':
            fragments = create_fragment_internal_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['P', 'PE', 'PEP', 'PEPT', 'PEPTI', 'PEPTID', 'PEPTIDE'][::-1])

        for ion_type in 'abc':
            fragments = create_fragment_internal_sequences(sequence, ion_type)
            self.assertEqual(fragments, ['E', 'DE', 'IDE', 'TIDE', 'PTIDE', 'EPTIDE', 'PEPTIDE'][::-1])

        modified_sequence = "[-10]PEP(2)TIDE(100)"

        for ion_type in 'xyz':
            fragments = create_fragment_internal_sequences(modified_sequence, ion_type)
            self.assertEqual(fragments,
                             ['[-10]P', '[-10]PE', '[-10]PEP(2)', '[-10]PEP(2)T', '[-10]PEP(2)TI', '[-10]PEP(2)TID',
                              '[-10]PEP(2)TIDE(100)'][::-1])

        for ion_type in 'abc':
            fragments = create_fragment_internal_sequences(modified_sequence, ion_type)
            self.assertEqual(fragments,
                             ['E(100)', 'DE(100)', 'IDE(100)', 'TIDE(100)', 'P(2)TIDE(100)', 'EP(2)TIDE(100)',
                              '[-10]PEP(2)TIDE(100)'][::-1])


if __name__ == "__main__":
    unittest.main()
