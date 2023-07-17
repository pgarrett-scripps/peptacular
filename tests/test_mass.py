import numpy as np

from peptacular.mass import fragment_sequence, calculate_mass_array

import unittest

class TestPeptide(unittest.TestCase):
    def test_peptide_fragments(self):
        pyteomics_y_frags = [346.1608765557, 173.584076511235, 249.10811270685, 125.05769458680999, 120.06551961887999,
                             60.536398042825]
        pyteomics_b_frags = [98.06004031562, 49.533658391195004, 227.10263340359, 114.05495493517999, 328.150311872,
                             164.578794169385]
        y_frags = list(fragment_sequence(sequence='PET', types=('y'), max_charge=2))
        b_frags = list(fragment_sequence(sequence='PET', types=('b'), max_charge=2))

        self.assertEqual(len(pyteomics_y_frags), len(y_frags))
        for f in pyteomics_y_frags:
            self.assertTrue(f in y_frags)

        self.assertEqual(len(pyteomics_b_frags), len(b_frags))
        for f in pyteomics_b_frags:
            self.assertTrue(f in b_frags)


    def test_peptide_fragments_mod(self):
        pyteomics_y_frags = [348.1608765557, 174.584076511235, 249.10811270685, 125.05769458680999, 120.06551961887999,
                             60.536398042825]
        pyteomics_b_frags = [100.06004031562, 50.533658391195004, 229.10263340359, 115.05495493517999, 330.150311872,
                             165.578794169385]

        y_frags = list(fragment_sequence(sequence='P(2)ET', types=('y'), max_charge=2))
        b_frags = list(fragment_sequence(sequence='P(2)ET', types=('b'), max_charge=2))

        self.assertEqual(len(pyteomics_y_frags), len(y_frags))
        for f in pyteomics_y_frags:
            self.assertTrue(f in y_frags)

        self.assertEqual(len(pyteomics_b_frags), len(b_frags))
        for f in pyteomics_b_frags:
            self.assertTrue(f in b_frags)


    def test_convert_to_mass_array(self):
        FRAGS = np.array([97.05276384885, 129.04259308796998, 97.05276384885, 101.04767846841,
                          113.08406397713001, 115.02694302383001, 129.04259308796998], dtype=np.float32)
        mass_arr = calculate_mass_array('PEPTIDE')
        for m in FRAGS:
            self.assertTrue(m in mass_arr)

        FRAGS = np.array([197.05276384885, 129.04259308796998, 197.05276384885, 101.04767846841,
                          113.08406397713001, 115.02694302383001, 229.04259308796998], dtype=np.float32)
        mass_arr = calculate_mass_array('(100)PEP(100)TIDE(100)')
        for m in FRAGS:
            self.assertTrue(m in mass_arr)

    def test_fragment2(self):
        pyteomics = {'a': [70.06512569606001, 199.10771878403003, 300.15539725243997],
                     'b': [98.06004031562, 227.10263340359, 328.150311872],
                     'c': [115.08658941662999, 244.1291825046, 345.17686097301],
                     'x': [372.14014111112, 275.08737726226997, 146.0447841743],
                     'y': [346.1608765557, 249.10811270685, 120.06551961887999],
                     'z': [329.13432745469, 232.08156360584, 103.03897051786998],
                     }

        for ion_type in 'abcxyz':
            frags = sorted(list(fragment_sequence(sequence='PET', types=(ion_type), max_charge=1)))
            pyteomics_frags = sorted(pyteomics[ion_type])
            for f, pf in zip(frags, pyteomics_frags):
                self.assertAlmostEqual(f, pf, 6)

if __name__ == "__main__":
    unittest.main()
