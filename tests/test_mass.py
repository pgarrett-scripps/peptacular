from peptacular.constants import ISOTOPIC_ATOMIC_MASSES, AVERAGE_ATOMIC_MASSES, PROTON_MASS
from peptacular.mass import mz, mass

import unittest


class TestMass(unittest.TestCase):
    def test_calculate_mz_with_unmodified_peptide(self):
        sequence = 'PEPTIDE'
        places = 2

        self.assertAlmostEqual(799.359964, mz(sequence, charge=0, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(800.367241, mz(sequence, charge=1, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(400.687259, mz(sequence, charge=2, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(267.460598, mz(sequence, charge=3, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(200.847268, mz(sequence, charge=4, ion_type='y', monoisotopic=True), places)

        # Average mass is off by 0.004 Da
        self.assertAlmostEqual(799.822520, mz(sequence, charge=0, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(800.829796, mz(sequence, charge=1, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(400.918536, mz(sequence, charge=2, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(267.614783, mz(sequence, charge=3, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(200.962906, mz(sequence, charge=4, ion_type='y', monoisotopic=False), places)

    def test_calculate_mz_with_modified_peptide(self):
        sequence = '(15)P(-10)EPTIDE(100)'
        places = 2

        self.assertAlmostEqual(799.359964 + 105, mz(sequence, charge=0, ion_type='y', monoisotopic=True),
                               places)
        self.assertAlmostEqual(800.367241 + 105 / 1, mz(sequence, charge=1, ion_type='y', monoisotopic=True),
                               places)
        self.assertAlmostEqual(400.687259 + 105 / 2, mz(sequence, charge=2, ion_type='y', monoisotopic=True),
                               places)
        self.assertAlmostEqual(267.460598 + 105 / 3, mz(sequence, charge=3, ion_type='y', monoisotopic=True),
                               places)
        self.assertAlmostEqual(200.847268 + 105 / 4, mz(sequence, charge=4, ion_type='y', monoisotopic=True),
                               places)

        self.assertAlmostEqual(799.822520 + 105, mz(sequence, charge=0, ion_type='y', monoisotopic=False),
                               places)
        self.assertAlmostEqual(800.829796 + 105 / 1, mz(sequence, charge=1, ion_type='y', monoisotopic=False),
                               places)
        self.assertAlmostEqual(400.918536 + 105 / 2, mz(sequence, charge=2, ion_type='y', monoisotopic=False),
                               places)
        self.assertAlmostEqual(267.614783 + 105 / 3, mz(sequence, charge=3, ion_type='y', monoisotopic=False),
                               places)
        self.assertAlmostEqual(200.962906 + 105 / 4, mz(sequence, charge=4, ion_type='y', monoisotopic=False),
                               places)

    def test_calculate_mass_with_modified_peptide(self):
        sequence = 'PEPTIDE'
        places = 2

        self.assertAlmostEqual(799.359964, mass(sequence, charge=0, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(799.359964 + PROTON_MASS * 1,
                               mass(sequence, charge=1, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(799.359964 + PROTON_MASS * 2,
                               mass(sequence, charge=2, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(799.359964 + PROTON_MASS * 3,
                               mass(sequence, charge=3, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(799.359964 + PROTON_MASS * 4,
                               mass(sequence, charge=4, ion_type='y', monoisotopic=True), places)

        self.assertAlmostEqual(799.822520, mass(sequence, charge=0, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(799.822520 + PROTON_MASS * 1,
                               mass(sequence, charge=1, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(799.822520 + PROTON_MASS * 2,
                               mass(sequence, charge=2, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(799.822520 + PROTON_MASS * 3,
                               mass(sequence, charge=3, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(799.822520 + PROTON_MASS * 4,
                               mass(sequence, charge=4, ion_type='y', monoisotopic=False), places)

    def test_calculate_mass_with_unmodified_peptide(self):
        sequence = '(15)P(-10)EPTIDE(100)'
        places = 2

        self.assertAlmostEqual(799.359964 + 105, mass(sequence, charge=0, ion_type='y', monoisotopic=True),
                               places)
        self.assertAlmostEqual(799.359964 + PROTON_MASS * 1 + 105,
                               mass(sequence, charge=1, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(799.359964 + PROTON_MASS * 2 + 105,
                               mass(sequence, charge=2, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(799.359964 + PROTON_MASS * 3 + 105,
                               mass(sequence, charge=3, ion_type='y', monoisotopic=True), places)
        self.assertAlmostEqual(799.359964 + PROTON_MASS * 4 + 105,
                               mass(sequence, charge=4, ion_type='y', monoisotopic=True), places)

        self.assertAlmostEqual(799.822520 + 105, mass(sequence, charge=0, ion_type='y', monoisotopic=False),
                               places)
        self.assertAlmostEqual(799.822520 + PROTON_MASS * 1 + 105,
                               mass(sequence, charge=1, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(799.822520 + PROTON_MASS * 2 + 105,
                               mass(sequence, charge=2, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(799.822520 + PROTON_MASS * 3 + 105,
                               mass(sequence, charge=3, ion_type='y', monoisotopic=False), places)
        self.assertAlmostEqual(799.822520 + PROTON_MASS * 4 + 105,
                               mass(sequence, charge=4, ion_type='y', monoisotopic=False), places)
