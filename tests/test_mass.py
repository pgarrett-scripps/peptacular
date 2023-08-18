from peptacular.constants import MONO_ISOTOPIC_ATOMIC_MASSES, AVERAGE_ATOMIC_MASSES
from peptacular.mass import calculate_mz, calculate_mass

import unittest


class TestMass(unittest.TestCase):
    def test_calculate_mz_with_unmodified_peptide(self):
        sequence = 'PEPTIDE'
        places = 2

        self.assertAlmostEqual(799.359964, calculate_mz(sequence, charge=0, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(800.367241, calculate_mz(sequence, charge=1, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(400.687259, calculate_mz(sequence, charge=2, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(267.460598, calculate_mz(sequence, charge=3, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(200.847268, calculate_mz(sequence, charge=4, ion_type='y', monoisotopic=True, uwpr_mass=True), places)

        # Average mass is off by 0.004 Da
        self.assertAlmostEqual(799.822520, calculate_mz(sequence, charge=0, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(800.829796, calculate_mz(sequence, charge=1, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(400.918536, calculate_mz(sequence, charge=2, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(267.614783, calculate_mz(sequence, charge=3, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(200.962906, calculate_mz(sequence, charge=4, ion_type='y', monoisotopic=False, uwpr_mass=True), places)

    def test_calculate_mz_with_modified_peptide(self):
        sequence = '(15)P(-10)EPTIDE(100)'
        places = 2

        self.assertAlmostEqual(799.359964 + 105, calculate_mz(sequence, charge=0, ion_type='y', monoisotopic=True, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(800.367241 + 105 / 1, calculate_mz(sequence, charge=1, ion_type='y', monoisotopic=True, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(400.687259 + 105 / 2, calculate_mz(sequence, charge=2, ion_type='y', monoisotopic=True, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(267.460598 + 105 / 3, calculate_mz(sequence, charge=3, ion_type='y', monoisotopic=True, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(200.847268 + 105 / 4, calculate_mz(sequence, charge=4, ion_type='y', monoisotopic=True, uwpr_mass=True),
                               places)

        self.assertAlmostEqual(799.822520 + 105, calculate_mz(sequence, charge=0, ion_type='y', monoisotopic=False, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(800.829796 + 105 / 1, calculate_mz(sequence, charge=1, ion_type='y', monoisotopic=False, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(400.918536 + 105 / 2, calculate_mz(sequence, charge=2, ion_type='y', monoisotopic=False, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(267.614783 + 105 / 3, calculate_mz(sequence, charge=3, ion_type='y', monoisotopic=False, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(200.962906 + 105 / 4, calculate_mz(sequence, charge=4, ion_type='y', monoisotopic=False, uwpr_mass=True),
                               places)

    def test_calculate_mass_with_modified_peptide(self):
        sequence = 'PEPTIDE'
        places = 2

        self.assertAlmostEqual(799.359964, calculate_mass(sequence, charge=0, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(799.359964 + MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * 1,
                               calculate_mass(sequence, charge=1, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(799.359964 + MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * 2,
                               calculate_mass(sequence, charge=2, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(799.359964 + MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * 3,
                               calculate_mass(sequence, charge=3, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(799.359964 + MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * 4,
                               calculate_mass(sequence, charge=4, ion_type='y', monoisotopic=True, uwpr_mass=True), places)

        self.assertAlmostEqual(799.822520, calculate_mass(sequence, charge=0, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(799.822520 + AVERAGE_ATOMIC_MASSES['PROTON'] * 1,
                               calculate_mass(sequence, charge=1, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(799.822520 + AVERAGE_ATOMIC_MASSES['PROTON'] * 2,
                               calculate_mass(sequence, charge=2, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(799.822520 + AVERAGE_ATOMIC_MASSES['PROTON'] * 3,
                               calculate_mass(sequence, charge=3, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(799.822520 + AVERAGE_ATOMIC_MASSES['PROTON'] * 4,
                               calculate_mass(sequence, charge=4, ion_type='y', monoisotopic=False, uwpr_mass=True), places)

    def test_calculate_mass_with_unmodified_peptide(self):
        sequence = '(15)P(-10)EPTIDE(100)'
        places = 2

        self.assertAlmostEqual(799.359964 + 105, calculate_mass(sequence, charge=0, ion_type='y', monoisotopic=True, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(799.359964 + MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * 1 + 105,
                               calculate_mass(sequence, charge=1, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(799.359964 + MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * 2 + 105,
                               calculate_mass(sequence, charge=2, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(799.359964 + MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * 3 + 105,
                               calculate_mass(sequence, charge=3, ion_type='y', monoisotopic=True, uwpr_mass=True), places)
        self.assertAlmostEqual(799.359964 + MONO_ISOTOPIC_ATOMIC_MASSES['PROTON'] * 4 + 105,
                               calculate_mass(sequence, charge=4, ion_type='y', monoisotopic=True, uwpr_mass=True), places)

        self.assertAlmostEqual(799.822520 + 105, calculate_mass(sequence, charge=0, ion_type='y', monoisotopic=False, uwpr_mass=True),
                               places)
        self.assertAlmostEqual(799.822520 + AVERAGE_ATOMIC_MASSES['PROTON'] * 1 + 105,
                               calculate_mass(sequence, charge=1, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(799.822520 + AVERAGE_ATOMIC_MASSES['PROTON'] * 2 + 105,
                               calculate_mass(sequence, charge=2, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(799.822520 + AVERAGE_ATOMIC_MASSES['PROTON'] * 3 + 105,
                               calculate_mass(sequence, charge=3, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
        self.assertAlmostEqual(799.822520 + AVERAGE_ATOMIC_MASSES['PROTON'] * 4 + 105,
                               calculate_mass(sequence, charge=4, ion_type='y', monoisotopic=False, uwpr_mass=True), places)
