import unittest

from peptacular.fragment import Fragment
from peptacular.score import compute_fragment_matches, FragmentMatch, match_spectra


class TestScore(unittest.TestCase):

    def test_fragment_match(self):
        # Test basic functionality with 'th' tolerance
        fragments = [Fragment(sequence='PEPT', mz=100, charge=1, ion_type='y', number=1, internal=False, parent_number=4),
                     Fragment(sequence='PEPT', mz=200, charge=1, ion_type='b', number=1, internal=False, parent_number=4)]
        mz_spectrum = [99.9, 100.1, 200, 200.05, 300.1, 300.2]
        intensity_spectrum = [100, 200, 300, 400, 500, 600]

        result = compute_fragment_matches(fragments, mz_spectrum, intensity_spectrum, 0.2, 'th')
        assert result == [FragmentMatch(fragments[0], 99.9, 100), FragmentMatch(fragments[0], 100.1, 200),
                          FragmentMatch(fragments[1], 200, 300), FragmentMatch(fragments[1], 200.05, 400)]


    def test_match_spectra(self):
        # Test basic functionality with 'th' tolerance
        mz_spectrum1 = [100, 200, 300]
        mz_spectrum2 = [99.9, 100.1, 200, 200.05, 300.1, 300.2]
        result = match_spectra(mz_spectrum1, mz_spectrum2, 0.2, 'th')
        assert result == [(0, 2), (2, 4), (4, 6)]

        # Test basic functionality with 'th' tolerance
        mz_spectrum1 = [100, 200, 300]
        mz_spectrum2 = [99.9, 100.1, 200, 200.05, 300.1]
        result = match_spectra(mz_spectrum1, mz_spectrum2, 0.2, 'th')
        assert result == [(0, 2), (2, 4), (4, 5)]

        # Test with 'ppm' tolerance
        mz_spectrum1 = [100, 200, 300]
        mz_spectrum2 = [99.9, 100.1, 200, 200.05, 300.1, 300.2]
        result = match_spectra(mz_spectrum1, mz_spectrum2, 1000, 'ppm')  # 0.1 ppm tolerance
        assert result == [(0, 2), (2, 4), (4, 6)]

        # Test with no matches
        mz_spectrum1 = [50, 150, 250]
        mz_spectrum2 = [100, 200, 300]
        result = match_spectra(mz_spectrum1, mz_spectrum2, 0.05, 'th')
        assert result == [None, None, None]

        # Test with empty spectra
        mz_spectrum1 = []
        mz_spectrum2 = [100, 200, 300]
        result = match_spectra(mz_spectrum1, mz_spectrum2, 0.1, 'th')
        assert result == []

        mz_spectrum1 = [100, 200, 300]
        mz_spectrum2 = []
        result = match_spectra(mz_spectrum1, mz_spectrum2, 0.1, 'th')
        assert result == [None, None, None]
