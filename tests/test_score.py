import unittest

import peptacular as pt


class TestScore(unittest.TestCase):

    def test_match_spectra(self):
        # Test basic functionality with 'th' tolerance
        mz_spectrum1 = [100, 200, 300]
        mz_spectrum2 = [99.9, 100.1, 200, 200.05, 300.1, 300.2]
        result = pt.get_matched_indices(mz_spectrum1, mz_spectrum2, 0.2, 'th')
        assert result == [(0, 2), (2, 4), (4, 6)]

        # Test basic functionality with 'th' tolerance
        mz_spectrum1 = [100, 200, 300]
        mz_spectrum2 = [99.9, 100.1, 200, 200.05, 300.1]
        result = pt.get_matched_indices(mz_spectrum1, mz_spectrum2, 0.2, 'th')
        assert result == [(0, 2), (2, 4), (4, 5)]

        # Test with 'ppm' tolerance
        mz_spectrum1 = [100, 200, 300]
        mz_spectrum2 = [99.9, 100.1, 200, 200.05, 300.1, 300.2]
        result = pt.get_matched_indices(mz_spectrum1, mz_spectrum2, 1000, 'ppm')  # 0.1 ppm tolerance
        assert result == [(0, 2), (2, 4), (4, 6)]

        # Test with no matches
        mz_spectrum1 = [50, 150, 250]
        mz_spectrum2 = [100, 200, 300]
        result = pt.get_matched_indices(mz_spectrum1, mz_spectrum2, 0.05, 'th')
        assert result == [None, None, None]

        # Test with empty spectra
        mz_spectrum1 = []
        mz_spectrum2 = [100, 200, 300]
        result = pt.get_matched_indices(mz_spectrum1, mz_spectrum2, 0.1, 'th')
        assert result == []

        mz_spectrum1 = [100, 200, 300]
        mz_spectrum2 = []
        result = pt.get_matched_indices(mz_spectrum1, mz_spectrum2, 0.1, 'th')
        assert result == [None, None, None]
