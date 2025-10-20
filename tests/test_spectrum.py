"""
Simple tests for spectra compression and decompression to/from strings.
These tests verify the basic functionality: compress spectra to a string and decompress back.
"""
import unittest
import peptacular as pt



class TestSpectrumCreation(unittest.TestCase):
    """Simple tests for creating Spectrum objects."""
    
    def test_create_spectrum_from_tuple_of_lists(self):
        """Create Spectrum from tuple of (mzs, intensities)."""
        mzs = [100.0, 200.0, 300.0]
        intensities = [1000.0, 2000.0, 3000.0]
        
        spectrum = pt.Spectrum((mzs, intensities))
        
        self.assertEqual(len(spectrum), 3)
        self.assertEqual(spectrum.mzs, mzs)
        self.assertEqual(spectrum.intensities, intensities)
    
    def test_create_spectrum_from_tuple_with_charges(self):
        """Create Spectrum from tuple of (mzs, intensities, charges)."""
        mzs = [100.0, 200.0, 300.0]
        intensities = [1000.0, 2000.0, 3000.0]
        charges = [2, 3, None]
        
        spectrum = pt.Spectrum((mzs, intensities, charges))
        
        self.assertEqual(len(spectrum), 3)
        self.assertEqual(spectrum.mzs, mzs)
        self.assertEqual(spectrum.intensities, intensities)
        self.assertEqual(spectrum.charges, charges)
    
    def test_create_spectrum_from_peak_list(self):
        """Create Spectrum from list of Peak objects."""
        peaks = [
            pt.SpectrumPeak(100.0, 1000.0),
            pt.SpectrumPeak(200.0, 2000.0),
            pt.SpectrumPeak(300.0, 3000.0)
        ]
        
        spectrum = pt.Spectrum(peaks)
        
        self.assertEqual(len(spectrum), 3)
        self.assertEqual(spectrum.peaks, peaks)
    
    def test_create_spectrum_from_tuple_list(self):
        """Create Spectrum from list of (mz, intensity) tuples."""
        tuples = [(100.0, 1000.0), (200.0, 2000.0), (300.0, 3000.0)]
        
        spectrum = pt.Spectrum(tuples)
        
        self.assertEqual(len(spectrum), 3)
        self.assertEqual(spectrum.mzs, [100.0, 200.0, 300.0])
        self.assertEqual(spectrum.intensities, [1000.0, 2000.0, 3000.0])
    
    def test_create_spectrum_from_tuple_list_with_charges(self):
        """Create Spectrum from list of (mz, intensity, charge) tuples."""
        tuples = [(100.0, 1000.0, 2), (200.0, 2000.0, 3), (300.0, 3000.0, None)]
        
        spectrum = pt.Spectrum(tuples)
        
        self.assertEqual(len(spectrum), 3)
        self.assertEqual(spectrum.charges, [2, 3, None])
    
    def test_create_empty_spectrum(self):
        """Create empty Spectrum."""
        spectrum = pt.Spectrum([])
        
        self.assertEqual(len(spectrum), 0)
        self.assertFalse(spectrum)
        self.assertEqual(spectrum.mzs, [])
        self.assertEqual(spectrum.intensities, [])
    
    def test_spectrum_bool(self):
        """Test Spectrum bool conversion."""
        empty_spectrum = pt.Spectrum([])
        non_empty_spectrum = pt.Spectrum([(100.0, 1000.0)])
        
        self.assertFalse(empty_spectrum)
        self.assertTrue(non_empty_spectrum)


class TestSpectraStringCompression(unittest.TestCase):
    """Simple tests to verify spectra can be compressed to/from strings."""
    
    def test_compress_returns_string(self):
        """Verify compression returns a string."""
        mzs = [100.0, 200.0, 300.0]
        intensities = [1000.0, 2000.0, 3000.0]
        spectra = (mzs, intensities)
        
        compressed = pt.compress_spectra(spectra)
        
        # Check it's a string
        self.assertIsInstance(compressed, str)
        self.assertTrue(len(compressed) > 0)
    
    def test_decompress_from_string(self):
        """Verify decompression from string works."""
        mzs = [100.0, 200.0, 300.0]
        intensities = [1000.0, 2000.0, 3000.0]
        spectra = (mzs, intensities)
        
        # Compress to string
        compressed = pt.compress_spectra(spectra)
        
        # Decompress from string
        result = pt.decompress_spectra(compressed)
        
        # Check we got data back
        self.assertEqual(len(result), 2)
        decomp_mzs, decomp_intensities = result
        self.assertEqual(len(decomp_mzs), 3)
        self.assertEqual(len(decomp_intensities), 3)
    
    def test_round_trip_preserves_data(self):
        """Verify data is preserved through compress/decompress cycle."""
        mzs = [100.5, 200.75, 300.25]
        intensities = [1000.5, 2000.75, 3000.25]
        spectra = (mzs, intensities)
        
        # Compress and decompress
        compressed = pt.compress_spectra(spectra)
        decomp_mzs, decomp_intensities = pt.decompress_spectra(compressed)
        
        # Check values match (with float tolerance)
        for orig, decomp in zip(mzs, decomp_mzs):
            self.assertAlmostEqual(orig, decomp, places=4)
        
        for orig, decomp in zip(intensities, decomp_intensities):
            self.assertAlmostEqual(orig, decomp, places=2)
    
    def test_empty_spectra(self):
        """Verify empty spectra can be compressed/decompressed."""
        spectra = ([], [])
        
        compressed = pt.compress_spectra(spectra)
        decomp_mzs, decomp_intensities = pt.decompress_spectra(compressed)
        
        self.assertEqual(decomp_mzs, [])
        self.assertEqual(decomp_intensities, [])
    
    def test_url_safe_string(self):
        """Verify url_safe option produces URL-safe string."""
        mzs = [100.0, 200.0]
        intensities = [1000.0, 2000.0]
        spectra = (mzs, intensities)
        
        compressed = pt.compress_spectra(spectra, url_safe=True)
        
        # URL-safe strings start with 'U'
        self.assertTrue(compressed.startswith('U'))
        
        # Should still decompress correctly
        decomp_mzs, decomp_intensities = pt.decompress_spectra(compressed)
        self.assertEqual(len(decomp_mzs), 2)
        self.assertEqual(len(decomp_intensities), 2)


if __name__ == '__main__':
    unittest.main()
