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



if __name__ == '__main__':
    unittest.main()
