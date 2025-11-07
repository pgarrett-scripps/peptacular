import unittest
import peptacular as pt

PRECISION_TOL = 4

class TestSpectraCompression(unittest.TestCase):
    
    def setUp(self):
        self.sample_mzs = [100.0, 150.5, 200.25, 300.75, 400.1]
        self.sample_intensities = [1000.0, 2500.5, 1200.25, 800.75, 3000.1]
        self.sample_spectra = (self.sample_mzs, self.sample_intensities)
    
    def test_basic_compression_round_trip(self):
        compressed = pt.compress_spectra(self.sample_spectra)
        decompressed = pt.decompress_spectra(compressed)
        
        mzs, intensities = decompressed
        self.assertEqual(len(mzs), len(self.sample_mzs))
        self.assertEqual(len(intensities), len(self.sample_intensities))
        
        for orig, decomp in zip(self.sample_mzs, mzs):
            self.assertAlmostEqual(orig, decomp, places=PRECISION_TOL)
        
        for orig, decomp in zip(self.sample_intensities, intensities):
            self.assertAlmostEqual(orig, decomp, places=2)
    
    def test_url_safe_compression(self):
        compressed = pt.compress_spectra(self.sample_spectra, url_safe=True)
        self.assertTrue(compressed.startswith('U'))
        
        decompressed = pt.decompress_spectra(compressed)
        mzs, intensities = decompressed
        
        self.assertEqual(len(mzs), len(self.sample_mzs))
        self.assertEqual(len(intensities), len(self.sample_intensities))
    
    def test_base85_compression(self):
        compressed = pt.compress_spectra(self.sample_spectra, url_safe=False)
        self.assertTrue(compressed.startswith('B'))
        
        decompressed = pt.decompress_spectra(compressed)
        mzs, intensities = decompressed
        
        self.assertEqual(len(mzs), len(self.sample_mzs))
        self.assertEqual(len(intensities), len(self.sample_intensities))
    
    def test_precision_compression(self):
        compressed = pt.compress_spectra(self.sample_spectra, intensity_precision=2)
        decompressed = pt.decompress_spectra(compressed)
        
        mzs, intensities = decompressed
        
        # MZ should be unchanged
        for orig, decomp in zip(self.sample_mzs, mzs):
            self.assertAlmostEqual(orig, decomp, places=PRECISION_TOL)
        
        # Intensities should be rounded to 2 decimal places
        for orig, decomp in zip(self.sample_intensities, intensities):
            self.assertEqual(round(orig, 2), round(decomp, 2))
    
    def test_empty_spectra(self):
        empty_spectra: tuple[list[float], list[float]] = ([], [])
        compressed = pt.compress_spectra(empty_spectra)
        decompressed = pt.decompress_spectra(compressed)
        
        self.assertEqual(decompressed, ([], []))
    
    def test_single_point_spectra(self):
        single_spectra = ([100.5], [1500.75])
        compressed = pt.compress_spectra(single_spectra)
        decompressed = pt.decompress_spectra(compressed)
        
        mzs, intensities = decompressed
        self.assertEqual(len(mzs), 1)
        self.assertEqual(len(intensities), 1)
        self.assertAlmostEqual(mzs[0], 100.5, places=PRECISION_TOL)
        self.assertAlmostEqual(intensities[0], 1500.75, places=PRECISION_TOL)
    
    def test_compressed_output_is_string(self):
        """Verify that compression returns a string."""
        compressed = pt.compress_spectra(self.sample_spectra)
        self.assertIsInstance(compressed, str)
        self.assertTrue(len(compressed) > 0)
    
    def test_compression_with_charges(self):
        """Test compression with charge information."""
        mzs = [100.0, 200.0, 300.0]
        intensities = [1000.0, 2000.0, 3000.0]
        charges = [2, 3, None]
        spectra_with_charges = (mzs, intensities, charges)
        
        compressed = pt.compress_spectra(spectra_with_charges)
        self.assertIsInstance(compressed, str)
        
        decompressed = pt.decompress_spectra(compressed)
        self.assertEqual(len(decompressed), 3)  # Should return 3-tuple
        
        decomp_mzs, decomp_intensities, decomp_charges = decompressed
        self.assertEqual(len(decomp_mzs), len(mzs))
        self.assertEqual(len(decomp_intensities), len(intensities))
        self.assertEqual(len(decomp_charges), len(charges))
        
        for orig, decomp in zip(mzs, decomp_mzs):
            self.assertAlmostEqual(orig, decomp, places=PRECISION_TOL)
        
        self.assertEqual(decomp_charges, charges)
    
    def test_different_compression_methods(self):
        """Test that different compression methods all work."""
        for method in ['gzip', 'zlib', 'brotli']:
            try:
                compressed = pt.compress_spectra(self.sample_spectra, compression=method)
                self.assertIsInstance(compressed, str)
                
                decompressed = pt.decompress_spectra(compressed)
                mzs, intensities = decompressed
                
                self.assertEqual(len(mzs), len(self.sample_mzs))
                for orig, decomp in zip(self.sample_mzs, mzs):
                    self.assertAlmostEqual(orig, decomp, places=PRECISION_TOL)
            except ImportError:
                if method == 'brotli':
                    self.skipTest(f"{method} not available")


class TestPeptideCompression(unittest.TestCase):
    
    def setUp(self):
        self.sample_peptides = [
            "PEPTIDE", "PEPTIDE/2", "PEPTIDE/3", "PEPTIDE/2",
            "REPFYD", "REPFYD/3", "MGLSDGEWQQVLNVWGK"
        ]
    
    def test_basic_compression_round_trip(self):
        compressed = pt.compress_peptides(self.sample_peptides)
        decompressed = pt.decompress_peptides(compressed)
        
        # Sort both lists since order might change
        self.assertEqual(sorted(self.sample_peptides), sorted(decompressed))
    
    def test_gzip_compression(self):
        compressed = pt.compress_peptides(self.sample_peptides, compression_method='gzip')
        self.assertTrue(compressed.startswith('SG') or compressed.startswith('UG'))
        
        decompressed = pt.decompress_peptides(compressed)
        self.assertEqual(sorted(self.sample_peptides), sorted(decompressed))
    
    def test_zlib_compression(self):
        compressed = pt.compress_peptides(self.sample_peptides, compression_method='zlib')
        self.assertTrue(compressed.startswith('SZ') or compressed.startswith('UZ'))
        
        decompressed = pt.decompress_peptides(compressed)
        self.assertEqual(sorted(self.sample_peptides), sorted(decompressed))
    
    def test_brotli_compression(self):
        try:
            compressed = pt.compress_peptides(self.sample_peptides, compression_method='brotli')
            self.assertTrue(compressed.startswith('SB') or compressed.startswith('UB'))
            
            decompressed = pt.decompress_peptides(compressed)
            self.assertEqual(sorted(self.sample_peptides), sorted(decompressed))
        except ImportError:
            self.skipTest("Brotli not available")
    
    def test_url_safe_encoding(self):
        compressed = pt.compress_peptides(self.sample_peptides, url_safe=True)
        self.assertTrue(compressed.startswith('U'))
        
        decompressed = pt.decompress_peptides(compressed)
        self.assertEqual(sorted(self.sample_peptides), sorted(decompressed))
    
    def test_standard_encoding(self):
        compressed = pt.compress_peptides(self.sample_peptides, url_safe=False)
        self.assertTrue(compressed.startswith('S'))
        
        decompressed = pt.decompress_peptides(compressed)
        self.assertEqual(sorted(self.sample_peptides), sorted(decompressed))
    
    def test_precision_parameter(self):
        compressed_p3 = pt.compress_peptides(self.sample_peptides, precision=3)
        compressed_p1 = pt.compress_peptides(self.sample_peptides, precision=1)
        
        decompressed_p3 = pt.decompress_peptides(compressed_p3)
        decompressed_p1 = pt.decompress_peptides(compressed_p1)
        
        self.assertEqual(sorted(self.sample_peptides), sorted(decompressed_p3))
        self.assertEqual(sorted(self.sample_peptides), sorted(decompressed_p1))
    
    def test_empty_peptide_list(self):
        empty_peptides: list[str] = []
        compressed = pt.compress_peptides(empty_peptides)
        decompressed = pt.decompress_peptides(compressed)
        
        self.assertEqual(decompressed, [])
    
    def test_single_peptide(self):
        single_peptide = ["PEPTIDE/2"]
        compressed = pt.compress_peptides(single_peptide)
        decompressed = pt.decompress_peptides(compressed)
        
        self.assertEqual(decompressed, single_peptide)
    
    def test_duplicate_peptides_handled(self):
        # Test that duplicate peptides are properly counted and reconstructed
        peptides_with_dups = ["PEPTIDE/2", "PEPTIDE/2", "PEPTIDE/2", "REPFYD"]
        compressed = pt.compress_peptides(peptides_with_dups)
        decompressed = pt.decompress_peptides(compressed)
        
        self.assertEqual(sorted(peptides_with_dups), sorted(decompressed))
        self.assertEqual(decompressed.count("PEPTIDE/2"), 3)
        self.assertEqual(decompressed.count("REPFYD"), 1)
    


if __name__ == '__main__':
    unittest.main()