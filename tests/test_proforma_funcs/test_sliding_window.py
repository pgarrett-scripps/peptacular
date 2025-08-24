
from types import GeneratorType
import unittest

import peptacular as pt


class TestSlidingWindow(unittest.TestCase):
    
    def test_basic_sliding_window(self):
        # Test slicing a ProFormaAnnotation
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        
        sliding_windows = annotation.sliding_windows(5)

        #ensure generator
        self.assertIsInstance(sliding_windows, GeneratorType)

        for i, window in enumerate(sliding_windows):
            if i == 0:
                self.assertEqual(window.sequence, "PEPTI")
            elif i == 1:
                self.assertEqual(window.sequence, "EPTID")
            elif i == 2:
                self.assertEqual(window.sequence, "PTIDE")
            else:
                break

    def test_reverse_sliding_window(self):
        # Test reverse sliding windows
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        
        sliding_windows = annotation.sliding_windows(5, reverse=True)

        #ensure generator
        self.assertIsInstance(sliding_windows, GeneratorType)

        for i, window in enumerate(sliding_windows):
            if i == 0:
                self.assertEqual(window.sequence, "PTIDE")
            elif i == 1:
                self.assertEqual(window.sequence, "EPTID")
            elif i == 2:
                self.assertEqual(window.sequence, "PEPTI")
            else:
                break