"""Tests for funcs.py module"""

import peptacular as pt
import unittest

class TestDeltaMassAdjustement(unittest.TestCase):
    """Test ModLabler class"""

    def test_delta_single_value(self):
        """Test delta mass adjustment with single float value"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass1 = annotation.mass()
        mass2 = annotation.mass(delta=15.99)
        self.assertAlmostEqual(mass1 + 15.99, mass2, places=6)

    def test_delta_list_values(self):
        """Test delta mass adjustment with list of float values"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass1 = annotation.mass()
        deltas = [-15.99, 15.99]
        mass2 = annotation.mass(delta=deltas)
        self.assertAlmostEqual(mass1, mass2, places=6)

    def test_delta_dict_values(self):
        """Test delta mass adjustment with dict of float values"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass1 = annotation.mass()
        deltas = {"C": 1} # add a single Carbon
        mass2 = annotation.mass(delta=deltas)
        self.assertAlmostEqual(mass1 + 12.0, mass2, places=6)

    def test_delta_mixed_values(self):
        """Test delta mass adjustment with mixed iterable values"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass1 = annotation.mass()
        deltas = [{"C": 1}, -12] 
        mass2 = annotation.mass(delta=deltas)
        self.assertAlmostEqual(mass1, mass2, places=6)

    def test_str_formula_delta(self):
        """Test delta mass adjustment with string formula"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass1 = annotation.mass()
        mass2 = annotation.mass(delta="H2O")  # adding water
        self.assertAlmostEqual(mass1 + 18.01056, mass2, places=5)

    def test_str_mod_delta(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        mass1 = annotation.mass()
        deltas = ["Oxidation"]
        mass2 = annotation.mass(delta=deltas)
        self.assertAlmostEqual(mass1 + 15.99491, mass2, places=5)