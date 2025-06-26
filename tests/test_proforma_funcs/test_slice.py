
import unittest

import peptacular as pt
from peptacular.proforma_dataclasses import Interval, Mod


class TestSlice(unittest.TestCase):
    
    def test_basic_slice(self):
        # Test slicing a ProFormaAnnotation
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        sliced_annotation = annotation.slice(0, 3)
        self.assertEqual(sliced_annotation.sequence, "PEP")
        self.assertEqual(annotation.sequence, "PEPTIDE")

    def test_basic_slice_inplace(self):
        # Test slicing a ProFormaAnnotation
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE")
        sliced_annotation = annotation.slice(0, 3, inplace=True)
        self.assertEqual(sliced_annotation.sequence, "PEP")
        self.assertEqual(annotation.sequence, "PEP")
        self.assertIs(sliced_annotation, annotation)

    """
    TESTS FOR: slicing n and c terms
    """
    def test_slice_nterm(self):
        # Test slicing a ProFormaAnnotation with N-terminal modification
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", nterm_mods=Mod("Acetyl"))
        sliced_annotation = annotation.slice(0, 3)
        self.assertEqual(sliced_annotation.sequence, "PEP")
        self.assertEqual(sliced_annotation.serialize(), "[Acetyl]-PEP")
        self.assertEqual(annotation.sequence, "PEPTIDE")

    def test_slice_cterm(self):
        # Test slicing a ProFormaAnnotation with C-terminal modification
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", cterm_mods=Mod("Amidated"))
        sliced_annotation = annotation.slice(3, 7)
        self.assertEqual(sliced_annotation.sequence, "TIDE")
        self.assertEqual(sliced_annotation.serialize(), "TIDE-[Amidated]")
        self.assertEqual(annotation.sequence, "PEPTIDE")

    def test_terms_middle_slice(self):
        # Test slicing a ProFormaAnnotation with modifications in the middle
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           nterm_mods=Mod("Acetyl"), 
                                           cterm_mods=Mod("Amidated"))
        sliced_annotation = annotation.slice(1, 5)
        self.assertEqual(sliced_annotation.sequence, "EPTI")
        self.assertEqual(sliced_annotation.serialize(), "EPTI")
        self.assertEqual(annotation.sequence, "PEPTIDE")

    def test_terms_middle_slice_with_keep_terms(self):
        # Test slicing a ProFormaAnnotation with modifications
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           nterm_mods=Mod("Acetyl"), 
                                           cterm_mods=Mod("Amidated"))
        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE-[Amidated]")
        sliced_annotation = annotation.slice(1, 5, keep_terms=True)
        self.assertEqual(sliced_annotation.sequence, "EPTI")
        self.assertEqual(sliced_annotation.serialize(), "[Acetyl]-EPTI-[Amidated]")
        self.assertEqual(annotation.sequence, "PEPTIDE")
        
    """
    TESTS FOR: slice_intervals=True (default)
    """
    def test_slice_with_intervals(self):
        # Test slicing a ProFormaAnnotation with intervals
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 3, False, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(EP)[Phosphorylated]TIDE")

        sliced_annotation = annotation.slice(0, 5)
        self.assertEqual(sliced_annotation.sequence, "PEPTI")
        self.assertEqual(sliced_annotation.serialize(), "P(EP)[Phosphorylated]TI")

    def test_slice_start_interval(self):
        # Test slicing a ProFormaAnnotation with a cutoff
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 3, False, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(EP)[Phosphorylated]TIDE")

        sliced_annotation = annotation.slice(2, 5)
        self.assertEqual(sliced_annotation.sequence, "PTI")
        self.assertEqual(sliced_annotation.serialize(), "(P)[Phosphorylated]TI")

    def test_slice_end_interval(self):
        # Test slicing a ProFormaAnnotation with a cutoff
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 3, False, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(EP)[Phosphorylated]TIDE")

        sliced_annotation = annotation.slice(0, 2)
        self.assertEqual(sliced_annotation.sequence, "PE")
        self.assertEqual(sliced_annotation.serialize(), "P(E)[Phosphorylated]")

    def test_slice_start_and_end_interval(self):
        # Test slicing a ProFormaAnnotation with a cutoff
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 5, False, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(EPTI)[Phosphorylated]DE")
        sliced_annotation = annotation.slice(2, 4)
        self.assertEqual(sliced_annotation.sequence, "PT")
        self.assertEqual(sliced_annotation.serialize(), "(PT)[Phosphorylated]")

    def test_slice_start_and_end_interval_ambiguous(self):
        # Test slicing a ProFormaAnnotation with a cutoff
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 5, True, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(?EPTI)[Phosphorylated]DE")
        sliced_annotation = annotation.slice(2, 4)
        self.assertEqual(sliced_annotation.sequence, "PT")
        self.assertEqual(sliced_annotation.serialize(), "(?PT)[Phosphorylated]")

    """
    TESTS FOR: slice_intervals=False
    Will throw error if slice cuts an interval
    """
    def test_slice_with_intervals(self):
        # Test slicing a ProFormaAnnotation with intervals
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 3, False, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(EP)[Phosphorylated]TIDE")

        sliced_annotation = annotation.slice(0, 5, slice_intervals=False)
        self.assertEqual(sliced_annotation.sequence, "PEPTI")
        self.assertEqual(sliced_annotation.serialize(), "P(EP)[Phosphorylated]TI")

    def test_slice_start_interval_with_error(self):
        # Test slicing a ProFormaAnnotation with a cutoff
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 3, False, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(EP)[Phosphorylated]TIDE")

        with self.assertRaises(ValueError):
            annotation.slice(2, 5, slice_intervals=False)

    def test_slice_end_interval_with_error(self):
        # Test slicing a ProFormaAnnotation with a cutoff
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 3, False, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(EP)[Phosphorylated]TIDE")

        with self.assertRaises(ValueError):
            annotation.slice(0, 2, slice_intervals=False)

    def test_slice_start_and_end_interval_with_error(self):
        # Test slicing a ProFormaAnnotation with a cutoff
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           intervals=Interval(1, 5, False, "Phosphorylated"))
        self.assertEqual(annotation.serialize(), "P(EPTI)[Phosphorylated]DE")
        
        with self.assertRaises(ValueError):
            annotation.slice(2, 4, slice_intervals=False)

    """
    Slicing with labile mods
    """
    def test_slice_with_labile_mods(self):
        # Test slicing a ProFormaAnnotation with labile modifications
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           labile_mods=Mod("LabileMod"))
        sliced_annotation = annotation.slice(0, 3)
        self.assertEqual(sliced_annotation.sequence, "PEP")
        self.assertEqual(sliced_annotation.serialize(), "{LabileMod}PEP")


    def test_slice_with_keep_labile_mods_false(self):
        # Test slicing a ProFormaAnnotation with labile modifications
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           labile_mods=Mod("LabileMod"))
        sliced_annotation = annotation.slice(0, 3, keep_labile=False)
        self.assertEqual(sliced_annotation.sequence, "PEP")
        self.assertEqual(sliced_annotation.serialize(), "PEP")

    """
    Slicing with other mods
    """
    def test_slice_with_static_mods(self):
        # Test slicing a ProFormaAnnotation with static modifications
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", 
                                           static_mods=Mod("57@C"))
        sliced_annotation = annotation.slice(0, 3)
        self.assertEqual(sliced_annotation.sequence, "PEP")
        self.assertEqual(sliced_annotation.serialize(), "<57@C>PEP")

    def test_slice_with_charge(self):
        # Test slicing a ProFormaAnnotation with charge
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", charge=2)
        sliced_annotation = annotation.slice(0, 3)
        self.assertEqual(sliced_annotation.sequence, "PEP")
        self.assertEqual(sliced_annotation.serialize(), "PEP/2")

    