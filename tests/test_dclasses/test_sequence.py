"""Tests for pt.SequenceElement and pt.SequenceRegion"""

import unittest
import peptacular as pt


class TestSequenceElement(unittest.TestCase):
    """Test cases for pt.SequenceElement"""

    def test_simple_amino_acid(self):
        """Test sequence element with just an amino acid"""
        se = pt.SequenceElement(amino_acid=pt.AminoAcid.M)
        self.assertEqual(se.amino_acid, pt.AminoAcid.M)
        self.assertEqual(len(se.modifications), 0)
        self.assertEqual(str(se), "M")

    def test_amino_acid_with_single_modification(self):
        """Test sequence element with one modification"""
        mod = pt.ModificationTags.from_string("Oxidation")
        se = pt.SequenceElement(amino_acid=pt.AminoAcid.M, modifications=(mod,))
        self.assertEqual(se.amino_acid, pt.AminoAcid.M)
        self.assertEqual(len(se.modifications), 1)
        self.assertEqual(str(se), "M[Oxidation]")


class TestSequenceRegion(unittest.TestCase):
    """Test cases for pt.SequenceRegion"""

    def test_simple_sequence_region(self):
        """Test sequence region with no modifications"""
        se1 = pt.SequenceElement(amino_acid=pt.AminoAcid.P)
        se2 = pt.SequenceElement(amino_acid=pt.AminoAcid.E)
        se3 = pt.SequenceElement(amino_acid=pt.AminoAcid.P)

        sr = pt.SequenceRegion(sequence=(se1, se2, se3), modifications=())

        self.assertEqual(len(sr.sequence), 3)
        self.assertEqual(len(sr.modifications), 0)
        self.assertEqual(str(sr), "(PEP)")

    def test_sequence_region_with_modification(self):
        """Test sequence region with one modification"""
        se1 = pt.SequenceElement(amino_acid=pt.AminoAcid.P)
        se2 = pt.SequenceElement(amino_acid=pt.AminoAcid.E)
        se3 = pt.SequenceElement(amino_acid=pt.AminoAcid.P)

        mod = pt.ModificationTags.from_string("Phospho")

        sr = pt.SequenceRegion(sequence=(se1, se2, se3), modifications=(mod,))

        self.assertEqual(len(sr.sequence), 3)
        self.assertEqual(len(sr.modifications), 1)
        self.assertEqual(str(sr), "(PEP)[Phospho]")
