"""Tests for SequenceElement and SequenceRegion"""

import unittest
from peptacular.components import (
    SequenceElement,
    SequenceRegion,
)
from peptacular.constants import AminoAcid

import peptacular as pt


class TestSequenceElement(unittest.TestCase):
    """Test cases for SequenceElement"""

    def test_simple_amino_acid(self):
        """Test sequence element with just an amino acid"""
        se = SequenceElement(amino_acid=AminoAcid.M)
        self.assertEqual(se.amino_acid, AminoAcid.M)
        self.assertEqual(len(se.modifications), 0)
        self.assertEqual(str(se), "M")

    def test_amino_acid_with_single_modification(self):
        """Test sequence element with one modification"""
        mod = pt.ModificationTags.from_string("Oxidation")
        se = SequenceElement(amino_acid=AminoAcid.M, modifications=(mod,))
        self.assertEqual(se.amino_acid, AminoAcid.M)
        self.assertEqual(len(se.modifications), 1)
        self.assertEqual(str(se), "M[Oxidation]")


class TestSequenceRegion(unittest.TestCase):
    """Test cases for SequenceRegion"""

    def test_simple_sequence_region(self):
        """Test sequence region with no modifications"""
        se1 = SequenceElement(amino_acid=AminoAcid.P)
        se2 = SequenceElement(amino_acid=AminoAcid.E)
        se3 = SequenceElement(amino_acid=AminoAcid.P)

        sr = SequenceRegion(sequence=(se1, se2, se3), modifications=())

        self.assertEqual(len(sr.sequence), 3)
        self.assertEqual(len(sr.modifications), 0)
        self.assertEqual(str(sr), "(PEP)")

    def test_sequence_region_with_modification(self):
        """Test sequence region with one modification"""
        se1 = SequenceElement(amino_acid=AminoAcid.P)
        se2 = SequenceElement(amino_acid=AminoAcid.E)
        se3 = SequenceElement(amino_acid=AminoAcid.P)

        mod = pt.ModificationTags.from_string("Phospho")

        sr = SequenceRegion(sequence=(se1, se2, se3), modifications=(mod,))

        self.assertEqual(len(sr.sequence), 3)
        self.assertEqual(len(sr.modifications), 1)
        self.assertEqual(str(sr), "(PEP)[Phospho]")
