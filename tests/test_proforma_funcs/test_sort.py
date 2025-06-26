# filepath: /workspaces/peptacular/tests/test_proforma_funcs/test_sort.py
import unittest

import peptacular as pt
from peptacular.proforma_dataclasses import Interval, Mod


class TestSort(unittest.TestCase):
    
    def test_basic_sort(self):
        annotation = pt.parse("PEPTIDE")
        sorted_annotation = annotation.sort()
        self.assertEqual(sorted_annotation.serialize(), "DEEIPPT")
        self.assertEqual(annotation.serialize(), "PEPTIDE")  # Original unchanged

    def test_basic_sort_inplace(self):
        # Test sorting in place
        annotation = pt.parse("PEPTIDE")
        sorted_annotation = annotation.sort(inplace=True)
        self.assertEqual(sorted_annotation.serialize(), "DEEIPPT")
        self.assertEqual(annotation.serialize(), "DEEIPPT")  # Original changed too

    def test_sort_single_amino_acid(self):
        annotation = pt.ProFormaAnnotation(sequence="A")
        sorted_annotation = annotation.sort()
        self.assertEqual(sorted_annotation.serialize(), "A")

    def test_sort_empty_sequence(self):
        annotation = pt.ProFormaAnnotation(sequence="")
        sorted_annotation = annotation.sort()
        self.assertEqual(sorted_annotation.serialize(), "")

    def test_sort_already_sorted(self):
        annotation = pt.parse("DEEIPPT")
        sorted_annotation = annotation.sort()
        self.assertEqual(sorted_annotation.serialize(), "DEEIPPT")

    def test_sort_reverse_order(self):
        annotation = pt.parse("PEPTIDE")
        sorted_annotation = annotation.sort(reverse=True)
        self.assertEqual(sorted_annotation.serialize(), "TPPIEED")

    """
    TESTS FOR: sorting with custom key function
    """

    def test_sort_with_reverse_alphabetical_key(self):
        annotation = pt.parse("ABCD")
        # Sort in reverse alphabetical order
        sorted_annotation = annotation.sort(key=lambda x: (-ord(x),))
        self.assertEqual(sorted_annotation.serialize(), "DCBA")

    def test_sort_with_custom_ordering_key(self):
        # Custom order: vowels first, then consonants
        def vowel_first_key(aa):
            vowels = "AEIOU"
            if aa in vowels:
                return (0, aa)  # Vowels get priority 0
            else:
                return (1, aa)  # Consonants get priority 1
        
        annotation = pt.parse("PEPTIDE")
        sorted_annotation = annotation.sort(key=vowel_first_key)
        # Should have E, E, I first (vowels), then D, P, P, T (consonants)
        self.assertEqual(sorted_annotation.serialize(), "EEIDPPT")

    """
    TESTS FOR: sorting with internal modifications
    """
    def test_sort_with_internal_mods(self):
        annotation = pt.parse("PE[Phospho]PTI[Methyl]DE")
        sorted_annotation = annotation.sort()
        # D, E, E, P, P, T, I -> positions change accordingly
        self.assertEqual(sorted_annotation.serialize(), "DE[Phospho]EI[Methyl]PPT")

    def test_sort_with_multiple_internal_mods_same_position(self):
        annotation = pt.parse("PE[Phospho][Methyl]PTIDE")
        sorted_annotation = annotation.sort()
        self.assertEqual(sorted_annotation.serialize(), "DE[Phospho][Methyl]EIPPT")

    def test_sort_with_multiple_internal_mods_different_positions(self):
        annotation = pt.parse("P[Mod1]E[Mod2]PTIDE")
        sorted_annotation = annotation.sort()
        self.assertEqual(sorted_annotation.serialize(), "DE[Mod2]EIP[Mod1]PT")

    """
    TESTS FOR: sorting with terminal modifications
    """
    def test_sort_with_nterm_mods(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE")
        sorted_annotation = annotation.sort()
        serialized = sorted_annotation.serialize()
        self.assertIn("[Acetyl]-", serialized)
        self.assertEqual(sorted_annotation.sequence, "DEEIPPT")

    def test_sort_with_cterm_mods(self):
        annotation = pt.parse("PEPTIDE-[Amide]")
        sorted_annotation = annotation.sort()
        serialized = sorted_annotation.serialize()
        self.assertIn("-[Amide]", serialized)
        self.assertEqual(sorted_annotation.sequence, "DEEIPPT")

    def test_sort_with_both_term_mods(self):
        annotation = pt.parse("[Acetyl]-PEPTIDE-[Amide]")
        sorted_annotation = annotation.sort()
        serialized = sorted_annotation.serialize()
        self.assertIn("[Acetyl]-", serialized)
        self.assertIn("-[Amide]", serialized)
        self.assertEqual(sorted_annotation.sequence, "DEEIPPT")

    """
    TESTS FOR: sorting with other modification types
    """
    def test_sort_with_labile_mods(self):
        annotation = pt.parse("{Glycan}PEPTIDE")
        sorted_annotation = annotation.sort()
        serialized = sorted_annotation.serialize()
        self.assertIn("{Glycan}", serialized)
        self.assertEqual(sorted_annotation.sequence, "DEEIPPT")

    def test_sort_with_static_mods(self):
        annotation = pt.parse("<57@C>PEPTIDE")
        sorted_annotation = annotation.sort()
        serialized = sorted_annotation.serialize()
        self.assertIn("<57@C>", serialized)
        self.assertEqual(sorted_annotation.sequence, "DEEIPPT")

    def test_sort_with_charge(self):
        annotation = pt.parse("PEPTIDE/2")
        sorted_annotation = annotation.sort()
        serialized = sorted_annotation.serialize()
        self.assertIn("/2", serialized)
        self.assertEqual(sorted_annotation.sequence, "DEEIPPT")

    def test_sort_with_unknown_mod(self):
        annotation = pt.parse("[Unknown]?PEPTIDE")
        sorted_annotation = annotation.sort()
        serialized = sorted_annotation.serialize()
        self.assertIn("[Unknown]?", serialized)
        self.assertEqual(sorted_annotation.sequence, "DEEIPPT")

    """
    TESTS FOR: edge cases and validation
    """
    def test_sort_preserves_amino_acid_composition(self):
        annotation = pt.parse("PEPTIDE")
        sorted_annotation = annotation.sort()
        # Check that same amino acids are present
        original_sorted = sorted(annotation.sequence)
        sorted_sorted = sorted(sorted_annotation.sequence)
        self.assertEqual(original_sorted, sorted_sorted)

    def test_sort_with_duplicate_amino_acids(self):
        annotation = pt.parse("AAABBBCCC")
        sorted_annotation = annotation.sort()
        self.assertEqual(sorted_annotation.sequence, "AAABBBCCC")

    def test_sort_with_key_none_explicit(self):
        annotation = pt.parse("PEPTIDE")
        sorted_annotation = annotation.sort(key=None)
        self.assertEqual(sorted_annotation.sequence, "DEEIPPT")


if __name__ == '__main__':
    unittest.main()