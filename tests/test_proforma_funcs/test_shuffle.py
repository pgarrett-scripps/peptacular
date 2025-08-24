# filepath: /workspaces/peptacular/tests/test_proforma_funcs/test_shuffle.py
import unittest

import peptacular as pt

class TestShuffle(unittest.TestCase):
    
    def test_basic_shuffle(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "ETIPEPD")

    def test_basic_shuffle_inplace(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "ETIPEPD")
        self.assertEqual(annotation.serialize(), "PEPTIDE")  # Original unchanged

    def test_shuffle_single_amino_acid(self):
        annotation = pt.ProFormaAnnotation(sequence="A")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "A")

    def test_shuffle_empty_sequence(self):
        annotation = pt.ProFormaAnnotation(sequence="")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "")

    """
    TESTS FOR: shuffling with internal modifications
    """
    def test_shuffle_with_internal_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "E[Phospho]TI[Methyl]PEPD")

    """
    TESTS FOR: shuffling with terminal modifications
    """

    def test_shuffle_with_both_term_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "[Acetyl]-ETIPEPD-[Amide]")

    """
    TESTS FOR: shuffling with other modification types
    """
    def test_shuffle_with_labile_mods(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "{Glycan}ETIPEPD")

    def test_shuffle_with_static_mods(self):
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "<57@C>ETIPEPD")

    def test_shuffle_with_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "ETIPEPD/2")

    def test_shuffle_with_unknown_mod(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        shuffled_annotation = annotation.shuffle(seed=42)
        self.assertEqual(shuffled_annotation.serialize(), "[Unknown]?ETIPEPD")

    def test_shuffle_preserves_amino_acid_composition(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shuffled_annotation = annotation.shuffle(seed=42)
        # Check that same amino acids are present
        original_sorted = sorted(annotation.sequence)
        shuffled_sorted = sorted(shuffled_annotation.sequence)
        self.assertEqual(original_sorted, shuffled_sorted)


if __name__ == '__main__':
    unittest.main()