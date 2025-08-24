# filepath: /workspaces/peptacular/tests/test_proforma_funcs/test_reverse.py
import unittest

import peptacular as pt


class TestReverse(unittest.TestCase):
    
    def test_basic_reverse(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.sequence, "EDITPEP")
        self.assertEqual(annotation.sequence, "PEPTIDE")  # Original unchanged

    def test_basic_reverse_inplace(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        reversed_annotation = annotation.reverse(inplace=True)
        self.assertEqual(reversed_annotation.sequence, "EDITPEP")
        self.assertEqual(annotation.sequence, "EDITPEP")  # Original changed too

    def test_reverse_single_amino_acid(self):
        annotation = pt.ProFormaAnnotation(sequence="A")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.sequence, "A")

    def test_reverse_empty_sequence(self):
        annotation = pt.ProFormaAnnotation(sequence="")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.sequence, "")

    """
    TESTS FOR: reversing with internal modifications
    """
    def test_reverse_with_internal_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "EDI[Methyl]TPE[Phospho]P")

    def test_reverse_with_multiple_internal_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho][Methyl]PTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "EDITPE[Phospho][Methyl]P")

    """
    TESTS FOR: reversing with terminal modifications
    """

    def test_reverse_with_both_term_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "[Acetyl]-EDITPEP-[Amide]")

    def test_reverse_with_swap_terms_true(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        reversed_annotation = annotation.reverse(swap_terms=True)
        self.assertEqual(reversed_annotation.serialize(), "[Amide]-EDITPEP-[Acetyl]")

    """
    TESTS FOR: reversing with other modification types
    """
    def test_reverse_with_labile_mods(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "{Glycan}EDITPEP")

    def test_reverse_with_static_mods(self):
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "<57@C>EDITPEP")

    def test_reverse_with_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "EDITPEP/2")

    def test_reverse_with_unknown_mod(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "[Unknown]?EDITPEP")


if __name__ == '__main__':
    unittest.main()