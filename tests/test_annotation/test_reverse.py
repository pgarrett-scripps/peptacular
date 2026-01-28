import unittest

import peptacular as pt


class TestReverse(unittest.TestCase):
    def test_basic_reverse(self):
        """Test basic sequence reversal"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        reversed_annotation = annotation.reverse()

        self.assertEqual(reversed_annotation.serialize(), "EDITPEP")
        self.assertEqual(annotation.serialize(), "PEPTIDE")  # Original unchanged

    def test_basic_reverse_inplace(self):
        """Test in-place sequence reversal"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        reversed_annotation = annotation.reverse(inplace=True)

        self.assertEqual(reversed_annotation.serialize(), "EDITPEP")
        self.assertEqual(annotation.serialize(), "EDITPEP")  # Original changed too

    def test_reverse_single_amino_acid(self):
        """Test reversing single amino acid"""
        annotation = pt.ProFormaAnnotation.parse("A")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "A")

    def test_reverse_empty_sequence(self):
        """Test reversing empty sequence"""
        annotation = pt.ProFormaAnnotation.parse("")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "")

    def test_reverse_with_internal_mods(self):
        """Test reversing with internal modifications"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "EDI[Methyl]TPE[Phospho]P")

    def test_reverse_with_multiple_internal_mods_same_position(self):
        """Test reversing with multiple modifications on same position"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho][Methyl]PTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "EDITPE[Phospho][Methyl]P")

    def test_reverse_with_intervals(self):
        """Test reversing with intervals"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        reversed_annotation = annotation.reverse()
        # Original interval was positions 1-3, reversed should be positions 3-5 (7-1-3=3, 7-1-1=5)
        self.assertEqual(reversed_annotation.serialize(), "EDIT(PE)[Phospho]P")

    def test_reverse_with_ambiguous_intervals(self):
        """Test reversing with ambiguous intervals"""
        annotation = pt.ProFormaAnnotation.parse("P(?EPTI)[Phospho]DE")
        reversed_annotation = annotation.reverse()
        # Original interval was positions 1-5, reversed should be positions 1-5 (7-1-5=1, 7-1-1=5)
        self.assertEqual(reversed_annotation.serialize(), "ED(?ITPE)[Phospho]P")

    def test_reverse_with_multiple_intervals(self):
        """Test reversing with multiple intervals"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TI(DE)[Methyl]")
        reversed_annotation = annotation.reverse()
        # First interval (1,3) becomes (4,6), second interval (4,6) becomes (1,3)
        self.assertEqual(
            reversed_annotation.serialize(), "(ED)[Methyl]IT(PE)[Phospho]P"
        )

    def test_reverse_with_labile_mods(self):
        """Test reversing with labile modifications"""
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "{Glycan}EDITPEP")

    def test_reverse_with_static_mods(self):
        """Test reversing with static modifications"""
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "<57@C>EDITPEP")

    def test_reverse_with_isotope_mods(self):
        """Test reversing with isotope modifications"""
        annotation = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "<15N>EDITPEP")

    def test_reverse_with_charge(self):
        """Test reversing with charge state"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "EDITPEP/2")

    def test_reverse_with_unknown_mod(self):
        """Test reversing with unknown modifications"""
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "[Unknown]?EDITPEP")

    def test_reverse_with_charge_adducts(self):
        """Test reversing with charge adducts"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/[Na:z+1]")
        reversed_annotation = annotation.reverse()
        self.assertEqual(reversed_annotation.serialize(), "EDITPEP/[Na:z+1]")

    def test_reverse_with_complex_modifications(self):
        """Test reversing with multiple types of modifications"""
        annotation = pt.ProFormaAnnotation.parse(
            "<15N>{Glycan}[Acetyl]-PE[Phospho]PTI[Methyl]DE-[Amide]/2"
        )
        reversed_annotation = annotation.reverse()
        self.assertEqual(
            reversed_annotation.serialize(),
            "<15N>{Glycan}[Acetyl]-EDI[Methyl]TPE[Phospho]P-[Amide]/2",
        )

    def test_reverse_preserves_original_when_not_inplace(self):
        """Test that original annotation is preserved when not in-place"""
        original_annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTIDE")
        original_serialized = original_annotation.serialize()

        reversed_annotation = original_annotation.reverse()

        self.assertEqual(original_annotation.serialize(), original_serialized)
        self.assertEqual(reversed_annotation.serialize(), "EDITPE[Phospho]P")


if __name__ == "__main__":
    unittest.main()
