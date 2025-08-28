import unittest
import peptacular as pt


class TestSplit(unittest.TestCase):

    def test_split_basic_sequence(self):
        """Test splitting a basic sequence without modifications"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)
        expected_sequences = ["P", "E", "P", "T", "I", "D", "E"]
        for i, split_annot in enumerate(split_annotations):
            self.assertEqual(split_annot.sequence, expected_sequences[i])

    def test_split_single_amino_acid(self):
        """Test splitting a single amino acid sequence"""
        annotation = pt.ProFormaAnnotation.parse("A")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 1)
        self.assertEqual(split_annotations[0].sequence, "A")

    def test_split_empty_sequence(self):
        """Test splitting an empty sequence"""
        annotation = pt.ProFormaAnnotation.parse("")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 0)

    def test_split_preserves_terminal_modifications(self):
        """Test that terminal modifications are preserved in appropriate split fragments"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # First fragment should have N-terminal modification
        self.assertEqual(split_annotations[0].serialize(), "[Acetyl]-P")

        # Last fragment should have C-terminal modification
        self.assertEqual(split_annotations[6].serialize(), "E-[Amide]")

        # Middle fragments should have no terminal modifications
        for i in range(1, 6):
            self.assertFalse(split_annotations[i].has_nterm_mods)
            self.assertFalse(split_annotations[i].has_cterm_mods)

    def test_split_preserves_internal_modifications(self):
        """Test that internal modifications are preserved in correct split fragments"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # Second fragment (E) should have the phospho modification
        self.assertEqual(split_annotations[1].serialize(), "E[Phospho]")

        # Fifth fragment (I) should have the methyl modification
        self.assertEqual(split_annotations[4].serialize(), "I[Methyl]")

        # Other fragments should have no internal modifications
        for i in [0, 2, 3, 5, 6]:
            self.assertFalse(split_annotations[i].has_internal_mods)

    def test_split_preserves_labile_modifications(self):
        """Test that labile modifications are preserved in all split fragments"""
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # All fragments should have the labile modification
        for split_annot in split_annotations:
            self.assertTrue(split_annot.has_labile_mods)
            self.assertEqual(split_annot.labile_mods[0], "Glycan")

    def test_split_preserves_static_modifications(self):
        """Test that static modifications are preserved in all split fragments"""
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # All fragments should have the static modification
        for split_annot in split_annotations:
            self.assertTrue(split_annot.has_static_mods)

    def test_split_preserves_isotope_modifications(self):
        """Test that isotope modifications are preserved in all split fragments"""
        annotation = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # All fragments should have the isotope modification
        for split_annot in split_annotations:
            self.assertTrue(split_annot.has_isotope_mods)

    def test_split_preserves_charge_state(self):
        """Test that charge state is preserved in all split fragments"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # All fragments should have the same charge
        for split_annot in split_annotations:
            self.assertEqual(split_annot.charge, 2)

    def test_split_preserves_unknown_modifications(self):
        """Test that unknown modifications are preserved in all split fragments"""
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # All fragments should have the unknown modification
        for split_annot in split_annotations:
            self.assertTrue(split_annot.has_unknown_mods)

    def test_split_with_multiple_modifications(self):
        """Test splitting with multiple types of modifications"""
        annotation = pt.ProFormaAnnotation.parse(
            "{Glycan}<15N>[Acetyl]-PE[Phospho]PTIDE-[Amide]/2"
        )
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # Check first fragment has all global mods plus N-term
        first = split_annotations[0]
        self.assertTrue(first.has_labile_mods)
        self.assertTrue(first.has_isotope_mods)
        self.assertTrue(first.has_nterm_mods)
        self.assertEqual(first.charge, 2)
        self.assertEqual(first.serialize(), "{Glycan}<15N>[Acetyl]-P/2")

        # Check fragment with internal modification
        second = split_annotations[1]
        self.assertTrue(second.has_internal_mods)
        self.assertEqual(second.serialize(), "{Glycan}<15N>E[Phospho]/2")

        # Check last fragment has C-term
        last = split_annotations[6]
        self.assertTrue(last.has_cterm_mods)
        self.assertEqual(last.serialize(), "{Glycan}<15N>E-[Amide]/2")

    def test_split_with_intervals(self):
        """Test splitting annotation with intervals"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        split_annotations = annotation.split()

        # Should have 5 fragments: P, (EP)[Phospho], T, I, DE
        self.assertEqual(len(split_annotations), 6)

        expected_fragments = ["P", "(EP)[Phospho]", "T", "I", "D", "E"]
        for i, expected in enumerate(expected_fragments):
            self.assertEqual(split_annotations[i].serialize(), expected)

    def test_split_with_ambiguous_intervals(self):
        """Test splitting annotation with ambiguous intervals"""
        annotation = pt.ProFormaAnnotation.parse("P(?EPTI)[Phospho]DE")
        split_annotations = annotation.split()

        # Should have 3 fragments: P, (?EPTI)[Phospho], DE
        self.assertEqual(len(split_annotations), 4)

        expected_fragments = ["P", "(?EPTI)[Phospho]", "D", "E"]
        for i, expected in enumerate(expected_fragments):
            self.assertEqual(split_annotations[i].serialize(), expected)

    def test_split_original_annotation_unchanged(self):
        """Test that splitting doesn't modify the original annotation"""
        original_sequence = "PE[Phospho]PTIDE"
        annotation = pt.ProFormaAnnotation.parse(original_sequence)
        original_serialized = annotation.serialize()

        split_annotations = annotation.split()

        # Original should be unchanged
        self.assertEqual(annotation.serialize(), original_serialized)
        self.assertEqual(len(split_annotations), 7)

    def test_split_with_stacked_modifications(self):
        """Test splitting with stacked modifications on same residue"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho][Methyl]PTIDE")
        split_annotations = annotation.split()

        self.assertEqual(len(split_annotations), 7)

        # Second fragment should have both modifications
        second = split_annotations[1]
        self.assertTrue(second.has_internal_mods)
        # Check that both modifications are present
        internal_mods = second.internal_mods[0]  # Position 0 in the split fragment
        self.assertEqual(len(internal_mods), 2)


if __name__ == "__main__":
    unittest.main()
