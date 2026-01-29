import unittest

import peptacular as pt


class TestJoin(unittest.TestCase):
    def test_join_basic_sequence(self):
        """Test joining basic split sequence back together"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.sequence, "PEPTIDE")
        self.assertEqual(joined.serialize(), "PEPTIDE")

    def test_join_empty_list_raises_error(self):
        """Test that joining empty list raises ValueError"""
        with self.assertRaises(ValueError):
            pt.ProFormaAnnotation.join([])

    def test_join_single_annotation(self):
        """Test joining single annotation returns copy"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        joined = pt.ProFormaAnnotation.join([annotation])

        self.assertEqual(joined.serialize(), "PEPTIDE")
        # Should be a copy, not the same object
        self.assertIsNot(joined, annotation)

    def test_join_preserves_terminal_modifications(self):
        """Test joining preserves N-term from first and C-term from last"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "[Acetyl]-PEPTIDE-[Amide]")

    def test_join_preserves_internal_modifications(self):
        """Test joining preserves internal modifications at correct positions"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "PE[Phospho]PTI[Methyl]DE")

    def test_join_preserves_labile_modifications(self):
        """Test joining preserves labile modifications"""
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "{Glycan}PEPTIDE")

    def test_join_preserves_static_modifications(self):
        """Test joining preserves static modifications"""
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "<57@C>PEPTIDE")

    def test_join_preserves_isotope_modifications(self):
        """Test joining preserves isotope modifications"""
        annotation = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "<15N>PEPTIDE")

    def test_join_preserves_charge_state(self):
        """Test joining preserves charge state"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "PEPTIDE/2")

    def test_join_preserves_unknown_modifications(self):
        """Test joining preserves unknown modifications"""
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "[Unknown]?PEPTIDE")

    def test_join_with_intervals(self):
        """Test joining preserves intervals at correct positions"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "P(EP)[Phospho]TIDE")

    def test_join_with_ambiguous_intervals(self):
        """Test joining preserves ambiguous intervals"""
        annotation = pt.ProFormaAnnotation.parse("P(?EPTI)[Phospho]DE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "P(?EPTI)[Phospho]DE")

    def test_join_with_multiple_modifications(self):
        """Test joining complex annotation with multiple modification types"""
        annotation = pt.ProFormaAnnotation.parse("<15N>{Glycan}[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "<15N>{Glycan}[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")

    def test_join_with_stacked_modifications(self):
        """Test joining preserves stacked modifications on same residue"""
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho][Methyl]PTIDE")
        split_annotations = annotation.split()
        joined = pt.ProFormaAnnotation.join(split_annotations)

        self.assertEqual(joined.serialize(), "PE[Phospho][Methyl]PTIDE")

    def test_join_partial_split(self):
        """Test joining partial list of split annotations"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        split_annotations = annotation.split()

        # Join only first 3 fragments
        partial_joined = pt.ProFormaAnnotation.join(split_annotations[:3])
        self.assertEqual(partial_joined.sequence, "PEP")

    def test_join_preserves_modifications_in_partial_join(self):
        """Test joining partial annotations preserves appropriate modifications"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        split_annotations = annotation.split()

        # Join first 3 fragments - should keep N-term but lose C-term
        partial_joined = pt.ProFormaAnnotation.join(split_annotations[:3])
        self.assertEqual(partial_joined.serialize(), "[Acetyl]-PE[Phospho]P")

        # Join last 3 fragments - should keep C-term but lose N-term
        partial_joined = pt.ProFormaAnnotation.join(split_annotations[-3:])
        self.assertEqual(partial_joined.serialize(), "IDE-[Amide]")

    def test_join_split_roundtrip_invariant(self):
        """Test that join(split(annotation)) == annotation for various cases"""
        test_cases = [
            "PEPTIDE",
            "[Acetyl]-PEPTIDE",
            "PEPTIDE-[Amide]",
            "[Acetyl]-PEPTIDE-[Amide]",
            "PE[Phospho]PTIDE",
            "{Glycan}PEPTIDE",
            "<15N>PEPTIDE",
            "PEPTIDE/2",
            "P(EP)[Phospho]TIDE",
            "P(?EPTI)[Phospho]DE",
            "<15N>{Glycan}[Acetyl]-PE[Phospho]PTIDE-[Amide]/2",
            "PE[Phospho][Methyl]PTIDE",
        ]

        for test_case in test_cases:
            with self.subTest(test_case=test_case):
                annotation = pt.ProFormaAnnotation.parse(test_case)
                split_annotations = annotation.split()
                joined = pt.ProFormaAnnotation.join(split_annotations)
                self.assertEqual(joined.serialize(), test_case)

    def test_join_different_annotations(self):
        """Test joining annotations that weren't created from split"""
        ann1 = pt.ProFormaAnnotation.parse("[Acetyl]-PEP")
        ann2 = pt.ProFormaAnnotation.parse("TI[Phospho]")
        ann3 = pt.ProFormaAnnotation.parse("DE-[Amide]")

        joined = pt.ProFormaAnnotation.join([ann1, ann2, ann3])
        self.assertEqual(joined.sequence, "PEPTIDE")
        self.assertEqual(joined.serialize(), "[Acetyl]-PEPTI[Phospho]DE-[Amide]")


if __name__ == "__main__":
    unittest.main()
