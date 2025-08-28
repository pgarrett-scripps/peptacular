import unittest

import peptacular as pt


class TestCoverage(unittest.TestCase):

    def test_coverage_basic_exact_match(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [pt.ProFormaAnnotation.parse("PEPTIDE")]
        result = target.coverage(annotations)

        self.assertEqual(result, [1, 1, 1, 1, 1, 1, 1])

    def test_coverage_basic_no_match(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [pt.ProFormaAnnotation.parse("DIFFERENT")]
        result = target.coverage(annotations)

        self.assertEqual(result, [0, 0, 0, 0, 0, 0, 0])

    def test_coverage_subsequence_match(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [pt.ProFormaAnnotation.parse("PEP")]
        result = target.coverage(annotations)

        self.assertEqual(result, [1, 1, 1, 0, 0, 0, 0])

    def test_coverage_multiple_subsequences(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [
            pt.ProFormaAnnotation.parse("PEP"),
            pt.ProFormaAnnotation.parse("TIDE"),
        ]
        result = target.coverage(annotations)

        self.assertEqual(result, [1, 1, 1, 1, 1, 1, 1])

    def test_coverage_overlapping_subsequences(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [
            pt.ProFormaAnnotation.parse("PEPT"),
            pt.ProFormaAnnotation.parse("TIDE"),
        ]
        result = target.coverage(annotations)

        self.assertEqual(result, [1, 1, 1, 1, 1, 1, 1])

    def test_coverage_accumulate_true(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [
            pt.ProFormaAnnotation.parse("PEP"),
            pt.ProFormaAnnotation.parse("PEP"),
        ]  # Same subsequence twice
        result = target.coverage(annotations, accumulate=True)

        self.assertEqual(result, [2, 2, 2, 0, 0, 0, 0])

    def test_coverage_accumulate_false(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [
            pt.ProFormaAnnotation.parse("PEP"),
            pt.ProFormaAnnotation.parse("PEP"),
        ]  # Same subsequence twice
        result = target.coverage(annotations, accumulate=False)

        self.assertEqual(result, [1, 1, 1, 0, 0, 0, 0])

    def test_coverage_ignore_mods_true(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [pt.ProFormaAnnotation.parse("PE[Phospho]PTIDE")]
        result = target.coverage(annotations, ignore_mods=True)

        self.assertEqual(result, [1, 1, 1, 1, 1, 1, 1])

    def test_coverage_ignore_mods_false(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [pt.ProFormaAnnotation.parse("PE[Phospho]PTIDE")]
        result = target.coverage(annotations, ignore_mods=False)

        self.assertEqual(result, [0, 0, 0, 0, 0, 0, 0])

    def test_coverage_with_ambiguous_intervals(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        # Create annotation with ambiguous interval
        annotation = pt.ProFormaAnnotation.parse("P(?EP)TIDE")

        result = target.coverage([annotation], ignore_ambiguity=False)

        # Positions 1 and 2 should have 0 coverage due to ambiguity
        self.assertEqual(result, [1, 0, 0, 1, 1, 1, 1])

    def test_coverage_ignore_ambiguity_true(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        # Create annotation with ambiguous interval
        annotation = pt.ProFormaAnnotation.parse("P(?EP)TIDE")

        result = target.coverage([annotation], ignore_ambiguity=True)

        # All positions should have coverage when ignoring ambiguity
        self.assertEqual(result, [1, 1, 1, 1, 1, 1, 1])

    def test_coverage_empty_annotations_list(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        result = target.coverage([])

        self.assertEqual(result, [0, 0, 0, 0, 0, 0, 0])

    def test_coverage_partial_overlaps(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [
            pt.ProFormaAnnotation.parse("EPT"),
            pt.ProFormaAnnotation.parse("TID"),
        ]
        result = target.coverage(annotations)

        self.assertEqual(result, [0, 1, 1, 1, 1, 1, 0])

    def test_coverage_multiple_overlapping_with_accumulate(self):
        target = pt.ProFormaAnnotation.parse("PEPTIDE")
        annotations = [
            pt.ProFormaAnnotation.parse("PEPT"),
            pt.ProFormaAnnotation.parse("EPT"),
            pt.ProFormaAnnotation.parse("PT"),
        ]
        result = target.coverage(annotations, accumulate=True)

        # P=1, E=2, P=3, T=2, I=0, D=0, E=0
        self.assertEqual(result, [1, 2, 3, 3, 0, 0, 0])

    def test_coverage_with_modifications_in_target(self):
        target = pt.ProFormaAnnotation.parse("PE[Phospho]PTIDE")
        annotations = [pt.ProFormaAnnotation.parse("PEPTIDE")]
        result = target.coverage(annotations, ignore_mods=True)

        self.assertEqual(result, [1, 1, 1, 1, 1, 1, 1])


if __name__ == "__main__":
    unittest.main()
