import unittest
import peptacular as pt


class TestSlice(unittest.TestCase):
    def test_basic_slice(self):
        """Test basic slicing functionality"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")

        # Test normal slicing
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.sequence, "PEP")
        self.assertEqual(annotation.sequence, "PEPTIDE")  # Original unchanged

        # Test in-place slicing
        annotation.slice(0, 3, inplace=True)
        self.assertEqual(annotation.sequence, "PEP")

    def test_slice_preserves_nterm_modification(self):
        """Test slicing from start preserves N-terminal modification"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.serialize(), "[Acetyl]-PEP")

    def test_slice_preserves_cterm_modification(self):
        """Test slicing to end preserves C-terminal modification"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amidated]")
        sliced = annotation.slice(3, 7)
        self.assertEqual(sliced.serialize(), "TIDE-[Amidated]")

    def test_slice_removes_both_terminal_modifications_when_middle_slice(self):
        """Test middle slice removes both terminal modifications"""
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amidated]")
        sliced = annotation.slice(1, 5)
        self.assertEqual(sliced.serialize(), "EPTI")

    def test_slice_preserves_complete_interval(self):
        """Test slicing preserves interval when completely within slice"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")
        self.assertEqual(annotation.serialize(), "P(EP)[Phosphorylated]TIDE")

        sliced = annotation.slice(0, 5)
        self.assertEqual(sliced.serialize(), "P(EP)[Phosphorylated]TI")

    def test_slice_preserves_complete_interval_full_sequence(self):
        """Test slicing preserves interval when taking full sequence"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")

        sliced = annotation.slice(0, 7)
        self.assertEqual(sliced.serialize(), "P(EP)[Phosphorylated]TIDE")

    def test_slice_raises_error_when_cutting_interval_start(self):
        """Test slicing raises ValueError when interval start gets cut"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")

        with self.assertRaises(ValueError):
            annotation.slice(2, 5)

    def test_slice_raises_error_when_cutting_interval_end(self):
        """Test slicing raises ValueError when interval end gets cut"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")

        with self.assertRaises(ValueError):
            annotation.slice(0, 2)

    def test_slice_raises_error_when_cutting_both_interval_ends(self):
        """Test slicing raises ValueError when both interval ends get cut"""
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phosphorylated]TIDE")

        with self.assertRaises(ValueError):
            annotation.slice(2, 4)

    def test_slice_preserves_complete_ambiguous_interval(self):
        """Test slicing preserves ambiguous interval when completely within slice"""
        annotation = pt.ProFormaAnnotation.parse("P(?EPTI)[Phosphorylated]DE")
        self.assertEqual(annotation.serialize(), "P(?EPTI)[Phosphorylated]DE")

        sliced = annotation.slice(0, 7)
        self.assertEqual(sliced.serialize(), "P(?EPTI)[Phosphorylated]DE")

    def test_slice_raises_error_when_cutting_ambiguous_interval(self):
        """Test slicing raises ValueError when ambiguous interval gets cut"""
        annotation = pt.ProFormaAnnotation.parse("P(?EPTI)[Phosphorylated]DE")

        with self.assertRaises(ValueError):
            annotation.slice(2, 4)

    def test_slice_preserves_labile_modifications(self):
        """Test slicing preserves labile modifications"""
        annotation = pt.ProFormaAnnotation.parse("{LabileMod}PEPTIDE")
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.serialize(), "{LabileMod}PEP")

    def test_slice_preserves_static_modifications(self):
        """Test slicing preserves static modifications"""
        annotation = pt.ProFormaAnnotation(
            sequence="PEPTIDE", static_mods=pt.Mod("57@C")
        )
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.serialize(), "<57@C>PEP")

    def test_slice_preserves_charge_state(self):
        """Test slicing preserves charge state"""
        annotation = pt.ProFormaAnnotation(sequence="PEPTIDE", charge=2)
        sliced = annotation.slice(0, 3)
        self.assertEqual(sliced.serialize(), "PEP/2")


if __name__ == "__main__":
    unittest.main()
