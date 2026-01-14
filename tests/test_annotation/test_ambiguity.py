import unittest

import peptacular as pt


class TestCondenseAmbiguityToXNotation(unittest.TestCase):
    def test_basic_ambiguous_interval(self):
        """Test condensing basic ambiguous interval to X notation"""
        annotation = pt.ProFormaAnnotation.parse("PEP(?TIDE)[Phospho]")
        condensed = annotation.condense_ambiguity_to_xnotation()

        self.assertEqual(condensed.sequence, "PEPX")  # PEPX
        self.assertAlmostEqual(annotation.mass(), condensed.mass(), places=2)

    def test_inplace_modification(self):
        """Test in-place modification of annotation"""
        annotation = pt.ProFormaAnnotation.parse("PEP(?TIDE)[Phospho]")
        condensed = annotation.condense_ambiguity_to_xnotation(inplace=True)

        self.assertEqual(condensed.sequence, "PEPX")  # PEPX
        self.assertEqual(id(condensed), id(annotation))

    def test_no_intervals_unchanged(self):
        """Test that annotation without intervals is unchanged"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE[Phospho]")
        condensed = annotation.condense_ambiguity_to_xnotation()

        self.assertEqual(condensed.serialize(), "PEPTIDE[Phospho]")

    def test_multiple_ambiguous_intervals(self):
        """Test condensing multiple ambiguous intervals"""
        annotation = pt.ProFormaAnnotation.parse("(?PEP)[Acetyl]TI(?DE)[Phospho]")
        condensed = annotation.condense_ambiguity_to_xnotation()

        self.assertEqual(condensed.sequence, "XTIX")  # PEPX
        self.assertAlmostEqual(annotation.mass(), condensed.mass(), places=2)

    def test_preserves_non_interval_modifications(self):
        """Test that non-interval modifications are preserved"""
        annotation = pt.ProFormaAnnotation.parse(
            "{Oxidation}[Acetyl]-PEP(?TIDE)[Phospho]/2"
        )
        condensed = annotation.condense_ambiguity_to_xnotation()

        self.assertEqual(condensed.sequence, "PEPX")  # PEPX
        self.assertAlmostEqual(annotation.mass(), condensed.mass(), places=2)


class TestCondenseModsToIntervals(unittest.TestCase):
    def test_basic_modification_condensation(self):
        """Test condensing internal modification into existing interval"""
        annotation = pt.ProFormaAnnotation.parse("[2]-PEPT(I[Phospho]D)E-[1]")
        condensed = annotation.condense_mods_to_intervals()

        # T[Methyl] should move into the interval (TI)
        self.assertEqual(condensed.serialize(), "[2]-PEPT(ID)[Phospho]E-[1]")
        self.assertEqual(id(annotation), id(condensed))

    def test_inplace_false(self):
        """Test that inplace=False returns a copy"""
        annotation = pt.ProFormaAnnotation.parse("PEPT(I[Phospho]D)E")
        condensed = annotation.condense_mods_to_intervals(inplace=False)

        # T[Methyl] should move into the interval (TI)
        self.assertEqual(condensed.serialize(), "PEPT(ID)[Phospho]E")
        self.assertNotEqual(id(annotation), id(condensed))

    def test_no_intervals_unchanged(self):
        """Test that annotation without intervals is unchanged"""
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE[Phospho]")
        condensed = annotation.condense_mods_to_intervals()

        self.assertEqual(condensed.serialize(), "PEPTIDE[Phospho]")

    def test_no_internal_mods_unchanged(self):
        """Test that annotation without internal mods is unchanged"""
        annotation = pt.ProFormaAnnotation.parse("PEP(TI)[Phospho]DE")
        condensed = annotation.condense_mods_to_intervals()

        self.assertEqual(condensed.serialize(), "PEP(TI)[Phospho]DE")

    def test_modification_outside_interval_unchanged(self):
        """Test that internal mods outside intervals remain unchanged"""
        annotation = pt.ProFormaAnnotation.parse("PE[Methyl]P(TI)[Phospho]DE")
        condensed = annotation.condense_mods_to_intervals()

        self.assertEqual(condensed.serialize(), "PE[Methyl]P(TI)[Phospho]DE")

    def test_multiple_intervals_and_mods(self):
        """Test condensing multiple mods into multiple intervals"""
        annotation = pt.ProFormaAnnotation.parse("PEPT(I[Phospho]D)EPEPT(I[Phospho]D)E")
        condensed = annotation.condense_mods_to_intervals()

        # T[Methyl] should move into the interval (TI)
        self.assertEqual(condensed.serialize(), "PEPT(ID)[Phospho]EPEPT(ID)[Phospho]E")

    def test_multiple_mods_same_interval(self):
        """Test condensing multiple mods into same interval"""
        annotation = pt.ProFormaAnnotation.parse("PEPT[Phospho]I[Methyl]DE")
        # Assuming (TI) interval exists at positions 3-4
        annotation = pt.ProFormaAnnotation.parse("PEP(T[Phospho]I[Methyl])DE")
        condensed = annotation.condense_mods_to_intervals()

        self.assertEqual(condensed.serialize(), "PEP(TI)[Phospho][Methyl]DE")


if __name__ == "__main__":
    unittest.main()
