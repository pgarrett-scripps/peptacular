from collections import Counter
import unittest
import peptacular as pt


class TestPopMods(unittest.TestCase):
    def test_pop_labile_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        popped_mods = Counter(map(str, annotation.pop_labile_mods()))

        self.assertIsNotNone(popped_mods)
        assert popped_mods is not None  # Type guard
        self.assertEqual(len(popped_mods), 1)
        self.assertFalse(annotation.has_labile_mods)
        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_labile_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_labile_mods()
        assert len(popped_mods) == 0

    def test_pop_labile_mods_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}{HexNAc}PEPTIDE")
        popped_mods = Counter(map(str, annotation.pop_labile_mods()))

        self.assertEqual(len(popped_mods), 2)
        self.assertFalse(annotation.has_labile_mods)

    def test_pop_unknown_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        popped_mods = Counter(map(str, annotation.pop_unknown_mods()))
        self.assertEqual(len(popped_mods), 1)
        self.assertFalse(annotation.has_unknown_mods)

    def test_pop_unknown_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_unknown_mods()
        self.assertEqual(len(popped_mods), 0)

    def test_pop_nterm_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = Counter(map(str, annotation.pop_nterm_mods()))

        self.assertEqual(len(popped_mods), 1)
        self.assertFalse(annotation.has_nterm_mods)

    def test_pop_nterm_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_nterm_mods()
        assert len(popped_mods) == 0

    def test_pop_nterm_mods_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl][Phospho]-PEPTIDE")
        popped_mods = Counter(map(str, annotation.pop_nterm_mods()))

        self.assertEqual(len(popped_mods), 2)
        self.assertFalse(annotation.has_nterm_mods)

    def test_pop_cterm_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        popped_mods = Counter(map(str, annotation.pop_cterm_mods()))

        self.assertEqual(len(popped_mods), 1)
        self.assertFalse(annotation.has_cterm_mods)

    def test_pop_cterm_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_cterm_mods()
        assert len(popped_mods) == 0

    def test_pop_internal_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        popped_mods = annotation.pop_internal_mods()
        popped_mods = {
            pos: Counter(map(str, mods)) for pos, mods in popped_mods.items()
        }

        self.assertEqual(len(popped_mods), 2)

    def test_pop_internal_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_internal_mods()
        self.assertEqual(popped_mods, {})

    def test_pop_intervals_with_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        popped_intervals = annotation.pop_intervals()

        assert popped_intervals is not None  # Type guard
        self.assertEqual(len(popped_intervals), 1)
        self.assertEqual(popped_intervals[0].start, 1)
        self.assertEqual(popped_intervals[0].end, 3)

    def test_pop_intervals_no_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_intervals = annotation.pop_intervals()
        self.assertEqual(popped_intervals, [])

    def test_pop_charge_with_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        popped_charge = annotation.pop_charge()

        self.assertEqual(popped_charge, 2)
        self.assertFalse(annotation.has_charge)

    def test_pop_charge_no_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_charge = annotation.pop_charge()
        self.assertEqual(popped_charge, None)

    def test_pop_charge_negative(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/-2")
        popped_charge = annotation.pop_charge()
        self.assertEqual(popped_charge, -2)

    def test_pop_charge_adducts_with_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/[Na:z+1]")
        popped_adducts = annotation.pop_charge()

        assert isinstance(popped_adducts, pt.Mods)
        self.assertEqual(len(popped_adducts), 1)

    def test_pop_charge_adducts_no_adducts(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        popped_adducts = annotation.pop_charge()
        assert isinstance(popped_adducts, int)
        self.assertEqual(popped_adducts, 2)

    def test_pop_charge_adducts_multiple(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/[Na:z+1,K:z+1]")
        popped_adducts = annotation.pop_charge()
        assert isinstance(popped_adducts, pt.Mods)
        self.assertEqual(len(popped_adducts), 2)

    def test_pop_isotope_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("<15N>PEPTIDE")
        popped_mods = Counter(map(str, annotation.pop_isotope_mods()))

        self.assertEqual(len(popped_mods), 1)

    def test_pop_isotope_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_isotope_mods()
        self.assertEqual(len(popped_mods), 0)

    def test_pop_static_mods_with_mods(self):
        annotation = pt.ProFormaAnnotation.parse("<[57]@C>PEPTIDE")
        popped_mods = Counter(map(str, annotation.pop_static_mods()))

        self.assertEqual(len(popped_mods), 1)

    def test_pop_static_mods_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        popped_mods = annotation.pop_static_mods()
        assert len(popped_mods) == 0

    def test_pop_mod_by_type_enum(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        popped_mods = annotation.pop_mods(mod_types=pt.ModType.NTERM)
        assert isinstance(popped_mods, dict)
        self.assertEqual(len(popped_mods), 1)
        self.assertIn(pt.ModType.NTERM, popped_mods)
        nterm_mods = popped_mods[pt.ModType.NTERM]
        assert isinstance(nterm_mods, pt.Mods)
        self.assertEqual(len(nterm_mods), 1)

    def test_pop_mods_all_default(self):
        annotation = pt.ProFormaAnnotation.parse(
            "<[57]@C>{Glycan}[Acetyl]-PE[Phospho]PTIDE-[Amide]/2"
        )
        popped_mods = annotation.pop_mods()

        expected_keys = ["labile", "static", "nterm", "internal", "cterm", "charge"]
        for key in expected_keys:
            self.assertIn(key, popped_mods)

        self.assertEqual(annotation.serialize(), "PEPTIDE")

    def test_pop_mods_specific_type(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mod_types="nterm")

        self.assertIn("nterm", popped_mods)
        self.assertEqual(len(popped_mods["nterm"]), 1)  # tuple length
        self.assertTrue(annotation.has_internal_mods)
        self.assertTrue(annotation.has_cterm_mods)

    def test_pop_mods_multiple_types(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mod_types=["nterm", "cterm"])

        self.assertIn("nterm", popped_mods)
        self.assertIn("cterm", popped_mods)
        self.assertTrue(annotation.has_internal_mods)

    def test_pop_mods_enum_type(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mod_types=pt.ModType.INTERNAL)

        self.assertIn("internal", popped_mods)
        self.assertTrue(annotation.has_nterm_mods)
        self.assertTrue(annotation.has_cterm_mods)

    def test_pop_mods_mixed_types(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PE[Phospho]PTIDE-[Amide]/2")
        popped_mods = annotation.pop_mods(mod_types=["nterm", pt.ModType.CHARGE])

        self.assertIn("nterm", popped_mods)
        self.assertIn("charge", popped_mods)
        self.assertEqual(popped_mods["charge"], 2)

    def test_pop_mods_empty_annotation(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")

        popped_mods = annotation.pop_mods()
        for key, value in popped_mods.items():
            if key == "charge":
                self.assertEqual(value, None)
            elif key == "interval":  # Note: "interval" not "intervals"
                self.assertEqual(len(value), 0)
            elif key == "internal":
                # internal mods returns {} when empty
                self.assertEqual(len(value), 0)
            else:
                # All other mod types return empty tuples
                self.assertEqual(len(value), 0)

    def test_pop_mods_empty_list(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        popped_mods = annotation.pop_mods(mod_types=[])

        self.assertEqual(annotation.serialize(), "[Acetyl]-PEPTIDE-[Amide]")
        self.assertEqual(popped_mods, {})

    def test_multiple_pops_return_none(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")

        first_pop = annotation.pop_nterm_mods()
        self.assertIsNotNone(first_pop)
        self.assertEqual(len(first_pop), 1)  # tuple length

        second_pop = annotation.pop_nterm_mods()
        assert len(second_pop) == 0

    def test_pop_mods_with_multipliers(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTIDE")
        popped_mods = annotation.pop_internal_mods()
        popped_mods = {
            pos: Counter(map(str, mods)) for pos, mods in popped_mods.items()
        }

        self.assertEqual(len(popped_mods), 1)


if __name__ == "__main__":
    unittest.main()
