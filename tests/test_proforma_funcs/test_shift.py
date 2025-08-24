# filepath: /workspaces/peptacular/tests/test_proforma_funcs/test_shift.py
import unittest

import peptacular as pt


class TestShift(unittest.TestCase):
    
    def test_basic_shift(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shifted_annotation = annotation.shift(2)
        self.assertEqual(shifted_annotation.sequence, "PTIDEPE")
        self.assertEqual(annotation.sequence, "PEPTIDE")  # Original unchanged

    def test_basic_shift_inplace(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shifted_annotation = annotation.shift(2, inplace=True)
        self.assertEqual(shifted_annotation.sequence, "PTIDEPE")
        self.assertEqual(annotation.sequence, "PTIDEPE")  # Original changed too

    def test_shift_zero(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shifted_annotation = annotation.shift(0)
        self.assertEqual(shifted_annotation.sequence, "PEPTIDE")
        self.assertEqual(annotation.sequence, "PEPTIDE")

    def test_shift_full_length(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shifted_annotation = annotation.shift(7)  # length of PEPTIDE
        self.assertEqual(shifted_annotation.sequence, "PEPTIDE")

    def test_shift_negative(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shifted_annotation = annotation.shift(-2)
        self.assertEqual(shifted_annotation.serialize(), "DEPEPTI")


    def test_shift_large_number(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE")
        shifted_annotation = annotation.shift(9)
        self.assertEqual(shifted_annotation.serialize(), "PTIDEPE")

    """
    TESTS FOR: shifting with internal modifications
    """
    def test_shift_with_internal_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho]PTI[Methyl]DE")
        shifted_annotation = annotation.shift(2)
        self.assertEqual(shifted_annotation.serialize(), "PTI[Methyl]DEPE[Phospho]")

    def test_shift_with_multiple_internal_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE[Phospho][Methyl]PTIDE")
        shifted_annotation = annotation.shift(1)
        self.assertEqual(shifted_annotation.serialize(), "E[Phospho][Methyl]PTIDEP")

    """
    TESTS FOR: shifting with intervals
    """
    def test_shift_with_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")
        shifted_annotation = annotation.shift(2)
        self.assertEqual(shifted_annotation.serialize(), "(P)[#a1]TIDEP(E)[Phospho#a1]")

    def test_shift_with_ambiguous_intervals(self):
        annotation = pt.ProFormaAnnotation.parse("P(?EPT)[Phospho]IDE")        
        shifted_annotation = annotation.shift(1)
        self.assertEqual(shifted_annotation.serialize(), "(?EPT)[Phospho]IDEP")

    def test_shift_with_interval_no_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PE(?PTI)DE")        
        shifted_annotation = annotation.shift(3)
        self.assertEqual(shifted_annotation.serialize(), "(?TI)DEPE(?P)")

    """
    TESTS FOR: shifting with terminal modifications
    """
    def test_shift_with_nterm_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE")
        shifted_annotation = annotation.shift(2)
        self.assertEqual(shifted_annotation.serialize(), "[Acetyl]-PTIDEPE")

    def test_shift_with_cterm_mods(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE-[Amide]")
        shifted_annotation = annotation.shift(3)
        self.assertEqual(shifted_annotation.serialize(), "TIDEPEP-[Amide]")

    def test_shift_with_both_term_mods(self):
        annotation = pt.ProFormaAnnotation.parse("[Acetyl]-PEPTIDE-[Amide]")
        shifted_annotation = annotation.shift(1)
        self.assertEqual(shifted_annotation.serialize(), "[Acetyl]-EPTIDEP-[Amide]")

    """
    TESTS FOR: shifting with other modification types
    """
    def test_shift_with_labile_mods(self):
        annotation = pt.ProFormaAnnotation.parse("{Glycan}PEPTIDE")
        shifted_annotation = annotation.shift(2)
        self.assertEqual(shifted_annotation.serialize(), "{Glycan}PTIDEPE")

    def test_shift_with_static_mods(self):
        # Test shifting with static modifications (should be preserved)
        annotation = pt.ProFormaAnnotation.parse("<57@C>PEPTIDE")
        shifted_annotation = annotation.shift(2)
        self.assertEqual(shifted_annotation.serialize(), "<57@C>PTIDEPE")

    def test_shift_with_charge(self):
        annotation = pt.ProFormaAnnotation.parse("PEPTIDE/2")
        shifted_annotation = annotation.shift(2)
        self.assertEqual(shifted_annotation.serialize(), "PTIDEPE/2")

    def test_shift_with_unknown_mod(self):
        annotation = pt.ProFormaAnnotation.parse("[Unknown]?PEPTIDE")
        shifted_annotation = annotation.shift(2)
        self.assertEqual(shifted_annotation.serialize(), "[Unknown]?PTIDEPE")

    def test_shift_single_amino_acid(self):
        annotation = pt.ProFormaAnnotation(sequence="A")
        shifted_annotation = annotation.shift(5)
        self.assertEqual(shifted_annotation.sequence, "A")

    def test_shift_single_amino_acid2(self):
        annotation = pt.ProFormaAnnotation(sequence="")
        shifted_annotation = annotation.shift(5)
        self.assertEqual(shifted_annotation.sequence, "")

    """
    TESTS FOR: shift_intervals parameter
    """
    def test_shift_with_shift_intervals_true(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")        
        shifted_annotation = annotation.shift(2, shift_intervals=True)
        self.assertEqual(shifted_annotation.serialize(), "(P)[#a1]TIDEP(E)[Phospho#a1]")

    def test_shift_with_shift_intervals_false(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")        
        shifted_annotation = annotation.shift(2, shift_intervals=False)
        self.assertEqual(shifted_annotation.serialize(), "P(TI)[Phospho]DEPE")

    """
    TESTS FOR: slice_intervals parameter
    """
    def test_shift_with_slice_intervals_true(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")        
        shifted_annotation = annotation.shift(2, slice_intervals=True)
        self.assertEqual(shifted_annotation.serialize(), "(P)[#a1]TIDEP(E)[Phospho#a1]")

    def test_shift_with_slice_intervals_true_multi_mods(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho][Oxidation]TIDE")        
        shifted_annotation = annotation.shift(2, slice_intervals=True)
        self.assertEqual(shifted_annotation.serialize(), "(P)[#a1][#a2]TIDEP(E)[Phospho#a1][Oxidation#a2]")

    def test_shift_with_slice_intervals_true_ambiguous(self):
        annotation = pt.ProFormaAnnotation.parse("P(?EP)TIDE")        
        shifted_annotation = annotation.shift(2, slice_intervals=True)
        self.assertEqual(shifted_annotation.serialize(), "(?P)TIDEP(?E)")

    def test_shift_with_slice_intervals_false_raises_error(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")        
        with self.assertRaises(ValueError) as _:
            annotation.shift(2, slice_intervals=False)

    def test_shift_slice_intervals_false_no_error_when_no_wrap(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")        
        shifted_annotation = annotation.shift(1, slice_intervals=False)
        self.assertEqual(shifted_annotation.serialize(), "(EP)[Phospho]TIDEP")

    def test_shift_slice_intervals_false_no_error_when_over_wrap(self):
        annotation = pt.ProFormaAnnotation.parse("P(EP)[Phospho]TIDE")        
        shifted_annotation = annotation.shift(3, slice_intervals=False)
        self.assertEqual(shifted_annotation.serialize(), "TIDEP(EP)[Phospho]")
