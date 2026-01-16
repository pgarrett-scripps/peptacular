import unittest

from peptacular.annotation.frag import Fragment
from peptacular.fragment import IonType


class TestFragmentMzPAF(unittest.TestCase):
    def test_simple_peptide_to_mzpaf(self):
        frag = Fragment(
            ion_type=IonType.PRECURSOR,
            position=None,
            mass=1000.0,
            charge_state=1,
            monoisotopic=True,
            sequence="PEPTIDE",
        )
        mzpaf = frag.to_mzpaf(include_annotation=True).serialize()
        self.assertEqual(mzpaf, "p")

        mzpaf = frag.to_mzpaf(include_annotation=False).serialize()
        self.assertEqual(mzpaf, "p")

        mzpaf = frag.to_mzpaf(include_annotation=False, confidence=0.95).serialize()
        self.assertEqual(mzpaf, "p*0.95")

        mzpaf = frag.to_mzpaf(
            include_annotation=False,
            confidence=0.95,
            mass_error=0.1,
            mass_error_type="da",
        ).serialize()
        self.assertEqual(mzpaf, "p/0.1*0.95")

        mzpaf = frag.to_mzpaf(
            include_annotation=False,
            confidence=0.95,
            mass_error=0.1,
            mass_error_type="ppm",
        ).serialize()
        self.assertEqual(mzpaf, "p/0.1ppm*0.95")

    def test_fragment_with_loss_and_isotope(self):
        frag = Fragment(
            ion_type=IonType.B,
            position=3,
            mass=500.0,
            charge_state=2,
            monoisotopic=True,
            sequence="PEP",
            _losses={"H2O": 1},
            _isotopes={"13C": 2},
        )
        mzpaf = frag.to_mzpaf().serialize()
        self.assertEqual(mzpaf, "b3{PEP}Formula:H2O+2i13C^2")

    def test_fragment_with_charge_adduct(self):
        frag = Fragment(
            ion_type=IonType.Y,
            position=5,
            mass=700.0,
            charge_state=1,
            monoisotopic=True,
            sequence="PEPTI",
            _charge_adducts=("H",),
        )
        mzpaf = frag.to_mzpaf().serialize()
        self.assertEqual(mzpaf, "y5{PEPTI}[M+H]")

    def test_fragment_with_multiple_losses(self):
        frag = Fragment(
            ion_type=IonType.Y,
            position=7,
            mass=900.0,
            charge_state=1,
            monoisotopic=True,
            sequence="PEPTIDER",
            _losses={"H2O": 1, "NH3": 1},
        )
        mzpaf = frag.to_mzpaf().serialize()
        self.assertEqual(mzpaf, "y7{PEPTIDER}-H2O-NH3")

    def test_fragment_with_adduct_and_charge(self):
        frag = Fragment(
            ion_type=IonType.B,
            position=2,
            mass=600.0,
            charge_state=2,
            monoisotopic=True,
            sequence="PE",
            _charge_adducts=("Na", "H"),
        )
        mzpaf = frag.to_mzpaf()
        self.assertEqual(mzpaf.charge, 2)
        self.assertIsInstance(mzpaf.adducts, list)
        self.assertTrue(any("Na" in str(a) for a in mzpaf.adducts))
        self.assertTrue(any("H" in str(a) for a in mzpaf.adducts))
        # Check serialization
        self.assertTrue("[M+Na+H]^2" in mzpaf.serialize())

    def test_fragment_with_isotope_and_charge(self):
        frag = Fragment(
            ion_type=IonType.X,
            position=4,
            mass=700.0,
            charge_state=3,
            monoisotopic=True,
            sequence="PEPT",
            _isotopes={"13C": 3},
        )
        mzpaf = frag.to_mzpaf()
        self.assertEqual(mzpaf.charge, 3)
        self.assertIsInstance(mzpaf.isotope, list)
        self.assertTrue(any(getattr(i, "count", None) == 3 for i in mzpaf.isotope))
        # Check serialization
        self.assertTrue("+3i^3" in mzpaf.serialize())

    def test_fragment_with_multiple_adducts_and_losses(self):
        frag = Fragment(
            ion_type=IonType.Y,
            position=8,
            mass=950.0,
            charge_state=2,
            monoisotopic=True,
            sequence="PEPTIDERK",
            _charge_adducts=("Na", "K"),
            _losses={"H2O": 1, "SO3": 1},
        )
        mzpaf = frag.to_mzpaf().serialize()
        self.assertEqual("y8{PEPTIDERK}-H2O-SO3[M+Na+K]^2", mzpaf)

    def test_fragment_with_no_sequence(self):
        frag = Fragment(
            ion_type=IonType.B,
            position=1,
            mass=400.0,
            charge_state=1,
            monoisotopic=True,
            sequence=None,
        )
        mzpaf = frag.to_mzpaf()
        self.assertEqual(mzpaf.charge, 1)
        # Check serialization
        self.assertTrue(mzpaf.serialize().startswith("b1"))

    def test_fragment_with_unknown_ion_type(self):
        frag = Fragment(
            ion_type=None,
            position=None,
            mass=500.0,
            charge_state=1,
            monoisotopic=True,
            sequence="PEPTIDE",
        )
        mzpaf = frag.to_mzpaf()
        self.assertEqual(mzpaf.charge, 1)
        # Check serialization
        self.assertTrue(mzpaf.serialize().startswith("?"))

    def test_fragment_with_immonium_ion(self):
        frag = Fragment(
            ion_type=IonType.IMMONIUM,
            position="Y",
            mass=136.0,
            charge_state=1,
            monoisotopic=True,
            sequence=None,
        )
        mzpaf = frag.to_mzpaf()
        self.assertEqual(mzpaf.ion_type.__class__.__name__, "ImmoniumIon")
        self.assertEqual(mzpaf.charge, 1)
        # Check serialization
        self.assertTrue(mzpaf.serialize().startswith("IY"))

    def test_fragment_with_internal_fragment(self):
        frag = Fragment(
            ion_type=IonType.INTERNAL,
            position=(2, 5),
            mass=600.0,
            charge_state=1,
            monoisotopic=True,
            sequence="EPTI",
        )
        mzpaf = frag.to_mzpaf()
        self.assertEqual(mzpaf.ion_type.__class__.__name__, "InternalFragment")
        self.assertEqual(mzpaf.charge, 1)
        # Check serialization
        self.assertTrue(mzpaf.serialize().startswith("m2:5"))


if __name__ == "__main__":
    unittest.main()
