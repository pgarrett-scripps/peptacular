import unittest

import peptacular as pt


class TestFragment(unittest.TestCase):
    def test_get_losses_single_match_max_losses_1(self):
        """Test get_losses with single amino acid and max_losses=1."""
        result = pt.get_losses("AA", {"A": [-10, -5]}, 1)
        expected: set[float] = {0.0, -5, -10}
        self.assertEqual(result, expected)

    def test_get_losses_single_match_max_losses_2(self):
        """Test get_losses with single amino acid and max_losses=2."""
        result = pt.get_losses("AA", {"A": [-10, -5]}, 2)
        expected: set[float] = {0.0, -15, -10, -5, -20}
        self.assertEqual(result, expected)

    def test_get_losses_single_match_max_losses_3(self):
        """Test get_losses with single amino acid and max_losses=3."""
        result = pt.get_losses("AA", {"A": [-10, -5]}, 3)
        expected: set[float] = {0.0, -15, -10, -25, -5, -20}
        self.assertEqual(result, expected)

    def test_fragment_modified_sequence_count(self):
        """Test fragment count with modified sequence."""
        result: list[float] = list(pt.fragment("[1.0]-P[2.0]E[3.0]-[4.0]", "y", 1))
        self.assertEqual(len(result), 2)

    def test_fragment_b_ions_mz_charge_1(self):
        """Test B-ion m/z values with charge 1."""
        result: list[float] = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["b"],
                charges=1,
                return_type="mz",
                precision=3,
            )
        )
        expected = [459.209, 330.166, 215.139, 102.055]
        self.assertEqual(result, expected)

    def test_fragment_b_ions_mz_charge_2(self):
        """Test B-ion m/z values with charge 2."""
        result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["b"],
                charges=2,
                return_type="mz",
                precision=3,
            )
        )
        expected = [230.108, 165.587, 108.073, 51.531]
        self.assertEqual(result, expected)

    def test_fragment_immonium_ions_mz(self):
        """Test immonium ions m/z values."""
        result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["i"],
                charges=1,
                return_type="mz",
                precision=3,
            )
        )
        expected = [74.06, 86.096, 88.039, 102.055]
        self.assertEqual(result, expected)

    def test_fragment_internal_ions_mz_label(self):
        """Test internal fragment ions with mz-label return type."""
        result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["by"],
                charges=1,
                return_type="mz-label",
                precision=3,
            )
        )
        expected = [(114.091, "+by1-2"), (229.118, "+by1-3"), (116.034, "+by2-3")]
        self.assertEqual(result, expected)

    def test_fragment_regex_losses(self):
        """Test fragment with regex-based losses."""
        result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["b"],
                charges=2,
                return_type="mz",
                precision=3,
                losses={"E": [-10]},
            )
        )
        expected = [230.108, 225.108, 165.587, 108.073, 51.531]
        self.assertEqual(result, expected)

    def test_fragment_multiple_losses_max_1(self):
        """Test fragment with multiple losses and max_losses=1."""
        losses = {"A": [-10, -5]}
        result = list(
            pt.fragment("AA", ["b"], 1, return_type="mz", precision=3, losses=losses)
        )
        expected = [143.082, 138.082, 133.082, 72.044, 67.044, 62.044]
        self.assertEqual(result, expected)

    def test_fragment_multiple_losses_label_max_1(self):
        """Test fragment with multiple losses returning labels and max_losses=1."""
        losses = {"A": [-10, -5]}
        result = list(
            pt.fragment(
                "AA", ["b"], 1, return_type="mz-label", precision=3, losses=losses
            )
        )
        lables = [label for _, label in result]
        expected = ["+b2", "+b2(-5)", "+b2(-10)", "+b1", "+b1(-5)", "+b1(-10)"]
        self.assertEqual(lables, expected)

    def test_fragment_multiple_losses_max_2(self):
        """Test fragment with multiple losses and max_losses=2."""
        losses = {"A": [-10, -5]}
        result = list(
            pt.fragment(
                "AA",
                ["b"],
                1,
                return_type="mz",
                precision=3,
                losses=losses,
                max_losses=2,
            )
        )
        expected = [
            143.082,
            123.082,
            128.082,
            133.082,
            138.082,
            72.044,
            57.044,
            67.044,
            62.044,
        ]
        self.assertEqual(result, expected)

    def test_fragment_multiple_losses_label_max_2(self):
        """Test fragment with multiple losses returning labels and max_losses=2."""
        losses = {"A": [-10, -5]}
        result = list(
            pt.fragment(
                "AA",
                ["b"],
                1,
                return_type="mz-label",
                precision=3,
                losses=losses,
                max_losses=2,
            )
        )
        labels = [label for _, label in result]
        expected = [
            "+b2",
            "+b2(-20)",
            "+b2(-15)",
            "+b2(-10)",
            "+b2(-5)",
            "+b1",
            "+b1(-15)",
            "+b1(-5)",
            "+b1(-10)",
        ]
        self.assertEqual(labels, expected)

    def test_fragment_water_ammonia_loss_mz(self):
        """Test fragment with water and ammonia losses m/z values."""
        result = list(
            pt.fragment(
                "AQE",
                ["b"],
                1,
                return_type="mz",
                precision=3,
                water_loss=True,
                ammonia_loss=True,
            )
        )
        expected = [329.146, 311.135, 312.119, 200.103, 183.076, 72.044]
        self.assertEqual(result, expected)

    def test_fragment_water_ammonia_loss_label(self):
        """Test fragment with water and ammonia losses returning labels."""
        result = list(
            pt.fragment(
                "AQE",
                ["b"],
                1,
                return_type="mz-label",
                precision=3,
                water_loss=True,
                ammonia_loss=True,
            )
        )
        labels = [label for _, label in result]
        expected = [
            "+b3",
            "+b3(-18.011)",
            "+b3(-17.027)",
            "+b2",
            "+b2(-17.027)",
            "+b1",
        ]
        self.assertEqual(labels, expected)

    def test_fragment_multiple_ion_types(self):
        """Test fragment with multiple ion types."""
        result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["b", "y"],
                charges=1,
                return_type="mz",
                precision=3,
            )
        )
        # Should contain both b and y ions
        self.assertTrue(len(result) > 4)  # More than just b or y alone

    def test_fragment_multiple_charges(self):
        """Test fragment with multiple charge states."""
        result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["b"],
                charges=[1, 2],
                return_type="mz",
                precision=3,
            )
        )
        # Should contain fragments for both charge states
        self.assertEqual(len(result), 8)  # 4 fragments * 2 charges

    def test_fragment_isotopes(self):
        """Test fragment with isotope offsets."""
        result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["b"],
                charges=1,
                isotopes=[0, 1],
                return_type="mz",
                precision=3,
            )
        )
        # Should contain fragments for both isotope states
        self.assertEqual(len(result), 8)  # 4 fragments * 2 isotopes

    def test_fragment_empty_sequence(self):
        """Test fragment with empty sequence."""
        result = list(
            pt.fragment(sequence="", ion_types=["b"], charges=1, return_type="mz")
        )
        self.assertEqual(result, [])

    def test_fragment_single_amino_acid(self):
        """Test fragment with single amino acid."""
        result = list(
            pt.fragment(
                sequence="A", ion_types=["b"], charges=1, return_type="mz", precision=3
            )
        )
        self.assertEqual(len(result), 1)
        self.assertIsInstance(result[0], float)

    def test_fragment_average_mass(self):
        """Test fragment with average masses instead of monoisotopic."""
        mono_result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["b"],
                charges=1,
                monoisotopic=True,
                return_type="mz",
                precision=3,
            )
        )
        avg_result = list(
            pt.fragment(
                sequence="TIDE",
                ion_types=["b"],
                charges=1,
                monoisotopic=False,
                return_type="mz",
                precision=3,
            )
        )
        # Results should be different
        self.assertNotEqual(mono_result, avg_result)
        self.assertEqual(len(mono_result), len(avg_result))

    def test_fragment_complex_modified_sequence(self):
        """Test fragment with complex ProForma modifications."""
        modified_seq = "[Acetyl]-M[Oxidation]PEPTIDE[Phospho]K-[Amidated]"
        result = list(
            pt.fragment(
                sequence=modified_seq,
                ion_types=["y"],
                charges=1,
                return_type="mz",
                precision=3,
            )
        )
        # Should successfully fragment modified sequence
        self.assertTrue(len(result) > 0)
        for mz in result:
            self.assertIsInstance(mz, float)

    def test_fragment_water_loss_only(self):
        """Test fragment with only water loss enabled."""
        result = list(
            pt.fragment(
                "STED", ["b"], 1, return_type="mz", precision=3, water_loss=True
            )
        )
        # Should include fragments with and without water loss
        self.assertTrue(len(result) > 4)  # More than base fragments

    def test_fragment_ammonia_loss_only(self):
        """Test fragment with only ammonia loss enabled."""
        result = list(
            pt.fragment(
                "RKNQ", ["b"], 1, return_type="mz", precision=3, ammonia_loss=True
            )
        )
        # Should include fragments with and without ammonia loss
        self.assertTrue(len(result) > 4)  # More than base fragments

    def test_fragment_precision_none(self):
        """Test fragment with no precision specified."""
        result = list(
            pt.fragment(sequence="TIDE", ion_types=["b"], charges=1, return_type="mz")
        )
        self.assertTrue(len(result) > 0)
        for mz in result:
            self.assertIsInstance(mz, float)

    def test_fragment_complex_losses_mapping(self):
        """Test fragment with complex losses mapping."""
        losses = {"K": [-10.0, -5.0], "R": [-15.0], "[ST]": [-20.0]}
        result = list(
            pt.fragment(
                sequence="KSTR",
                ion_types=["b"],
                charges=1,
                return_type="mz",
                precision=3,
                losses=losses,
            )
        )
        # Should include base fragments plus various losses
        self.assertTrue(len(result) > 4)

    def test_fragment_y_ions_basic(self):
        """Test basic Y-ion fragmentation."""
        result = list(
            pt.fragment(
                sequence="PEPTIDE",
                ion_types=["y"],
                charges=1,
                return_type="mz",
                precision=3,
            )
        )
        self.assertEqual(len(result), 7)  # 7 y-ions for 7-residue peptide
        for mz in result:
            self.assertIsInstance(mz, float)
            self.assertGreater(mz, 0)

    def test_fragment_c_ions(self):
        """Test C-ion fragmentation."""
        result = list(
            pt.fragment(
                sequence="PEPTIDE",
                ion_types=["c"],
                charges=1,
                return_type="mz",
                precision=3,
            )
        )
        self.assertTrue(len(result) > 0)
        for mz in result:
            self.assertIsInstance(mz, float)

    def test_fragment_z_ions(self):
        """Test Z-ion fragmentation."""
        result = list(
            pt.fragment(
                sequence="PEPTIDE",
                ion_types=["z"],
                charges=1,
                return_type="mz",
                precision=3,
            )
        )
        self.assertTrue(len(result) > 0)
        for mz in result:
            self.assertIsInstance(mz, float)


if __name__ == "__main__":
    unittest.main()
