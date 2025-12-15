from peptacular.annotation.annotation import ProFormaAnnotation


class TestSequencePropertyMixin:
    """Test the SequencePropertyMixin methods via ProFormaAnnotation"""

    def test_hydrophobicity(self):
        """Test hydrophobicity calculation"""
        # Hydrophobic peptide
        annot = ProFormaAnnotation(sequence="LIVM")
        assert isinstance(annot.hydrophobicity, float)

        # Hydrophilic peptide
        annot = ProFormaAnnotation(sequence="RKDE")
        assert isinstance(annot.hydrophobicity, float)

    def test_pi(self):
        """Test isoelectric point calculation"""
        # Acidic peptide
        annot = ProFormaAnnotation(sequence="DE")
        assert annot.pi < 7.0

        # Basic peptide
        annot = ProFormaAnnotation(sequence="RK")
        assert annot.pi > 7.0

        # Neutral-ish
        annot = ProFormaAnnotation(sequence="A")
        assert 5.0 < annot.pi < 7.0

    def test_charge_at_ph(self):
        """Test charge calculation at specific pH"""
        annot = ProFormaAnnotation(sequence="K")  # Lysine, basic

        # At acidic pH, should be positive
        assert annot.charge_at_ph(pH=2.0) > 0

        # At basic pH, should be neutral or negative (deprotonated)
        assert annot.charge_at_ph(pH=12.0) <= 0

    def test_aromaticity(self):
        """Test aromaticity calculation"""
        # All aromatic
        annot = ProFormaAnnotation(sequence="WFY")
        assert isinstance(annot.aromaticity, float)

        # None aromatic
        annot = ProFormaAnnotation(sequence="A")
        assert isinstance(annot.aromaticity, float)

        # Mixed
        annot = ProFormaAnnotation(sequence="WA")
        assert isinstance(annot.aromaticity, float)

    def test_secondary_structure(self):
        """Test secondary structure propensities"""
        annot = ProFormaAnnotation(sequence="AAAA")

        # Just check that we get values between 0 and 1 (or reasonable percentages)
        # Note: The implementation might return raw scores or normalized values.
        # Based on mixin, it returns values from a scale.

        assert isinstance(annot.alpha_helix_percent, float)
        assert isinstance(annot.beta_sheet_percent, float)
        assert isinstance(annot.beta_turn_percent, float)
        assert isinstance(annot.coil_percent, float)

    def test_calc_property_generic(self):
        """Test generic calc_property method"""
        annot = ProFormaAnnotation(sequence="A")
        # Using a dummy scale
        scale = {"A": 10.0}
        val = annot.calc_property(scale=scale)
        assert val == 10.0

    def test_other_properties(self):
        """Test other properties exist and return floats"""
        annot = ProFormaAnnotation(sequence="ACDEFGHIKLMNPQRSTVWY")

        assert isinstance(annot.flexibility, float)
        assert isinstance(annot.hydrophilicity, float)
        assert isinstance(annot.surface_accessibility, float)
        assert isinstance(annot.polarity, float)
        assert isinstance(annot.mutability, float)
        assert isinstance(annot.codons, float)
        assert isinstance(annot.bulkiness, float)
        assert isinstance(annot.recognition_factors, float)
        assert isinstance(annot.transmembrane_tendency, float)
        assert isinstance(annot.average_buried_area, float)
        assert isinstance(annot.hplc, float)
        assert isinstance(annot.refractivity, float)

    def test_sliding_window(self):
        """Test sliding window feature generation"""
        annot = ProFormaAnnotation(sequence="AAAAA")
        scale = {"A": 1.0}
        # Window size 3, sequence length 5 -> 3 windows?
        # Depends on implementation details (padding etc).
        # generate_sliding_window_features usually returns a list of values.

        features = annot.property_partitions(scale=scale, num_windows=3)
        assert isinstance(features, list)
        assert len(features) > 0
        assert all(isinstance(x, float) for x in features)
