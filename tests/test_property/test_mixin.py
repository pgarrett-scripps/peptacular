import peptacular as pt


class TestSequencePropertyMixin:
    """Test the SequencePropertyMixin methods via pt.ProFormaAnnotation"""

    def test_hydrophobicity(self):
        """Test hydrophobicity calculation"""
        # Hydrophobic peptide
        annot = pt.ProFormaAnnotation(sequence="LIVM")
        assert isinstance(annot.prop.hydrophobicity, float)

        # Hydrophilic peptide
        annot = pt.ProFormaAnnotation(sequence="RKDE")
        assert isinstance(annot.prop.hydrophobicity, float)

    def test_pi(self):
        """Test isoelectric point calculation"""
        # Acidic peptide
        annot = pt.ProFormaAnnotation(sequence="DE")
        assert annot.prop.pi < 7.0

        # Basic peptide
        annot = pt.ProFormaAnnotation(sequence="RK")
        assert annot.prop.pi > 7.0

        # Neutral-ish
        annot = pt.ProFormaAnnotation(sequence="A")
        assert 5.0 < annot.prop.pi < 7.0

    def test_charge_at_ph(self):
        """Test charge calculation at specific pH"""
        annot = pt.ProFormaAnnotation(sequence="K")  # Lysine, basic

        # At acidic pH, should be positive
        assert annot.prop.charge_at_ph(pH=2.0) > 0

        # At basic pH, should be neutral or negative (deprotonated)
        assert annot.prop.charge_at_ph(pH=12.0) <= 0

    def test_aromaticity(self):
        """Test aromaticity calculation"""
        # All aromatic
        annot = pt.ProFormaAnnotation(sequence="WFY")
        assert isinstance(annot.prop.aromaticity, float)

        # None aromatic
        annot = pt.ProFormaAnnotation(sequence="A")
        assert isinstance(annot.prop.aromaticity, float)

        # Mixed
        annot = pt.ProFormaAnnotation(sequence="WA")
        assert isinstance(annot.prop.aromaticity, float)

    def test_secondary_structure(self):
        """Test secondary structure propensities"""
        annot = pt.ProFormaAnnotation(sequence="AAAA")

        # Just check that we get values between 0 and 1 (or reasonable percentages)
        # Note: The implementation might return raw scores or normalized values.
        # Based on mixin, it returns values from a scale.

        assert isinstance(annot.prop.alpha_helix_percent, float)
        assert isinstance(annot.prop.beta_sheet_percent, float)
        assert isinstance(annot.prop.beta_turn_percent, float)
        assert isinstance(annot.prop.coil_percent, float)

    def test_calc_property_generic(self):
        """Test generic calc_property method"""
        annot = pt.ProFormaAnnotation(sequence="A")
        # Using a dummy scale
        scale = {"A": 10.0}
        val = annot.prop.calc_property(scale=scale)
        assert val == 10.0

    def test_other_properties(self):
        """Test other properties exist and return floats"""
        annot = pt.ProFormaAnnotation(sequence="ACDEFGHIKLMNPQRSTVWY")

        assert isinstance(annot.prop.flexibility, float)
        assert isinstance(annot.prop.hydrophilicity, float)
        assert isinstance(annot.prop.surface_accessibility, float)
        assert isinstance(annot.prop.polarity, float)
        assert isinstance(annot.prop.mutability, float)
        assert isinstance(annot.prop.codons, float)
        assert isinstance(annot.prop.bulkiness, float)
        assert isinstance(annot.prop.recognition_factors, float)
        assert isinstance(annot.prop.transmembrane_tendency, float)
        assert isinstance(annot.prop.average_buried_area, float)
        assert isinstance(annot.prop.hplc, float)
        assert isinstance(annot.prop.refractivity, float)

    def test_sliding_window(self):
        """Test sliding window feature generation"""
        annot = pt.ProFormaAnnotation(sequence="AAAAA")
        scale = {"A": 1.0}
        # Window size 3, sequence length 5 -> 3 windows?
        # Depends on implementation details (padding etc).
        # generate_sliding_window_features usually returns a list of values.

        features = annot.prop.property_partitions(scale=scale, num_windows=3)
        assert isinstance(features, list)
        assert len(features) > 0
        assert all(isinstance(x, float) for x in features)
