"""
Tests using real-world ProForma modification examples.
"""

import peptacular as pt


class TestRealWorldExamples:
    """Tests using real-world ProForma modification examples"""

    def test_oxidation_variants(self):
        """Test different ways to represent oxidation"""
        # By name
        result1 = pt.parse_modification_tag("Oxidation")
        assert isinstance(result1, pt.TagName)

        # By UNIMOD
        result2 = pt.parse_modification_tag("UNIMOD:35")
        assert isinstance(result2, pt.TagAccession)

        # By mass
        result3 = pt.parse_modification_tag("+15.995")
        assert isinstance(result3, pt.TagMass)

        # By formula
        result4 = pt.parse_modification_tag("Formula:O")
        assert isinstance(result4, pt.ChargedFormula)

    def test_phosphorylation_variants(self):
        """Test different ways to represent phosphorylation"""
        # By name
        result1 = pt.parse_modification_tag("Phospho")
        assert isinstance(result1, pt.TagName)

        # By PSI-MOD
        result2 = pt.parse_modification_tag("MOD:00046")
        assert isinstance(result2, pt.TagAccession)

        # By mass
        result3 = pt.parse_modification_tag("+79.966")
        assert isinstance(result3, pt.TagMass)

    def test_crosslinker_modification(self):
        """Test cross-linking modification"""
        result = pt.parse_modification_tag("XLMOD:02001")
        assert isinstance(result, pt.TagAccession)
        assert result.cv == pt.CV.XL_MOD
