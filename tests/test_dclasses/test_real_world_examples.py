"""
Tests using real-world ProForma modification examples.
"""

import peptacular as pt


class TestRealWorldExamples:
    """Tests using real-world ProForma modification examples"""

    def test_oxidation_variants(self):
        """Test different ways to represent oxidation"""
        # By name
        result1 = pt.ModificationTags.from_string("Oxidation").tags[0]
        assert isinstance(result1, pt.TagName)

        # By UNIMOD
        result2 = pt.ModificationTags.from_string("UNIMOD:35").tags[0]
        assert isinstance(result2, pt.TagAccession)

        # By mass
        result3 = pt.ModificationTags.from_string("+15.995").tags[0]
        assert isinstance(result3, pt.TagMass)

        # By formula
        result4 = pt.ModificationTags.from_string("Formula:O").tags[0]
        assert isinstance(result4, pt.ChargedFormula)

    def test_phosphorylation_variants(self):
        """Test different ways to represent phosphorylation"""
        # By name
        result1 = pt.ModificationTags.from_string("Phospho").tags[0]
        assert isinstance(result1, pt.TagName)

        # By PSI-MOD
        result2 = pt.ModificationTags.from_string("MOD:00046").tags[0]
        assert isinstance(result2, pt.TagAccession)

        # By mass
        result3 = pt.ModificationTags.from_string("+79.966").tags[0]
        assert isinstance(result3, pt.TagMass)

    def test_crosslinker_modification(self):
        """Test cross-linking modification"""
        result = pt.ModificationTags.from_string("XLMOD:02001").tags[0]
        assert isinstance(result, pt.TagAccession)
        assert result.cv == pt.CV.XL_MOD
