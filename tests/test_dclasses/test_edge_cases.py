"""
Tests for edge cases and error conditions.
"""

import pytest

import peptacular as pt


class TestEdgeCases:
    """Tests for edge cases and error conditions"""

    def test_empty_string(self):
        """Test parsing empty string raises error"""
        with pytest.raises(ValueError, match="Empty modification string"):
            pt.ModificationTags.from_string("").tags[0]

    def test_whitespace_only(self):
        """Test parsing whitespace-only string raises error"""
        with pytest.raises(ValueError, match="Empty modification string"):
            pt.ModificationTags.from_string("   ").tags[0]

    def test_unknown_cv_prefix_treated_as_name(self):
        """Test that unknown CV prefix is treated as a name"""
        result = pt.ModificationTags.from_string("UNKNOWN:123").tags[0]
        assert isinstance(result, pt.TagName)
        assert result.name == "UNKNOWN:123"

    def test_invalid_formula_element(self):
        """Test parsing formula with invalid element"""
        with pytest.raises(ValueError, match="Unknown element symbol"):
            pt.ModificationTags.from_string("Formula:Zz").tags[0]

    def test_invalid_glycan_composition(self):
        """Test parsing invalid glycan composition"""
        with pytest.raises(ValueError, match="Could not parse glycan composition"):
            pt.ModificationTags.from_string("Glycan:InvalidMono").tags[0]
