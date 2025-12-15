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
            pt.parse_modification_tag("")

    def test_whitespace_only(self):
        """Test parsing whitespace-only string raises error"""
        with pytest.raises(ValueError, match="Empty modification string"):
            pt.parse_modification_tag("   ")

    def test_unknown_cv_prefix_treated_as_name(self):
        """Test that unknown CV prefix is treated as a name"""
        result = pt.parse_modification_tag("UNKNOWN:123")
        assert isinstance(result, pt.TagName)
        assert result.name == "UNKNOWN:123"

    def test_invalid_formula_element(self):
        """Test parsing formula with invalid element"""
        with pytest.raises(ValueError, match="Unknown element symbol"):
            pt.parse_modification_tag("Formula:Zz")

    def test_invalid_glycan_composition(self):
        """Test parsing invalid glycan composition"""
        with pytest.raises(ValueError, match="Could not parse glycan composition"):
            pt.parse_modification_tag("Glycan:InvalidMono")
