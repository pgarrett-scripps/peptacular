"""
Tests for parsing custom modification tags.
"""

import pytest

import peptacular as pt


class TestTagCustom:
    """Tests for pt.TagCustom.from_string() function"""

    def test_simple_custom_tag(self):
        """Test parsing simple custom tag"""
        result = pt.TagCustom.from_string("C:MyCustomMod")
        assert isinstance(result, pt.TagCustom)
        assert result.name == "MyCustomMod"

    def test_custom_tag_with_underscores(self):
        """Test custom tag with underscores"""
        result = pt.TagCustom.from_string("C:my_custom_mod")
        assert result.name == "my_custom_mod"

    def test_custom_tag_with_numbers(self):
        """Test custom tag with numbers"""
        result = pt.TagCustom.from_string("C:CustomMod123")
        assert result.name == "CustomMod123"

    def test_custom_tag_with_special_chars(self):
        """Test custom tag with special characters"""
        result = pt.TagCustom.from_string("C:Mod-v2.0")
        assert result.name == "Mod-v2.0"

    def test_missing_prefix_raises_error(self):
        """Test that missing C: prefix raises ValueError"""
        with pytest.raises(ValueError, match="must start with 'C:'"):
            pt.TagCustom.from_string("CustomMod")

    def test_empty_name_after_prefix(self):
        """Test custom tag with empty name after prefix"""
        # This should work - name is empty string
        result = pt.TagCustom.from_string("C:")
        assert result.name == ""

    def test_string_representation(self):
        """Test string representation of custom tag"""
        result = pt.TagCustom.from_string("C:MyMod")
        assert str(result) == "C:MyMod"

    def test_via_modification_string(self):
        """Test parsing custom tag via parse_modification_string"""
        result = pt.ModificationTags.from_string("C:MyCustomMod").tags[0]
        assert isinstance(result, pt.TagCustom)
        assert result.name == "MyCustomMod"
