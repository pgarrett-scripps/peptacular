"""
Tests for parsing INFO tag modifications.
"""

import peptacular as pt


class TestTagInfo:
    """Tests for parsing INFO tag modifications"""

    def test_simple_info(self):
        """Test parsing simple INFO tag"""
        result = pt.ModificationTags.from_string("INFO:custom annotation").tags[0]
        assert isinstance(result, pt.TagInfo)
        assert result.info == "custom annotation"

    def test_info_with_special_chars(self):
        """Test parsing INFO tag with special characters"""
        result = pt.ModificationTags.from_string("INFO:annotation-with_special.chars").tags[0]
        assert isinstance(result, pt.TagInfo)
        assert result.info == "annotation-with_special.chars"

    def test_case_insensitivity(self):
        """Test that INFO tag parsing is case-insensitive"""
        result1 = pt.ModificationTags.from_string("info:lowercase").tags[0]
        result2 = pt.ModificationTags.from_string("INFO:UPPERCASE").tags[0]
        assert isinstance(result1, pt.TagInfo)
        assert isinstance(result2, pt.TagInfo)
        assert result1.info == "lowercase"
        assert result2.info == "UPPERCASE"
