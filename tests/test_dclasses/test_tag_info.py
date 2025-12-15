"""
Tests for parsing INFO tag modifications.
"""

import peptacular as pt


class TestTagInfo:
    """Tests for parsing INFO tag modifications"""

    def test_simple_info(self):
        """Test parsing simple INFO tag"""
        result = pt.parse_modification_tag("INFO:custom annotation")
        assert isinstance(result, pt.TagInfo)
        assert result.info == "custom annotation"

    def test_info_with_special_chars(self):
        """Test parsing INFO tag with special characters"""
        result = pt.parse_modification_tag("INFO:annotation-with_special.chars")
        assert isinstance(result, pt.TagInfo)
        assert result.info == "annotation-with_special.chars"

    def test_case_insensitivity(self):
        """Test that INFO tag parsing is case-insensitive"""
        result1 = pt.parse_modification_tag("info:lowercase")
        result2 = pt.parse_modification_tag("INFO:UPPERCASE")
        assert isinstance(result1, pt.TagInfo)
        assert isinstance(result2, pt.TagInfo)
        assert result1.info == "lowercase"
        assert result2.info == "UPPERCASE"
