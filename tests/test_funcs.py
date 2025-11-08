"""Tests for funcs.py module"""

from peptacular.funcs import (
    ModLabler,
    get_label,
    get_number,
)


class TestModLabler:
    """Test ModLabler class"""

    def test_modlabler_num_increment(self):
        """Test numeric label increment"""
        labeler = ModLabler()
        assert labeler.num == "1"
        labeler.increment_num()
        assert labeler.num == "2"

    def test_modlabler_letter_increment(self):
        """Test letter label increment"""
        labeler = ModLabler()
        assert labeler.letter == "a"
        labeler.increment_letter()
        assert labeler.letter == "b"

    def test_modlabler_letter_multi_char(self):
        """Test letter label for indices beyond 'z'"""
        labeler = ModLabler()
        labeler.letter_index = 26
        assert labeler.letter == "aa"

    def test_modlabler_combined_label(self):
        """Test combined letter+number label"""
        labeler = ModLabler()
        assert labeler.label == "a1"
        labeler.increment_num()
        labeler.increment_letter()
        assert labeler.label == "b2"


class TestGetLabel:
    """Test get_label function"""

    def test_get_label_with_loss(self):
        """Test get_label with neutral loss"""
        result = get_label("b", 1, "2", -18.0, 0)
        assert result == "+b2(-18.0)"

    def test_get_label_with_isotope(self):
        """Test get_label with isotope"""
        result = get_label("y", 1, "3", 0.0, 2)
        assert result == "+y3**"

    def test_get_label_with_precision(self):
        """Test get_label with precision rounding"""
        result = get_label("b", 1, "2", -17.9876, 0, precision=2)
        assert result == "+b2(-17.99)"

    def test_get_label_double_charge(self):
        """Test get_label with double charge"""
        result = get_label("y", 2, "5", 0.0, 0)
        assert result == "++y5"


class TestGetNumber:
    """Test get_number function"""

    def test_get_number_forward_ion(self):
        """Test get_number for forward ion (b type)"""
        result = get_number("b", 10, 0, 5)
        assert result == "5"  # Returns string

    def test_get_number_backward_ion(self):
        """Test get_number for backward ion (y type)"""
        result = get_number("y", 10, 3, 10)
        assert result == "7"  # 10 - 3, returns string

    def test_get_number_internal_ion(self):
        """Test get_number for internal ion"""
        result = get_number("by", 10, 2, 7)
        assert result == "2-7"
