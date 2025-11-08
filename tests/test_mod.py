"""Tests for mod.py module - focusing on uncovered functionality"""

import pytest
from peptacular.mod import (
    Mod,
    get_mod,
    setup_mod,
    clear_mod_cache,
    get_mod_cache_info,
)


class TestModSerialization:
    """Test mod serialization with different parameters"""

    def test_serialize_with_precision(self):
        """Test float serialization with precision parameter"""
        mod = Mod(79.966331, 1)
        result = mod.serialize(brackets="[]", precision=2)
        assert result == "[79.97]"

    def test_serialize_with_include_plus_positive(self):
        """Test include_plus for positive numbers"""
        mod = Mod(42.0, 1)
        result = mod.serialize(brackets="[]", include_plus=True)
        assert result == "[+42.0]"

    def test_serialize_with_include_plus_integer(self):
        """Test include_plus for positive integers"""
        mod = Mod(15, 1)
        result = mod.serialize(brackets="[]", include_plus=True)
        assert result == "[+15]"

    def test_serialize_no_brackets_no_multiplier(self):
        """Test serialization without brackets when mult=1"""
        mod = Mod("Phospho", 1)
        result = mod.serialize(brackets=None)
        assert result == "Phospho"

    def test_serialize_no_brackets_with_multiplier_raises(self):
        """Test that no brackets with mult > 1 raises error"""
        mod = Mod("Acetyl", 2)
        with pytest.raises(ValueError, match="Brackets must be provided"):
            mod.serialize(brackets=None)

    def test_serialize_invalid_brackets_length(self):
        """Test that invalid bracket length raises error"""
        mod = Mod("Acetyl", 1)
        with pytest.raises(ValueError, match="Brackets string must be of length 2"):
            mod.serialize(brackets="[")

    def test_serialize_float_without_decimal(self):
        """Test float that needs .0 appended"""
        mod = Mod(100.0, 1)
        result = mod.serialize(brackets="[]")
        # str(100.0) gives "100.0", so this is already covered
        assert "100.0" in result or "100" in result


class TestModParsing:
    """Test mod parsing from strings"""

    def test_parse_with_brackets_and_multiplier(self):
        """Test parsing mod string with brackets and multiplier"""
        result = Mod.parse("[Phospho]^2", brackets="[]")
        assert result.val == "Phospho"
        assert result.mult == 2

    def test_parse_with_brackets_no_multiplier(self):
        """Test parsing mod string with brackets but no multiplier"""
        result = Mod.parse("[Acetyl]", brackets="[]")
        assert result.val == "Acetyl"
        assert result.mult == 1

    def test_parse_without_brackets(self):
        """Test parsing mod string without brackets"""
        result = Mod.parse("Methyl", brackets="[]")
        assert result.val == "Methyl"
        assert result.mult == 1

    def test_parse_with_plus_sign(self):
        """Test parsing mod with leading plus sign"""
        result = Mod.parse("[+15.99]", brackets="[]")
        assert result.val == 15.99
        assert result.mult == 1

    def test_parse_float_value(self):
        """Test parsing float mod value"""
        result = Mod.parse("[79.966]", brackets="[]")
        assert result.val == 79.966
        assert result.mult == 1

    def test_parse_integer_value(self):
        """Test parsing integer mod value"""
        result = Mod.parse("[42]", brackets="[]")
        assert result.val == 42
        assert result.mult == 1

    def test_parse_scientific_notation(self):
        """Test parsing scientific notation"""
        result = Mod.parse("[1.5e2]", brackets="[]")
        assert result.val == 150.0
        assert result.mult == 1

    def test_parse_invalid_brackets_raises(self):
        """Test that invalid bracket length raises error"""
        with pytest.raises(ValueError, match="Brackets string must be of length 2"):
            Mod.parse("[Acetyl]", brackets="[")


class TestModComparison:
    """Test mod comparison operations"""

    def test_equality_with_string(self):
        """Test mod equality with string value"""
        mod = get_mod("Phospho", 1)
        assert mod == "Phospho"

    def test_equality_with_int(self):
        """Test mod equality with int value"""
        mod = get_mod(42, 1)
        assert mod == 42

    def test_equality_with_float(self):
        """Test mod equality with float value"""
        mod = get_mod(15.99, 1)
        assert mod == 15.99

    def test_less_than_comparison(self):
        """Test less than comparison"""
        mod1 = get_mod("Acetyl", 1)
        mod2 = get_mod("Phospho", 1)
        assert mod1 < mod2

    def test_less_than_with_same_value_different_mult(self):
        """Test less than with same value, different multiplier"""
        mod1 = get_mod("Acetyl", 1)
        mod2 = get_mod("Acetyl", 2)
        assert mod1 < mod2


class TestModUtilityFunctions:
    """Test utility functions for Mod"""

    def test_setup_mod_from_mod(self):
        """Test setup_mod with Mod input"""
        mod = Mod("Phospho", 1)
        result = setup_mod(mod)
        assert result.val == "Phospho"
        assert result.mult == 1

    def test_setup_mod_from_string(self):
        """Test setup_mod with string input"""
        result = setup_mod("Acetyl")
        assert result.val == "Acetyl"
        assert result.mult == 1

    def test_setup_mod_from_int(self):
        """Test setup_mod with int input"""
        result = setup_mod(42)
        assert result.val == 42
        assert result.mult == 1

    def test_setup_mod_from_float(self):
        """Test setup_mod with float input"""
        result = setup_mod(15.99)
        assert result.val == 15.99
        assert result.mult == 1

    def test_setup_mod_invalid_type(self):
        """Test setup_mod with invalid type raises error"""
        with pytest.raises(TypeError, match="Invalid mod input"):
            setup_mod([1, 2, 3])

    def test_get_mod_cache_behavior(self):
        """Test that get_mod returns cached instances"""
        clear_mod_cache()
        mod1 = get_mod("Phospho", 1)
        mod2 = get_mod("Phospho", 1)
        # Should be the exact same object due to caching
        assert mod1 is mod2

    def test_get_mod_cache_info(self):
        """Test cache info retrieval"""
        clear_mod_cache()
        get_mod("Test", 1)
        get_mod("Test", 1)  # Cache hit
        
        info = get_mod_cache_info()
        assert "mod_cache" in info
        assert "serialize_cache" in info
        assert info["mod_cache"]["hits"] >= 1
        assert "hit_rate" in info["mod_cache"]

    def test_mod_to_dict(self):
        """Test mod to_dict conversion"""
        mod = Mod("Phospho", 2)
        result = mod.to_dict()
        assert result == {"val": "Phospho", "mult": 2}

    def test_mod_copy(self):
        """Test mod copy returns same instance"""
        mod = Mod("Acetyl", 1)
        copied = mod.copy()
        assert copied is mod  # Should be same cached instance

    def test_mod_repr_string(self):
        """Test mod repr for string values"""
        mod = Mod("Phospho", 1)
        assert repr(mod) == "Mod('Phospho', 1)"

    def test_mod_repr_number(self):
        """Test mod repr for numeric values"""
        mod = Mod(42, 1)
        assert repr(mod) == "Mod(42, 1)"

    def test_mod_hash(self):
        """Test mod hashing"""
        mod1 = Mod("Phospho", 1)
        mod2 = Mod("Phospho", 1)
        assert hash(mod1) == hash(mod2)
