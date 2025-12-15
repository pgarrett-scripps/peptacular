"""
Tests for the ElementLookup class.
"""

import pytest
import warnings
import peptacular as pt

from peptacular.elements.lookup import _handle_key_input  # type: ignore


class TestHandleKeyInput:
    """Tests for the _handle_key_input helper function"""

    def test_tuple_key_valid(self):
        """Test valid tuple key parsing"""
        result = _handle_key_input(("C", 12))
        assert result == ("C", 12)

    def test_tuple_key_with_none(self):
        """Test tuple key with None mass number"""
        result = _handle_key_input(("C", None))
        assert result == ("C", None)

    def test_string_symbol_only(self):
        """Test string with symbol only"""
        result = _handle_key_input("C")
        assert result == ("C", None)

    def test_string_with_mass_prefix(self):
        """Test string with mass number prefix like '13C'"""
        result = _handle_key_input("13C")
        assert result == ("C", 13)

    def test_string_deuterium(self):
        """Test deuterium symbol 'D'"""
        result = _handle_key_input("D")
        assert result == ("D", None)

    def test_string_2H(self):
        """Test '2H' notation for deuterium"""
        result = _handle_key_input("2H")
        assert result == ("H", 2)

    def test_invalid_tuple_length(self):
        """Test tuple with wrong number of elements"""
        with pytest.raises(ValueError, match="must have exactly 2 elements"):
            _handle_key_input(("C", 12, "extra"))  # type: ignore

    def test_invalid_symbol_type(self):
        """Test tuple with non-string symbol"""
        with pytest.raises(TypeError, match="Symbol must be str"):
            _handle_key_input((123, 12))  # type: ignore

    def test_invalid_mass_type(self):
        """Test tuple with invalid mass number type"""
        with pytest.raises(TypeError, match="Mass number must be int or None"):
            _handle_key_input(("C", "twelve"))  # type: ignore

    def test_empty_string(self):
        """Test empty string raises error"""
        with pytest.raises(ValueError, match="cannot be empty string"):
            _handle_key_input("")

    def test_lowercase_symbol(self):
        """Test symbol not starting with uppercase"""
        with pytest.raises(ValueError, match="must start with uppercase"):
            _handle_key_input("carbon")

    def test_invalid_isotope_notation(self):
        """Test isotope notation with only digits"""
        with pytest.raises(ValueError, match="no element symbol found"):
            _handle_key_input("123")

    def test_invalid_key_type(self):
        """Test completely invalid key type"""
        with pytest.raises(TypeError, match="Key must be tuple"):
            _handle_key_input(123)  # type: ignore


class TestElementLookupBasics:
    """Tests for basic ElementLookup functionality"""

    def test_lookup_by_tuple_exact_isotope(self):
        """Test looking up exact isotope by tuple"""
        carbon_12 = pt.ELEMENT_LOOKUP[("C", 12)]
        assert carbon_12.symbol == "C"
        assert carbon_12.mass_number == 12
        assert carbon_12.mass == pytest.approx(12.0, abs=0.1)  # type: ignore

    def test_lookup_by_tuple_monoisotopic(self):
        """Test looking up monoisotopic isotope by tuple with None"""
        carbon_mono = pt.ELEMENT_LOOKUP[("C", None)]
        assert carbon_mono.symbol == "C"
        assert carbon_mono.mass_number == 12
        assert carbon_mono.mass == pytest.approx(12.0, abs=0.1)  # type: ignore

    def test_lookup_by_string_symbol(self):
        """Test looking up by symbol string returns monoisotopic"""
        carbon = pt.ELEMENT_LOOKUP["C"]
        assert carbon.symbol == "C"
        assert carbon.mass_number == 12

    def test_lookup_by_string_with_mass(self):
        """Test looking up by string with mass prefix like '13C'"""
        carbon_13 = pt.ELEMENT_LOOKUP["13C"]
        assert carbon_13.symbol == "C"
        assert carbon_13.mass_number == 13
        assert carbon_13.mass == pytest.approx(13.003, abs=0.01)  # type: ignore

    def test_lookup_deuterium_by_D(self):
        """Test looking up deuterium by 'D' symbol"""
        deuterium = pt.ELEMENT_LOOKUP["D"]
        assert deuterium.symbol == "D"
        assert deuterium.mass_number == 2

    def test_lookup_deuterium_by_2H(self):
        """Test looking up deuterium by '2H' notation"""
        deuterium = pt.ELEMENT_LOOKUP["2H"]
        # Note: In the database, H-2 is stored as 'D' (deuterium)
        assert deuterium.symbol == "D"
        assert deuterium.mass_number == 2

    def test_lookup_tritium(self):
        """Test looking up tritium"""
        tritium = pt.ELEMENT_LOOKUP["T"]
        assert tritium.symbol == "T"
        assert tritium.mass_number == 3

    def test_lookup_nonexistent_element(self):
        """Test looking up non-existent element raises KeyError"""
        with pytest.raises(KeyError):
            _ = pt.ELEMENT_LOOKUP["Xx"]

    def test_lookup_nonexistent_monoisotopic(self):
        """Test looking up non-existent monoisotopic entry"""
        with pytest.raises(KeyError, match="Monoisotopic isotope"):
            _ = pt.ELEMENT_LOOKUP[("Zz", None)]


class TestElementLookupAutoGeneration:
    """Tests for automatic isotope generation"""

    def test_auto_generate_carbon_16(self):
        """Test auto-generating Carbon-16 (not in standard database)"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            carbon_16 = pt.ELEMENT_LOOKUP["16C"]

            # Check warning was issued
            assert len(w) == 1
            assert "not found in database" in str(w[0].message)
            assert "Generated automatically" in str(w[0].message)

            # Check the isotope properties
            assert carbon_16.symbol == "C"
            assert carbon_16.mass_number == 16
            assert carbon_16.abundance == 0.0

            # Mass should be approximately 16.0 (12.0 + 4 * neutron_mass)
            expected_mass = 12.0 + 4 * pt.ElementLookup.NEUTRON_MASS
            assert carbon_16.mass == pytest.approx(expected_mass, abs=0.01)  # type: ignore

    def test_auto_generate_cached(self):
        """Test that auto-generated isotope is cached"""
        # First access generates it
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            carbon_15_first = pt.ELEMENT_LOOKUP["15C"]
            assert len(w) == 1  # Warning issued

        # Second access uses cached version
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            carbon_15_second = pt.ELEMENT_LOOKUP["15C"]
            assert len(w) == 0  # No warning

        # Should be the same object
        assert carbon_15_first is carbon_15_second

    def test_auto_generate_lighter_isotope(self):
        """Test generating lighter isotope (subtracting neutrons)"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            carbon_11 = pt.ELEMENT_LOOKUP["11C"]

            assert len(w) == 1
            assert "subtracting" in str(w[0].message)

            assert carbon_11.mass_number == 11
            expected_mass = 12.0 - 1 * pt.ElementLookup.NEUTRON_MASS
            assert carbon_11.mass == pytest.approx(expected_mass, abs=0.01)  # type: ignore

    def test_auto_generate_nonexistent_element_fails(self):
        """Test that auto-generation fails for completely unknown element"""
        with pytest.raises(KeyError, match="No monoisotopic entry found"):
            _ = pt.ELEMENT_LOOKUP[("Zz", 123)]


class TestElementLookupContains:
    """Tests for the __contains__ method"""

    def test_contains_existing_symbol(self):
        """Test that existing symbol is found"""
        assert "C" in pt.ELEMENT_LOOKUP

    def test_contains_existing_isotope_string(self):
        """Test that existing isotope string is found"""
        assert "13C" in pt.ELEMENT_LOOKUP

    def test_contains_existing_tuple(self):
        """Test that existing tuple key is found"""
        assert ("C", 12) in pt.ELEMENT_LOOKUP

    def test_not_contains_nonexistent(self):
        """Test that non-existent element is not found"""
        assert "Xx" not in pt.ELEMENT_LOOKUP

    def test_not_contains_auto_generated(self):
        """Test that auto-generated isotopes are not in __contains__"""
        # This should NOT trigger auto-generation
        assert ("C", 999) not in pt.ELEMENT_LOOKUP

    def test_contains_invalid_key(self):
        """Test that invalid key returns False"""
        assert 123 not in pt.ELEMENT_LOOKUP  # type: ignore
        assert "" not in pt.ELEMENT_LOOKUP


class TestElementLookupMethods:
    """Tests for ElementLookup convenience methods"""

    def test_get_monoisotopic(self):
        """Test get_monoisotopic method"""
        carbon = pt.ELEMENT_LOOKUP.get_monoisotopic("C")
        assert carbon.symbol == "C"
        assert carbon.mass_number == 12

    def test_get_isotope(self):
        """Test get_isotope method"""
        carbon_13 = pt.ELEMENT_LOOKUP.get_isotope("C", 13)
        assert carbon_13.symbol == "C"
        assert carbon_13.mass_number == 13

    def test_get_isotope_auto_generate_disabled(self):
        """Test get_isotope with auto_generate=False"""
        with pytest.raises(KeyError, match="auto_generate=False"):
            pt.ELEMENT_LOOKUP.get_isotope("C", 999, auto_generate=False)

    def test_get_all_isotopes(self):
        """Test get_all_isotopes method"""
        carbon_isotopes = pt.ELEMENT_LOOKUP.get_all_isotopes("C")
        assert len(carbon_isotopes) > 0
        assert all(iso.symbol == "C" for iso in carbon_isotopes)
        assert all(iso.mass_number is not None for iso in carbon_isotopes)

        # Should be sorted by mass number
        mass_numbers = [iso.mass_number for iso in carbon_isotopes]
        assert mass_numbers == sorted(mass_numbers)

    def test_get_all_isotopes_exclude_generated(self):
        """Test that get_all_isotopes excludes generated isotopes by default"""
        # First generate a fake isotope
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _ = pt.ELEMENT_LOOKUP["99C"]  # This will be auto-generated

        # Get all isotopes - should not include the generated one
        carbon_isotopes = pt.ELEMENT_LOOKUP.get_all_isotopes(
            "C", include_generated=False
        )
        mass_numbers = [iso.mass_number for iso in carbon_isotopes]
        assert 99 not in mass_numbers

    def test_get_all_isotopes_include_generated(self):
        """Test that get_all_isotopes can include generated isotopes"""
        # Generate a fake isotope
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _ = pt.ELEMENT_LOOKUP["98C"]

        # Get all isotopes including generated
        carbon_isotopes = pt.ELEMENT_LOOKUP.get_all_isotopes(
            "C", include_generated=True
        )
        mass_numbers = [iso.mass_number for iso in carbon_isotopes]
        assert 98 in mass_numbers

    def test_get_all_isotopes_nonexistent(self):
        """Test get_all_isotopes for non-existent element"""
        with pytest.raises(KeyError, match="No isotopes found"):
            pt.ELEMENT_LOOKUP.get_all_isotopes("Xx")

    def test_get_elements(self):
        """Test get_elements method"""
        elements = pt.ELEMENT_LOOKUP.get_elements()
        assert isinstance(elements, list)
        assert "C" in elements
        assert "H" in elements
        assert "N" in elements
        assert "O" in elements

        # Should be sorted
        assert elements == sorted(elements)

    def test_mass_monoisotopic(self):
        """Test mass method with monoisotopic=True"""
        carbon_mass = pt.ELEMENT_LOOKUP.mass("C", monoisotopic=True)
        assert carbon_mass == pytest.approx(12.0, abs=0.1)  # type: ignore

    def test_mass_average(self):
        """Test mass method with monoisotopic=False"""
        carbon_avg = pt.ELEMENT_LOOKUP.mass("C", monoisotopic=False)
        # Average mass should be slightly higher than 12 due to C-13
        assert carbon_avg > 12.0
        assert carbon_avg == pytest.approx(12.011, abs=0.001)  # type: ignore

    def test_mass_specific_isotope_ignores_monoisotopic_flag(self):
        """Test that specific isotope mass ignores monoisotopic parameter"""
        # When requesting specific isotope, should always return exact mass
        c13_mass_mono = pt.ELEMENT_LOOKUP.mass("13C", monoisotopic=True)
        c13_mass_avg = pt.ELEMENT_LOOKUP.mass("13C", monoisotopic=False)

        # Both should return the same exact isotope mass
        assert c13_mass_mono == c13_mass_avg
        assert c13_mass_mono == pytest.approx(13.003, abs=0.01)  # type: ignore

    def test_mass_tuple_specific_isotope(self):
        """Test mass method with tuple for specific isotope"""
        c13_mass = pt.ELEMENT_LOOKUP.mass(("C", 13))
        assert c13_mass == pytest.approx(13.003, abs=0.01)  # type: ignore


class TestElementLookupSpecialCases:
    """Tests for special cases and edge conditions"""

    def test_len(self):
        """Test __len__ method"""
        length = len(pt.ELEMENT_LOOKUP)
        assert length > 0
        assert isinstance(length, int)

    def test_repr(self):
        """Test __repr__ method"""
        repr_str = repr(pt.ELEMENT_LOOKUP)
        assert "ElementLookup" in repr_str
        assert "entries" in repr_str
        assert "elements" in repr_str

    def test_multiple_formats_same_isotope(self):
        """Test that different formats return same isotope"""
        c12_tuple = pt.ELEMENT_LOOKUP[("C", 12)]
        c12_string = pt.ELEMENT_LOOKUP["12C"]

        # Should be the same object or at least same values
        assert c12_tuple.symbol == c12_string.symbol
        assert c12_tuple.mass_number == c12_string.mass_number
        assert c12_tuple.mass == c12_string.mass

    def test_neutron_mass_constant(self):
        """Test that neutron mass constant is reasonable"""
        assert pt.ElementLookup.NEUTRON_MASS == pytest.approx(1.0087, abs=0.001)  # type: ignore

    def test_element_info_attributes(self):
        """Test that returned ElementInfo has all expected attributes"""
        carbon = pt.ELEMENT_LOOKUP["C"]
        assert hasattr(carbon, "number")
        assert hasattr(carbon, "symbol")
        assert hasattr(carbon, "mass_number")
        assert hasattr(carbon, "mass")
        assert hasattr(carbon, "abundance")
        assert hasattr(carbon, "average_mass")

        assert isinstance(carbon.number, int)
        assert isinstance(carbon.symbol, str)
        assert isinstance(carbon.mass_number, int)
        assert isinstance(carbon.mass, float)
        assert isinstance(carbon.abundance, float)
        assert isinstance(carbon.average_mass, float)
