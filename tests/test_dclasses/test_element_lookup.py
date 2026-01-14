"""
Tests for the ElementLookup class.
"""

import pytest

import peptacular as pt
from peptacular.elements.lookup import _handle_key_input


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
        assert result == ("H", 2)

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
        assert carbon_12.mass == pytest.approx(12.0, abs=0.1)

    def test_lookup_by_tuple_monoisotopic(self):
        """Test looking up monoisotopic isotope by tuple with None"""
        carbon_mono = pt.ELEMENT_LOOKUP[("C", None)]
        assert carbon_mono.symbol == "C"
        assert carbon_mono.mass_number is None
        assert carbon_mono.mass == pytest.approx(12.0, abs=0.1)

    def test_lookup_by_string_symbol(self):
        """Test looking up by symbol string returns monoisotopic"""
        carbon = pt.ELEMENT_LOOKUP["C"]
        assert carbon.symbol == "C"
        assert carbon.mass_number is None

    def test_lookup_by_string_with_mass(self):
        """Test looking up by string with mass prefix like '13C'"""
        carbon_13 = pt.ELEMENT_LOOKUP["13C"]
        assert carbon_13.symbol == "C"
        assert carbon_13.mass_number == 13
        assert carbon_13.mass == pytest.approx(13.003, abs=0.01)

    def test_lookup_deuterium_by_D(self):
        """Test looking up deuterium by 'D' symbol"""
        deuterium = pt.ELEMENT_LOOKUP["D"]
        assert deuterium.symbol == "H"
        assert deuterium.mass_number == 2

    def test_lookup_deuterium_by_2H(self):
        """Test looking up deuterium by '2H' notation"""
        deuterium = pt.ELEMENT_LOOKUP["2H"]
        # Note: In the database, H-2 is stored as 'D' (deuterium)
        assert deuterium.symbol == "H"
        assert deuterium.mass_number == 2

    def test_lookup_tritium(self):
        """Test looking up tritium"""
        tritium = pt.ELEMENT_LOOKUP["T"]
        assert tritium.symbol == "H"
        assert tritium.mass_number == 3

    def test_lookup_nonexistent_element(self):
        """Test looking up non-existent element raises KeyError"""
        with pytest.raises(KeyError):
            _ = pt.ELEMENT_LOOKUP["Xx"]

    def test_lookup_nonexistent_monoisotopic(self):
        """Test looking up non-existent monoisotopic entry"""
        with pytest.raises(KeyError):
            _ = pt.ELEMENT_LOOKUP[("Zz", None)]


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

    def test_not_contains_nonexistent_isotope(self):
        """Test that non-existent isotopes return False"""
        # This should NOT raise an error, just return False
        assert ("C", 999) not in pt.ELEMENT_LOOKUP
        assert "999C" not in pt.ELEMENT_LOOKUP

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

    def test_get_all_isotopes(self):
        """Test get_all_isotopes method"""
        carbon_isotopes = pt.ELEMENT_LOOKUP.get_all_isotopes("C")
        assert len(carbon_isotopes) > 0
        assert all(iso.symbol == "C" for iso in carbon_isotopes)
        assert all(iso.mass_number is not None for iso in carbon_isotopes)

        # Should be sorted by mass number
        mass_numbers = [iso.mass_number for iso in carbon_isotopes]
        assert mass_numbers == sorted([m for m in mass_numbers if m is not None])

    def test_get_all_isotopes_nonexistent(self):
        """Test get_all_isotopes for non-existent element"""
        with pytest.raises(KeyError):
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
        assert carbon_mass == pytest.approx(12.0, abs=0.1)

    def test_mass_average(self):
        """Test mass method with monoisotopic=False"""
        carbon_avg = pt.ELEMENT_LOOKUP.mass("C", monoisotopic=False)
        # Average mass should be slightly higher than 12 due to C-13
        assert carbon_avg > 12.0
        assert carbon_avg == pytest.approx(12.011, abs=0.001)

    def test_mass_specific_isotope_ignores_monoisotopic_flag(self):
        """Test that specific isotope mass ignores monoisotopic parameter"""
        # When requesting specific isotope, should always return exact mass
        c13_mass_mono = pt.ELEMENT_LOOKUP.mass("13C", monoisotopic=True)
        c13_mass_avg = pt.ELEMENT_LOOKUP.mass("13C", monoisotopic=False)

        # Both should return the same exact isotope mass
        assert c13_mass_mono == c13_mass_avg
        assert c13_mass_mono == pytest.approx(13.003, abs=0.01)

    def test_mass_tuple_specific_isotope(self):
        """Test mass method with tuple for specific isotope"""
        c13_mass = pt.ELEMENT_LOOKUP.mass(("C", 13))
        assert c13_mass == pytest.approx(13.003, abs=0.01)


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

    def test_isotope_abundance_attribute(self):
        """Test that isotopes have proper abundance values"""
        # Natural isotopes should have non-zero abundance
        c12 = pt.ELEMENT_LOOKUP["12C"]
        assert c12.abundance is not None
        assert c12.abundance > 0.98  # C-12 is very abundant

        c13 = pt.ELEMENT_LOOKUP["13C"]
        assert c13.abundance is not None
        assert 0.01 < c13.abundance < 0.02  # C-13 is about 1.1%

    def test_element_info_attributes(self):
        """Test that returned ElementInfo has all expected attributes"""
        carbon = pt.ELEMENT_LOOKUP["12C"]
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


class TestElementLookupNeutronOffsetsAndAbundances:
    """Tests for get_neutron_offsets_and_abundances method"""

    def test_carbon_neutron_offsets(self):
        """Test neutron offsets for carbon isotopes"""
        offsets = pt.ELEMENT_LOOKUP.get_neutron_offsets_and_abundances("C")

        # Should have at least C-12 and C-13
        assert len(offsets) >= 2

        # Check format is list of (offset, abundance) tuples
        for offset, abundance in offsets:
            assert isinstance(offset, int)
            assert isinstance(abundance, (float, type(None)))

        # C-12 is monoisotopic, so should have offset 0
        offsets_dict = {off: abund for off, abund in offsets}
        assert 0 in offsets_dict  # C-12 offset
        assert 1 in offsets_dict  # C-13 offset (one extra neutron)

    def test_hydrogen_neutron_offsets(self):
        """Test neutron offsets for hydrogen isotopes"""
        offsets = pt.ELEMENT_LOOKUP.get_neutron_offsets_and_abundances("H")

        # H-1 (protium) is monoisotopic
        offsets_dict = {off: abund for off, abund in offsets}
        assert 0 in offsets_dict  # H-1
        assert 1 in offsets_dict  # H-2 (deuterium, +1 neutron)

    def test_nitrogen_neutron_offsets(self):
        """Test neutron offsets for nitrogen isotopes"""
        offsets = pt.ELEMENT_LOOKUP.get_neutron_offsets_and_abundances("N")

        # N-14 is monoisotopic
        offsets_dict = {off: abund for off, abund in offsets}
        assert 0 in offsets_dict  # N-14
        assert 1 in offsets_dict  # N-15 (+1 neutron)

    def test_oxygen_neutron_offsets(self):
        """Test neutron offsets for oxygen isotopes"""
        offsets = pt.ELEMENT_LOOKUP.get_neutron_offsets_and_abundances("O")

        # O-16 is monoisotopic
        offsets_dict = {off: abund for off, abund in offsets}
        assert 0 in offsets_dict  # O-16
        assert 1 in offsets_dict  # O-17
        assert 2 in offsets_dict  # O-18 (+2 neutrons)


class TestElementLookupMassesAndAbundances:
    """Tests for get_masses_and_abundances method"""

    def test_carbon_masses_and_abundances(self):
        """Test masses and abundances for carbon"""
        masses = pt.ELEMENT_LOOKUP.get_masses_and_abundances("C")

        # Should have at least C-12 and C-13
        assert len(masses) >= 2

        # Check format is list of (mass, abundance) tuples
        for mass, abundance in masses:
            assert isinstance(mass, float)
            assert isinstance(abundance, (float, type(None)))
            assert mass > 0  # Mass should be positive

        # Masses should be sorted
        mass_values = [m for m, _ in masses]
        assert mass_values == sorted(mass_values)

        # C-12 should be around 12.0, C-13 around 13.003
        assert any(11.9 < m < 12.1 for m, _ in masses)  # C-12
        assert any(13.0 < m < 13.1 for m, _ in masses)  # C-13

    def test_hydrogen_masses_and_abundances(self):
        """Test masses and abundances for hydrogen"""
        masses = pt.ELEMENT_LOOKUP.get_masses_and_abundances("H")

        # Should have H-1, and H-2 (deuterium)
        assert len(masses) >= 2

        # H-1 should be around 1.008
        mass_values = [m for m, _ in masses]
        assert any(1.0 < m < 1.01 for m in mass_values)  # H-1
        assert any(2.0 < m < 2.02 for m in mass_values)  # H-2

    def test_abundances_sum_close_to_one(self):
        """Test that abundances for natural isotopes sum close to 1"""
        # For elements with only natural isotopes, abundances should sum to ~1
        for element in ["C", "H", "N", "O"]:
            masses = pt.ELEMENT_LOOKUP.get_masses_and_abundances(element)

            # Filter out radioactive isotopes (abundance=0)
            natural = [(m, a) for m, a in masses if a is not None and a > 0]

            if natural:
                total_abundance = sum(a for _, a in natural)
                # Should be close to 1.0 (allowing some rounding errors)
                assert 0.99 < total_abundance < 1.01, (
                    f"{element} abundances sum to {total_abundance}"
                )

    def test_nonexistent_element_raises(self):
        """Test that nonexistent element raises KeyError"""
        with pytest.raises(KeyError):
            pt.ELEMENT_LOOKUP.get_masses_and_abundances("Xx")


class TestParseComposition:
    """Tests for parse_composition helper function"""

    def test_simple_composition(self):
        """Test parsing simple composition"""
        from peptacular.elements import parse_composition

        comp_dict = {"C": 2, "H": 6}
        parsed = parse_composition(comp_dict)

        # Should have 2 entries
        assert len(parsed) == 2

        # Check that keys are ElementInfo objects
        for key, count in parsed.items():
            assert isinstance(key, pt.ElementInfo)
            assert isinstance(count, int)
            assert count > 0

    def test_isotope_composition(self):
        """Test parsing composition with specific isotopes"""
        from peptacular.elements import parse_composition

        comp_dict = {"13C": 2, "D": 4}
        parsed = parse_composition(comp_dict)

        # Should have 2 entries
        assert len(parsed) == 2

        # Check isotopes are correct
        symbols_and_masses = {
            (elem.symbol, elem.mass_number): count for elem, count in parsed.items()
        }

        assert ("C", 13) in symbols_and_masses
        assert symbols_and_masses[("C", 13)] == 2

        # 2H is stored as H with mass number 2
        assert ("H", 2) in symbols_and_masses
        assert symbols_and_masses[("H", 2)] == 4

    def test_mixed_composition(self):
        """Test parsing composition with mix of regular and isotope-specific"""
        from peptacular.elements import parse_composition

        comp_dict = {"C": 2, "13C": 1, "H": 3, "O": 1}
        parsed = parse_composition(comp_dict)

        # Should have 4 entries (C, 13C, H, O)
        assert len(parsed) == 4

    def test_invalid_element_raises(self):
        """Test that invalid element key raises error"""
        from peptacular.elements import parse_composition

        with pytest.raises(KeyError):
            parse_composition({"Xx": 1})

    def test_zero_count_allowed(self):
        """Test that zero counts are allowed"""
        from peptacular.elements import parse_composition

        comp_dict = {"C": 0, "H": 2}
        parsed = parse_composition(comp_dict)

        # Should have both entries even though C has count 0
        assert len(parsed) == 2

        # Find the C entry
        c_count = next(count for elem, count in parsed.items() if elem.symbol == "C")
        assert c_count == 0


class TestElementLookupEdgeCases:
    """Tests for edge cases and error handling"""

    def test_get_isotope_with_none_raises(self):
        """Test that get_isotope with None mass number raises ValueError"""
        with pytest.raises(ValueError, match="cannot be None"):
            pt.ELEMENT_LOOKUP.get_isotope("C", None)  # type: ignore

    def test_lookup_empty_string_raises(self):
        """Test that empty string raises appropriate error"""
        with pytest.raises(ValueError, match="cannot be empty"):
            _ = pt.ELEMENT_LOOKUP[""]

    def test_lookup_invalid_tuple_length(self):
        """Test that tuple with wrong length raises error"""
        with pytest.raises(ValueError, match="exactly 2 elements"):
            _ = pt.ELEMENT_LOOKUP[("C", 12, "extra")]  # type: ignore

    def test_monoisotopic_flag(self):
        """Test that is_monoisotopic flag is set correctly"""
        # C-12 should be monoisotopic
        c12 = pt.ELEMENT_LOOKUP["12C"]
        assert c12.is_monoisotopic is True

        # C-13 should not be monoisotopic
        c13 = pt.ELEMENT_LOOKUP["13C"]
        assert c13.is_monoisotopic is False

        # When requesting by symbol only
        c = pt.ELEMENT_LOOKUP["C"]
        assert c.is_monoisotopic is None

    def test_radioactive_isotopes(self):
        """Test that radioactive isotopes have abundance 0"""
        # C-14 is radioactive
        isotopes = pt.ELEMENT_LOOKUP.get_all_isotopes("C")
        c14 = next((iso for iso in isotopes if iso.mass_number == 14), None)

        if c14 is not None:
            assert c14.is_radioactive is True
            assert c14.abundance == 0.0

    def test_element_properties(self):
        """Test ElementInfo properties"""
        c12 = pt.ELEMENT_LOOKUP["12C"]

        # Check derived properties
        assert c12.proton_count == 6
        assert c12.neutron_count == 6
        assert c12.number == 6

        # String representation
        assert "C" in str(c12)

    def test_get_mass_method(self):
        """Test ElementInfo.get_mass method"""
        c = pt.ELEMENT_LOOKUP["C"]

        # Monoisotopic mass
        mono_mass = c.get_mass(monoisotopic=True)
        assert mono_mass == pytest.approx(12.0, abs=0.1)

        # Average mass (should be slightly higher due to C-13)
        avg_mass = c.get_mass(monoisotopic=False)
        assert avg_mass > mono_mass
        assert avg_mass == pytest.approx(12.011, abs=0.001)


class TestElementLookupComprehensive:
    """Comprehensive integration tests for ElementLookup"""

    def test_all_common_elements_accessible(self):
        """Test that all common biological elements are accessible"""
        common_elements = ["H", "C", "N", "O", "P", "S"]

        for elem in common_elements:
            info = pt.ELEMENT_LOOKUP[elem]
            assert info.symbol == elem or info.symbol in [
                "H",
                "D",
                "T",
            ]  # Handle hydrogen isotopes
            assert info.mass > 0
            assert info.average_mass > 0

    def test_heavy_isotope_labeling_elements(self):
        """Test that common heavy isotope labels are accessible"""
        heavy_isotopes = ["13C", "15N", "18O", "2H", "34S"]

        for iso in heavy_isotopes:
            info = pt.ELEMENT_LOOKUP[iso]
            assert info.mass > 0
            # These should not be monoisotopic (except maybe for special cases)

    def test_selenium_for_selenocysteine(self):
        """Test selenium (used in selenocysteine) is accessible"""
        se = pt.ELEMENT_LOOKUP["Se"]
        assert se.symbol == "Se"
        assert se.number == 34  # Atomic number of selenium
        assert se.mass > 0

    def test_iterator_functionality(self):
        """Test that ElementLookup is iterable"""
        count = 0
        for elem in pt.ELEMENT_LOOKUP:
            assert isinstance(elem, pt.ElementInfo)
            count += 1

        # Should have many elements
        assert count > 100

    def test_lookup_consistency(self):
        """Test that different lookup methods return consistent results"""
        # These should all return the same isotope
        c12_by_tuple = pt.ELEMENT_LOOKUP[("C", 12)]
        c12_by_string = pt.ELEMENT_LOOKUP["12C"]

        assert c12_by_tuple.symbol == c12_by_string.symbol
        assert c12_by_tuple.mass_number == c12_by_string.mass_number
        assert c12_by_tuple.mass == c12_by_string.mass
        assert c12_by_tuple.abundance == c12_by_string.abundance

    def test_monoisotopic_lookup_consistency(self):
        """Test that monoisotopic lookups are consistent"""
        # These should all return monoisotopic carbon (C-12)
        c_by_symbol = pt.ELEMENT_LOOKUP["C"]
        c_by_tuple = pt.ELEMENT_LOOKUP[("C", None)]

        assert c_by_symbol.mass_number == c_by_tuple.mass_number
        assert c_by_symbol.mass == c_by_tuple.mass

    def test_mass_method_consistency(self):
        """Test that mass method returns consistent values"""
        # Direct access
        c = pt.ELEMENT_LOOKUP["C"]

        # Via mass method
        mass_via_method = pt.ELEMENT_LOOKUP.mass("C", monoisotopic=True)

        assert c.mass == mass_via_method

    def test_common_element_atomic_numbers(self):
        """Test that common elements have correct atomic numbers"""
        expected = {
            "H": 1,
            "C": 6,
            "N": 7,
            "O": 8,
            "P": 15,
            "S": 16,
        }

        for symbol, expected_number in expected.items():
            elem = pt.ELEMENT_LOOKUP[symbol]
            assert elem.number == expected_number

    def test_isotope_mass_differences(self):
        """Test that isotope masses differ by approximately neutron mass"""
        c12 = pt.ELEMENT_LOOKUP["12C"]
        c13 = pt.ELEMENT_LOOKUP["13C"]

        # Difference should be approximately 1 neutron mass
        mass_diff = c13.mass - c12.mass
        assert 1.0 < mass_diff < 1.1  # Approximately 1 Da

    def test_deuterium_mass(self):
        """Test that deuterium has approximately correct mass"""
        d = pt.ELEMENT_LOOKUP["D"]

        # Deuterium mass should be around 2.014
        assert 2.01 < d.mass < 2.02

    def test_elements_list_contains_biology_elements(self):
        """Test that get_elements includes all biology-relevant elements"""
        elements = pt.ELEMENT_LOOKUP.get_elements()

        required = ["H", "C", "N", "O", "P", "S", "Se"]
        for elem in required:
            assert elem in elements

    def test_to_dict_method(self):
        """Test ElementInfo.to_dict() conversion"""
        c12 = pt.ELEMENT_LOOKUP["12C"]

        d = c12.to_dict()

        assert isinstance(d, dict)
        assert "number" in d
        assert "symbol" in d
        assert "mass_number" in d
        assert "mass" in d
        assert "abundance" in d
        assert "average_mass" in d

        assert d["number"] == 6
        assert d["symbol"] == "C"
        assert d["mass_number"] == 12

    def test_element_update_method(self):
        """Test ElementInfo.update() method"""
        c12 = pt.ELEMENT_LOOKUP["12C"]

        # Update abundance
        updated = c12.update(abundance=0.5)

        # Original should be unchanged
        assert c12.abundance != 0.5

        # Updated should have new value
        assert updated.abundance == 0.5

        # Other fields should be unchanged
        assert updated.symbol == c12.symbol
        assert updated.mass == c12.mass
        assert updated.mass_number == c12.mass_number
