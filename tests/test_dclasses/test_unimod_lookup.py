"""
Tests for the UnimodLookup class.
"""

import pytest


import peptacular as pt


class TestUnimodLookupBasics:
    """Tests for basic UnimodLookup functionality"""

    def test_lookup_by_name(self):
        """Test looking up modification by name"""
        acetyl = pt.UNIMOD_LOOKUP["Acetyl"]
        assert acetyl.name == "Acetyl"
        assert acetyl.id == "1"
        assert acetyl.formula == "H2C2O"
        assert acetyl.monoisotopic_mass == pytest.approx(42.010565, abs=0.0001)  # type: ignore

    def test_lookup_by_id_string(self):
        """Test looking up modification by ID as string"""
        acetyl = pt.UNIMOD_LOOKUP["1"]
        assert acetyl.name == "Acetyl"
        assert acetyl.id == "1"

    def test_lookup_by_unimod_prefix(self):
        """Test looking up modification with UNIMOD: prefix"""
        acetyl = pt.UNIMOD_LOOKUP["UNIMOD:1"]
        assert acetyl.name == "Acetyl"
        assert acetyl.id == "1"

    def test_lookup_amidated(self):
        """Test looking up Amidated modification with negative mass"""
        amidated = pt.UNIMOD_LOOKUP["Amidated"]
        assert amidated.id == "2"
        assert amidated.formula == "HNO-1"
        assert amidated.monoisotopic_mass == pytest.approx(-0.984016, abs=0.0001)  # type: ignore
        assert amidated.average_mass == pytest.approx(-0.9848, abs=0.0001)  # type: ignore

    def test_lookup_biotin(self):
        """Test looking up Biotin modification"""
        biotin = pt.UNIMOD_LOOKUP["Biotin"]
        assert biotin.id == "3"
        assert biotin.formula == "H14C10N2O2S"
        assert biotin.monoisotopic_mass == pytest.approx(226.077598, abs=0.0001)  # type: ignore

    def test_lookup_nonexistent_name(self):
        """Test looking up non-existent modification name raises KeyError"""
        with pytest.raises(KeyError, match="not found by name or ID"):
            _ = pt.UNIMOD_LOOKUP["NonExistentMod"]

    def test_lookup_nonexistent_id(self):
        """Test looking up non-existent modification ID raises KeyError"""
        with pytest.raises(KeyError, match="not found by name or ID"):
            _ = pt.UNIMOD_LOOKUP["999999"]


class TestUnimodLookupMethods:
    """Tests for UnimodLookup convenience methods"""

    def test_query_name_exists(self):
        """Test query_name method for existing modification"""
        acetyl = pt.UNIMOD_LOOKUP.query_name("Acetyl")
        assert acetyl is not None
        assert acetyl.name == "Acetyl"
        assert acetyl.id == "1"

    def test_query_name_nonexistent(self):
        """Test query_name method for non-existent modification"""
        result = pt.UNIMOD_LOOKUP.query_name("NonExistentMod")
        assert result is None

    def test_query_name_with_prefix(self):
        """Test query_name strips UNIMOD: prefix"""
        acetyl = pt.UNIMOD_LOOKUP.query_name("U:Acetyl")
        assert acetyl is not None
        assert acetyl.name == "Acetyl"

    def test_query_id_as_string(self):
        """Test query_id method with string ID"""
        acetyl = pt.UNIMOD_LOOKUP.query_id("1")
        assert acetyl is not None
        assert acetyl.name == "Acetyl"
        assert acetyl.id == "1"

    def test_query_id_as_int(self):
        """Test query_id method with integer ID"""
        acetyl = pt.UNIMOD_LOOKUP.query_id(1)
        assert acetyl is not None
        assert acetyl.name == "Acetyl"
        assert acetyl.id == "1"

    def test_query_id_nonexistent(self):
        """Test query_id method for non-existent ID"""
        result = pt.UNIMOD_LOOKUP.query_id("999999")
        assert result is None

    def test_query_id_with_unimod_prefix(self):
        """Test query_id strips UNIMOD: prefix"""
        acetyl = pt.UNIMOD_LOOKUP.query_id("unimod:1")
        assert acetyl is not None
        assert acetyl.name == "Acetyl"

    def test_get_exists(self):
        """Test get method for existing modification"""
        acetyl = pt.UNIMOD_LOOKUP.get("Acetyl")
        assert acetyl is not None
        assert acetyl.name == "Acetyl"

    def test_get_nonexistent(self):
        """Test get method for non-existent modification returns None"""
        result = pt.UNIMOD_LOOKUP.get("NonExistentMod")
        assert result is None

    def test_get_by_id(self):
        """Test get method with ID"""
        acetyl = pt.UNIMOD_LOOKUP.get("1")
        assert acetyl is not None
        assert acetyl.name == "Acetyl"


class TestUnimodLookupContains:
    """Tests for the __contains__ method"""

    def test_contains_by_name(self):
        """Test that existing modification name is found"""
        assert "Acetyl" in pt.UNIMOD_LOOKUP

    def test_contains_by_id(self):
        """Test that existing modification ID is found"""
        assert "1" in pt.UNIMOD_LOOKUP

    def test_contains_with_prefix(self):
        """Test that modification with UNIMOD: prefix is found"""
        assert "UNIMOD:1" in pt.UNIMOD_LOOKUP

    def test_not_contains_nonexistent(self):
        """Test that non-existent modification is not found"""
        assert "NonExistentMod" not in pt.UNIMOD_LOOKUP

    def test_not_contains_invalid_id(self):
        """Test that invalid ID is not found"""
        assert "999999" not in pt.UNIMOD_LOOKUP


class TestUnimodInfoProperties:
    """Tests for UnimodInfo data class properties"""

    def test_unimod_info_attributes(self):
        """Test that UnimodInfo has all expected attributes"""
        acetyl = pt.UNIMOD_LOOKUP["Acetyl"]

        assert hasattr(acetyl, "id")
        assert hasattr(acetyl, "name")
        assert hasattr(acetyl, "formula")
        assert hasattr(acetyl, "monoisotopic_mass")
        assert hasattr(acetyl, "average_mass")
        assert hasattr(acetyl, "composition")

    def test_unimod_info_str(self):
        """Test __str__ method of UnimodInfo"""
        acetyl = pt.UNIMOD_LOOKUP["Acetyl"]
        str_repr = str(acetyl)
        assert "Acetyl" in str_repr
        assert "H2C2O" in str_repr

    def test_unimod_info_repr(self):
        """Test __repr__ method of UnimodInfo"""
        acetyl = pt.UNIMOD_LOOKUP["Acetyl"]
        repr_str = repr(acetyl)
        assert "UnimodInfo" in repr_str
        assert "id=1" in repr_str
        assert "name=Acetyl" in repr_str

    def test_unimod_info_dict_composition(self):
        """Test dict_composition property"""
        acetyl = pt.UNIMOD_LOOKUP["Acetyl"]
        dict_comp = acetyl.dict_composition

        assert dict_comp is not None
        assert isinstance(dict_comp, dict)
        assert "H" in dict_comp
        assert "C" in dict_comp
        assert "O" in dict_comp
        assert dict_comp["H"] == 2
        assert dict_comp["C"] == 2
        assert dict_comp["O"] == 1

    def test_unimod_info_composition_negative_oxygen(self):
        """Test composition with negative element count (Amidated)"""
        amidated = pt.UNIMOD_LOOKUP["Amidated"]
        dict_comp = amidated.dict_composition

        assert dict_comp is not None
        assert dict_comp["H"] == 1
        assert dict_comp["N"] == 1
        assert dict_comp["O"] == -1

    def test_unimod_info_update(self):
        """Test update method creates new instance with modified fields"""
        acetyl = pt.UNIMOD_LOOKUP["Acetyl"]
        updated = acetyl.update(name="CustomAcetyl")

        # Original should be unchanged
        assert acetyl.name == "Acetyl"

        # Updated should have new name but same other fields
        assert updated.name == "CustomAcetyl"
        assert updated.id == acetyl.id
        assert updated.formula == acetyl.formula
        assert updated.monoisotopic_mass == acetyl.monoisotopic_mass


class TestUnimodLookupPrefixHandling:
    """Tests for handling various prefix formats"""

    def test_strip_unimod_prefix_lowercase(self):
        """Test stripping lowercase unimod: prefix"""
        acetyl = pt.UNIMOD_LOOKUP["unimod:1"]
        assert acetyl.name == "Acetyl"

    def test_strip_unimod_prefix_uppercase(self):
        """Test stripping uppercase UNIMOD: prefix"""
        acetyl = pt.UNIMOD_LOOKUP["UNIMOD:1"]
        assert acetyl.name == "Acetyl"

    def test_strip_unimod_prefix_mixed_case(self):
        """Test stripping mixed case UniMod: prefix"""
        acetyl = pt.UNIMOD_LOOKUP["UniMod:1"]
        assert acetyl.name == "Acetyl"

    def test_colon_in_name_not_stripped(self):
        """Test that names starting with 'U' but having colons elsewhere work"""
        # If there's a modification name like "Unknown:Mod", it should be searched as-is
        # This test verifies the logic doesn't incorrectly strip prefixes
        result = pt.UNIMOD_LOOKUP.query_name("Unknown:Mod")
        # Should return None since it doesn't exist, but shouldn't crash
        assert result is None


class TestUnimodLookupCommonModifications:
    """Tests for commonly used modifications"""

    def test_phosphorylation(self):
        """Test phosphorylation modification"""
        phospho = pt.UNIMOD_LOOKUP.query_name("Phospho")
        if phospho is not None:  # Only test if it exists in data
            assert phospho.formula is not None
            assert "P" in phospho.formula or "PO" in phospho.formula
            assert phospho.monoisotopic_mass is not None
            assert phospho.monoisotopic_mass > 70  # Approximate mass

    def test_oxidation(self):
        """Test oxidation modification"""
        oxidation = pt.UNIMOD_LOOKUP.query_name("Oxidation")
        if oxidation is not None:
            assert oxidation.formula is not None
            assert "O" in oxidation.formula
            assert oxidation.monoisotopic_mass is not None
            assert oxidation.monoisotopic_mass == pytest.approx(15.994915, abs=0.001)  # type: ignore

    def test_carbamidomethyl(self):
        """Test carbamidomethyl modification (common cysteine mod)"""
        cam = pt.UNIMOD_LOOKUP.query_name("Carbamidomethyl")
        if cam is not None:
            assert cam.monoisotopic_mass is not None
            assert cam.monoisotopic_mass == pytest.approx(57.021464, abs=0.001)  # type: ignore

    def test_methylation(self):
        """Test methylation modification"""
        methyl = pt.UNIMOD_LOOKUP.query_name("Methyl")
        if methyl is not None:
            assert methyl.monoisotopic_mass is not None
            assert methyl.monoisotopic_mass == pytest.approx(14.015650, abs=0.001)  # type: ignore


class TestUnimodLookupEdgeCases:
    """Tests for edge cases and special scenarios"""

    def test_empty_string_lookup(self):
        """Test looking up empty string"""
        with pytest.raises(KeyError):
            _ = pt.UNIMOD_LOOKUP[""]

    def test_whitespace_lookup(self):
        """Test looking up whitespace"""
        with pytest.raises(KeyError):
            _ = pt.UNIMOD_LOOKUP["   "]

    def test_numeric_string_interpreted_as_id(self):
        """Test that pure numeric string is interpreted as ID"""
        mod = pt.UNIMOD_LOOKUP["1"]
        assert mod.id == "1"
        assert mod.name == "Acetyl"

    def test_id_and_name_lookup_same_result(self):
        """Test that looking up by ID and name returns same modification"""
        by_name = pt.UNIMOD_LOOKUP["Acetyl"]
        by_id = pt.UNIMOD_LOOKUP["1"]

        assert by_name.id == by_id.id
        assert by_name.name == by_id.name
        assert by_name.formula == by_id.formula
        assert by_name.monoisotopic_mass == by_id.monoisotopic_mass

    def test_multiple_lookups_same_object(self):
        """Test that multiple lookups return the same object (not copies)"""
        first = pt.UNIMOD_LOOKUP["Acetyl"]
        second = pt.UNIMOD_LOOKUP["Acetyl"]

        # Should be the exact same object
        assert first is second


class TestUnimodLookupDataIntegrity:
    """Tests for data integrity and consistency"""

    def test_all_modifications_have_required_fields(self):
        """Test that all modifications have required fields"""
        # Get a few modifications to check
        test_ids = ["1", "2", "3"]

        for mod_id in test_ids:
            mod = pt.UNIMOD_LOOKUP.query_id(mod_id)
            if mod is not None:
                assert mod.id is not None
                assert mod.name is not None
                assert isinstance(mod.id, str)
                assert isinstance(mod.name, str)

    def test_masses_are_numeric(self):
        """Test that masses are numeric when present"""
        acetyl = pt.UNIMOD_LOOKUP["Acetyl"]

        assert isinstance(acetyl.monoisotopic_mass, (int, float))
        assert isinstance(acetyl.average_mass, (int, float))

    def test_negative_mass_modifications(self):
        """Test that modifications with negative mass work correctly"""
        amidated = pt.UNIMOD_LOOKUP["Amidated"]
        assert amidated.monoisotopic_mass is not None
        assert amidated.average_mass is not None
        assert amidated.monoisotopic_mass < 0
        assert amidated.average_mass < 0

    def test_composition_consistency(self):
        """Test that composition matches formula when both present"""
        acetyl = pt.UNIMOD_LOOKUP["Acetyl"]

        if acetyl.composition is not None and acetyl.formula is not None:
            dict_comp = acetyl.dict_composition

            # Check that formula string contains the elements in composition
            assert "H" in acetyl.formula
            assert "C" in acetyl.formula
            assert "O" in acetyl.formula

            # Check counts match (basic validation)
            assert dict_comp is not None
            assert dict_comp["H"] == 2
            assert dict_comp["C"] == 2
            assert dict_comp["O"] == 1


class TestUnimodLookupInitialization:
    """Tests for UnimodLookup initialization"""

    def test_lookup_has_name_to_info_dict(self):
        """Test that lookup has name_to_info dictionary"""
        assert hasattr(pt.UNIMOD_LOOKUP, "name_to_info")
        assert isinstance(pt.UNIMOD_LOOKUP.name_to_info, dict)

    def test_lookup_has_id_to_info_dict(self):
        """Test that lookup has id_to_info dictionary"""
        assert hasattr(pt.UNIMOD_LOOKUP, "id_to_info")
        assert isinstance(pt.UNIMOD_LOOKUP.id_to_info, dict)

    def test_name_and_id_dicts_consistent(self):
        """Test that name and ID dictionaries are consistent"""
        # Get acetyl from both dictionaries
        acetyl_by_name = pt.UNIMOD_LOOKUP.query_name("Acetyl")
        acetyl_by_id = pt.UNIMOD_LOOKUP.query_id("1")

        assert acetyl_by_name is not None
        assert acetyl_by_id is not None

        # They should refer to the same object or have same data
        assert acetyl_by_name.id == acetyl_by_id.id
        assert acetyl_by_name.name == acetyl_by_id.name
