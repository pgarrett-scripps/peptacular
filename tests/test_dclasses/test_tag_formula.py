"""
Tests for parsing formula modification tags.
"""

import pytest

import peptacular as pt


class TestChargedFormula:
    """Tests for parsing formula modifications"""

    def test_simple_formula(self):
        """Test parsing simple formula"""
        result = pt.ModificationTags.from_string("Formula:C2H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[1].element == pt.Element.H
        assert result.formula[1].occurance == 6
        assert result.charge is None

    def test_single_element(self):
        """Test parsing single element formula"""
        result = pt.ModificationTags.from_string("Formula:O").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 1
        assert result.formula[0].element == pt.Element.O
        assert result.formula[0].occurance == 1

    def test_formula_with_negative_count(self):
        """Test parsing formula with negative element counts"""
        result = pt.ModificationTags.from_string("Formula:H-2O-1").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == pt.Element.H
        assert result.formula[0].occurance == -2
        assert result.formula[1].element == pt.Element.O
        assert result.formula[1].occurance == -1

    def test_formula_with_charge(self):
        """Test parsing formula with charge state"""
        result = pt.ModificationTags.from_string("Formula:C2H6:z+2").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.charge == 2

    def test_formula_with_negative_charge(self):
        """Test parsing formula with negative charge"""
        result = pt.ModificationTags.from_string("Formula:O:z-1").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert result.charge == -1

    def test_formula_with_isotope(self):
        """Test parsing formula with isotope specification [13C2] (count inside bracket)"""
        result = pt.ModificationTags.from_string("Formula:[13C2]H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        assert result.formula[1].element == pt.Element.H
        assert result.formula[1].occurance == 6

    def test_complex_formula(self):
        """Test parsing complex formula"""
        result = pt.ModificationTags.from_string("Formula:C10H15N3O6S").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 5
        # Check carbon
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 10

    def test_formula_with_spaces(self):
        """Test parsing formula with spaces between element pairs (ProForma Rule 1)"""
        result = pt.ModificationTags.from_string("Formula:C12 H20 O2").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 3
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 12
        assert result.formula[1].element == pt.Element.H
        assert result.formula[1].occurance == 20
        assert result.formula[2].element == pt.Element.O
        assert result.formula[2].occurance == 2

    def test_formula_isotope_prefix_notation(self):
        """Test parsing formula with isotope prefix notation [13C2] (ProForma Rule 3)"""
        result = pt.ModificationTags.from_string("Formula:[13C2]H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        assert result.formula[1].element == pt.Element.H
        assert result.formula[1].occurance == 6

    def test_formula_isotope_single_count(self):
        """Test parsing isotope with default count of 1: [13C]"""
        result = pt.ModificationTags.from_string("Formula:[13C]H6").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 1
        assert result.formula[0].isotope == 13

    def test_formula_multiple_isotopes(self):
        """Test parsing formula with multiple isotope specifications [13C2][12C-2]H2N"""
        result = pt.ModificationTags.from_string("Formula:[13C2][12C-2]H2N").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 4
        # First carbon: 13C with count 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        # Second carbon: 12C with count -2
        assert result.formula[1].element == pt.Element.C
        assert result.formula[1].occurance == -2
        assert result.formula[1].isotope == 12
        # Hydrogen
        assert result.formula[2].element == pt.Element.H
        assert result.formula[2].occurance == 2
        # Nitrogen
        assert result.formula[3].element == pt.Element.N
        assert result.formula[3].occurance == 1

    def test_formula_isotope_replacement_example(self):
        """Test ProForma spec example: [13C2]C-2H2N (2 12C replaced by 2 13C)"""
        result = pt.ModificationTags.from_string("Formula:[13C2]C-2H2N").tags[0]
        assert isinstance(result, pt.ChargedFormula)
        assert len(result.formula) == 4
        # 13C with count 2
        assert result.formula[0].element == pt.Element.C
        assert result.formula[0].occurance == 2
        assert result.formula[0].isotope == 13
        # Natural C with count -2
        assert result.formula[1].element == pt.Element.C
        assert result.formula[1].occurance == -2
        assert result.formula[1].isotope is None

    def test_formula_zero_cardinality_not_allowed(self):
        """Test that zero cardinality raises error (ProForma Rule 2)"""
        with pytest.raises(ValueError):
            pt.ModificationTags.from_string("Formula:C0H2").tags[0]
