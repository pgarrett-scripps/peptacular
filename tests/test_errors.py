"""Tests for error classes"""

import pytest
from peptacular.errors import (
    AmbiguousAminoAcidError,
    AmbiguousModificationError,
    AmbiguousSequenceError,
    InvalidModificationMassError,
    InvalidSequenceError,
    UnknownModificationError,
    UnknownModificationMassError,
    UnknownAminoAcidError,
    InvalidDeltaMassError,
    InvalidCompositionError,
    DeltaMassCompositionError,
    InvalidChemFormulaError,
    InvalidGlycanFormulaError,
    ProFormaFormatError,
)


class TestErrors:
    """Test custom error classes"""

    def test_ambiguous_amino_acid_error(self):
        """Test AmbiguousAminoAcidError"""
        with pytest.raises(AmbiguousAminoAcidError) as exc_info:
            raise AmbiguousAminoAcidError("X", "Multiple possibilities")
        
        assert exc_info.value.aa == "X"
        assert exc_info.value.msg == "Multiple possibilities"
        assert "Ambiguous amino acid: X" in str(exc_info.value)

    def test_ambiguous_modification_error(self):
        """Test AmbiguousModificationError"""
        with pytest.raises(AmbiguousModificationError):
            raise AmbiguousModificationError("Ambiguous mod")

    def test_ambiguous_sequence_error(self):
        """Test AmbiguousSequenceError"""
        with pytest.raises(AmbiguousSequenceError):
            raise AmbiguousSequenceError("Ambiguous sequence")

    def test_invalid_modification_mass_error(self):
        """Test InvalidModificationMassError"""
        with pytest.raises(InvalidModificationMassError) as exc_info:
            raise InvalidModificationMassError("BadMod")
        
        assert exc_info.value.modification_str == "BadMod"
        assert "Cannot determine mass" in str(exc_info.value)

    def test_invalid_sequence_error(self):
        """Test InvalidSequenceError"""
        with pytest.raises(InvalidSequenceError):
            raise InvalidSequenceError("Invalid seq")

    def test_unknown_modification_error(self):
        """Test UnknownModificationError"""
        with pytest.raises(UnknownModificationError) as exc_info:
            raise UnknownModificationError("UnknownMod")
        
        assert exc_info.value.modification == "UnknownMod"
        assert "Unknown modification" in str(exc_info.value)

    def test_unknown_modification_mass_error(self):
        """Test UnknownModificationMassError"""
        with pytest.raises(UnknownModificationMassError) as exc_info:
            raise UnknownModificationMassError("123.456")
        
        assert exc_info.value.mass == "123.456"
        assert "Unknown modification" in str(exc_info.value)

    def test_unknown_amino_acid_error(self):
        """Test UnknownAminoAcidError"""
        with pytest.raises(UnknownAminoAcidError) as exc_info:
            raise UnknownAminoAcidError("Z")
        
        assert exc_info.value.amino_acid == "Z"
        assert "Unknown amino acid" in str(exc_info.value)

    def test_invalid_delta_mass_error(self):
        """Test InvalidDeltaMassError"""
        with pytest.raises(InvalidDeltaMassError) as exc_info:
            raise InvalidDeltaMassError("abc")
        
        assert exc_info.value.mass == "abc"
        assert "Invalid delta mass" in str(exc_info.value)

    def test_invalid_composition_error(self):
        """Test InvalidCompositionError"""
        with pytest.raises(InvalidCompositionError) as exc_info:
            raise InvalidCompositionError("H2O3X")
        
        assert exc_info.value.composition == "H2O3X"
        assert "Cannot retrieve composition" in str(exc_info.value)

    def test_delta_mass_composition_error(self):
        """Test DeltaMassCompositionError"""
        with pytest.raises(DeltaMassCompositionError) as exc_info:
            raise DeltaMassCompositionError("BadComp")
        
        assert exc_info.value.composition == "BadComp"
        assert "Cannot retrieve composition" in str(exc_info.value)

    def test_invalid_chem_formula_error(self):
        """Test InvalidChemFormulaError"""
        with pytest.raises(InvalidChemFormulaError) as exc_info:
            raise InvalidChemFormulaError("H2X3", "Unknown element X")
        
        assert exc_info.value.formula == "H2X3"
        assert exc_info.value.msg == "Unknown element X"
        assert "Error parsing chem formula" in str(exc_info.value)

    def test_invalid_glycan_formula_error(self):
        """Test InvalidGlycanFormulaError"""
        with pytest.raises(InvalidGlycanFormulaError) as exc_info:
            raise InvalidGlycanFormulaError("{Hex}X", "Invalid syntax")
        
        assert exc_info.value.formula == "{Hex}X"
        assert exc_info.value.msg == "Invalid syntax"
        assert "Error parsing glycan formula" in str(exc_info.value)

    def test_proforma_format_error_within_bounds(self):
        """Test ProFormaFormatError with valid index"""
        with pytest.raises(ProFormaFormatError) as exc_info:
            raise ProFormaFormatError("Invalid character", 5, "PEPTIDE")
        
        assert exc_info.value.msg == "Invalid character"
        assert ">>>D<<<" in str(exc_info.value)  # Index 5 is 'D'
        assert "At index: 5" in str(exc_info.value)

    def test_proforma_format_error_out_of_bounds(self):
        """Test ProFormaFormatError with out of bounds index"""
        with pytest.raises(ProFormaFormatError) as exc_info:
            raise ProFormaFormatError("Invalid character", 100, "PEPTIDE")
        
        # Should not crash, just use sequence as-is
        assert "PEPTIDE" in str(exc_info.value)
        assert "At index: 100" in str(exc_info.value)

    def test_proforma_format_error_negative_index(self):
        """Test ProFormaFormatError with negative index"""
        with pytest.raises(ProFormaFormatError) as exc_info:
            raise ProFormaFormatError("Invalid character", -1, "PEPTIDE")
        
        # Negative index is out of bounds
        assert "PEPTIDE" in str(exc_info.value)
