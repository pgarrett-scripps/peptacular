"""Tests for the modern mzPAF implementation (simple.py)."""

import pytest
from peptacular.mzpaf.frag_annot import (
    MzpafAnnotationParser,
    FragmentAnnotation,
    PeptideIon,
    InternalFragment,
    PrecursorIon,
    ImmoniumIon,
    ReferenceIon,
    FormulaIon,
    SMILESIon,
    NamedCompound,
    UnknownIon,
    NeutralLoss,
    IsotopeLabel,
    MassError,
    parse_mzpaf_annotation,
    format_mzpaf_annotation,
)


class TestPeptideIonParsing:
    """Test parsing of standard peptide fragment ions."""
    
    def test_simple_b_ion(self):
        result = parse_mzpaf_annotation("b2")
        assert len(result) == 1
        assert isinstance(result[0].ion_type, PeptideIon)
        assert result[0].ion_type.series == "b"
        assert result[0].ion_type.position == 2
        assert str(result[0]) == "b2"
    
    def test_simple_y_ion(self):
        result = parse_mzpaf_annotation("y5")
        assert len(result) == 1
        assert result[0].ion_type.series == "y"
        assert result[0].ion_type.position == 5
        assert str(result[0]) == "y5"
    
    def test_peptide_ion_with_sequence(self):
        result = parse_mzpaf_annotation("b2{LC[Carbamidomethyl]}")
        assert result[0].ion_type.sequence == "LC[Carbamidomethyl]"
        assert str(result[0]) == "b2{LC[Carbamidomethyl]}"
    
    def test_peptide_ion_with_neutral_loss(self):
        result = parse_mzpaf_annotation("y5-H2O")
        assert len(result[0].neutral_losses) == 1
        assert result[0].neutral_losses[0].name == "H2O"
        assert str(result[0]) == "y5-H2O"
    
    def test_peptide_ion_multiple_neutral_losses(self):
        result = parse_mzpaf_annotation("b7-H2O-NH3")
        assert len(result[0].neutral_losses) == 2
        assert str(result[0]) == "b7-H2O-NH3"


class TestInternalFragmentParsing:
    """Test parsing of internal fragment ions."""
    
    def test_simple_internal(self):
        result = parse_mzpaf_annotation("m3:6")
        assert isinstance(result[0].ion_type, InternalFragment)
        assert result[0].ion_type.start == 3
        assert result[0].ion_type.end == 6
        assert str(result[0]) == "m3:6"
    
    def test_internal_with_backbone(self):
        result = parse_mzpaf_annotation("mxa3:6")
        assert result[0].ion_type.n_terminal_cleavage == "x"
        assert result[0].ion_type.c_terminal_cleavage == "a"
        assert str(result[0]) == "mxa3:6"
    
    def test_internal_with_sequence(self):
        result = parse_mzpaf_annotation("m5:8{PEPT}")
        assert result[0].ion_type.sequence == "PEPT"
        assert str(result[0]) == "m5:8{PEPT}"


class TestPrecursorIonParsing:
    """Test parsing of precursor ions."""
    
    def test_simple_precursor(self):
        result = parse_mzpaf_annotation("p")
        assert isinstance(result[0].ion_type, PrecursorIon)
        assert str(result[0]) == "p"
    
    def test_precursor_with_neutral_loss(self):
        result = parse_mzpaf_annotation("p-NH3")
        assert len(result[0].neutral_losses) == 1
        assert str(result[0]) == "p-NH3"
    
    def test_precursor_with_sequence(self):
        result = parse_mzpaf_annotation("p{PEPTIDE}")
        assert result[0].ion_type.sequence == "PEPTIDE"
        assert str(result[0]) == "p{PEPTIDE}"


class TestImmoniumIonParsing:
    """Test parsing of immonium ions."""
    
    def test_simple_immonium(self):
        result = parse_mzpaf_annotation("IA")
        assert isinstance(result[0].ion_type, ImmoniumIon)
        assert result[0].ion_type.amino_acid == "A"
        assert str(result[0]) == "IA"
    
    def test_immonium_with_modification(self):
        result = parse_mzpaf_annotation("IA[Oxidation]")
        assert result[0].ion_type.modification == "Oxidation"
        assert str(result[0]) == "IA[Oxidation]"


class TestReferenceIonParsing:
    """Test parsing of reference ions."""
    
    def test_reference_ion(self):
        result = parse_mzpaf_annotation("r[TMT126]")
        assert isinstance(result[0].ion_type, ReferenceIon)
        assert result[0].ion_type.name == "TMT126"
        assert str(result[0]) == "r[TMT126]"


class TestFormulaIonParsing:
    """Test parsing of formula ions."""
    
    def test_formula_ion(self):
        result = parse_mzpaf_annotation("f{C2H5NO}")
        assert isinstance(result[0].ion_type, FormulaIon)
        assert result[0].ion_type.formula == "C2H5NO"
        assert str(result[0]) == "f{C2H5NO}"


class TestSMILESIonParsing:
    """Test parsing of SMILES ions."""
    
    def test_smiles_ion(self):
        result = parse_mzpaf_annotation("s{CC(C)C}")
        assert isinstance(result[0].ion_type, SMILESIon)
        assert result[0].ion_type.smiles == "CC(C)C"
        assert str(result[0]) == "s{CC(C)C}"


class TestNamedCompoundParsing:
    """Test parsing of named compound ions."""
    
    def test_named_compound(self):
        result = parse_mzpaf_annotation("_{lipid_a}")
        assert isinstance(result[0].ion_type, NamedCompound)
        assert result[0].ion_type.name == "lipid_a"
        assert str(result[0]) == "_{lipid_a}"


class TestUnannotatedParsing:
    """Test parsing of unannotated/unknown ions."""
    
    def test_simple_unannotated(self):
        result = parse_mzpaf_annotation("?")
        assert isinstance(result[0].ion_type, UnknownIon)
        assert result[0].ion_type.label is None
        assert str(result[0]) == "?"
    
    def test_unannotated_with_label(self):
        result = parse_mzpaf_annotation("?42")
        assert result[0].ion_type.label == 42
        assert str(result[0]) == "?42"


class TestModifiersParsing:
    """Test parsing of various modifiers (isotopes, adducts, charge, etc.)."""
    
    def test_isotope_simple(self):
        result = parse_mzpaf_annotation("b2+i")
        assert len(result[0].isotopes) == 1
        assert result[0].isotopes[0].offset == 1
        assert str(result[0]) == "b2+i"
    
    def test_isotope_with_element(self):
        result = parse_mzpaf_annotation("y5+i13C")
        assert result[0].isotopes[0].element == "C"
        assert result[0].isotopes[0].nucleon_count == 13
        assert str(result[0]) == "y5+i13C"
    
    def test_charge_state(self):
        result = parse_mzpaf_annotation("b5^2")
        assert result[0].charge == 2
        assert str(result[0]) == "b5^2"
    
    def test_analyte_reference(self):
        result = parse_mzpaf_annotation("2@b5")
        assert result[0].analyte_reference == 2
        assert str(result[0]) == "2@b5"
    
    def test_mass_error_da(self):
        result = parse_mzpaf_annotation("y3/0.5")
        assert result[0].mass_error.value == 0.5
        assert result[0].mass_error.unit == "Da"
        assert str(result[0]) == "y3/0.5"
    
    def test_mass_error_ppm(self):
        result = parse_mzpaf_annotation("b2/0.5ppm")
        assert result[0].mass_error.value == 0.5
        assert result[0].mass_error.unit == "ppm"
        assert str(result[0]) == "b2/0.5ppm"
    
    def test_confidence(self):
        result = parse_mzpaf_annotation("y5*0.95")
        assert result[0].confidence == 0.95
        assert str(result[0]) == "y5*0.95"
    
    def test_auxiliary_annotation(self):
        result = parse_mzpaf_annotation("&y3")
        assert result[0].is_auxiliary is True
        # Note: auxiliary prefix is added in __str__
        assert str(result[0]) == "&y3"


class TestComplexAnnotations:
    """Test complex annotations with multiple modifiers."""
    
    def test_full_annotation(self):
        result = parse_mzpaf_annotation("2@b5+i^2/0.5ppm*0.95")
        assert result[0].analyte_reference == 2
        assert result[0].ion_type.position == 5
        assert len(result[0].isotopes) == 1
        assert result[0].charge == 2
        assert result[0].mass_error.value == 0.5
        assert result[0].confidence == 0.95
        assert str(result[0]) == "2@b5+i^2/0.5ppm*0.95"
    
    def test_peptide_with_sequence_and_loss(self):
        result = parse_mzpaf_annotation("b2{LC[Carbamidomethyl]}-H2O")
        assert result[0].ion_type.sequence == "LC[Carbamidomethyl]"
        assert len(result[0].neutral_losses) == 1
        assert str(result[0]) == "b2{LC[Carbamidomethyl]}-H2O"


class TestMultipleAnnotations:
    """Test parsing multiple comma-separated annotations."""
    
    def test_two_annotations(self):
        result = parse_mzpaf_annotation("b2*0.6,y3*0.4")
        assert len(result) == 2
        assert result[0].confidence == 0.6
        assert result[1].confidence == 0.4
    
    def test_confidence_validation(self):
        # Total confidence should not exceed 1.0
        with pytest.raises(ValueError, match="exceeds 1.0"):
            parse_mzpaf_annotation("b2*0.7,y3*0.5")


class TestRoundTripSerialization:
    """Test that parsing and serializing produces the same string."""
    
    @pytest.mark.parametrize("annotation_str", [
        "b2",
        "y5-H2O",
        "m3:6",
        "mxa3:6",
        "mzc3:6",
        "myc3:6-H2O^2",
        "p-NH3",
        "IA[Oxidation]",
        "r[TMT126]",
        "f{C2H5NO}",
        "_{lipid_a}",
        "s{CC(C)C}",
        "?",
        "?42",
        "2@b5+i^2/0.5ppm*0.95",
        "b2{LC[Carbamidomethyl]}",
    ])
    def test_round_trip(self, annotation_str):
        result = parse_mzpaf_annotation(annotation_str)
        assert len(result) == 1
        reconstructed = str(result[0])
        assert reconstructed == annotation_str


class TestMassCalculation:
    """Test mass calculation functionality."""
    
    def test_formula_mass(self):
        result = parse_mzpaf_annotation("f{C2H5NO}")
        # Should be able to calculate mass
        mass = result[0].neutral_mass()
        assert mass > 0
    
    def test_peptide_mass_with_sequence(self):
        result = parse_mzpaf_annotation("b2{AC}")
        # Should calculate mass from sequence
        mass = result[0].neutral_mass()
        assert mass > 0
    
    def test_peptide_mass_without_sequence_raises(self):
        result = parse_mzpaf_annotation("b2")
        # Should raise without sequence
        with pytest.raises(ValueError, match="Cannot calculate mass without sequence"):
            result[0].neutral_mass()
    
    def test_reference_mass(self):
        result = parse_mzpaf_annotation("r[TMT126]")
        # Should get mass from reference database
        mass = result[0].neutral_mass()
        assert mass > 0
    
    def test_mz_calculation(self):
        result = parse_mzpaf_annotation("b2{AC}^2")
        # Should calculate m/z
        mz = result[0].mz()
        neutral = result[0].neutral_mass()
        # m/z should be less than neutral mass for charge > 1
        assert mz < neutral


class TestErrorHandling:
    """Test error handling and invalid inputs."""
    
    def test_invalid_annotation_raises(self):
        with pytest.raises(ValueError, match="Invalid mzPAF annotation"):
            parse_mzpaf_annotation("invalid")
    
    def test_wrap_errors(self):
        result = parse_mzpaf_annotation("invalid", wrap_errors=True)
        # Should return error annotation instead of raising
        assert len(result) == 1
        assert result[0].is_auxiliary is True
    
    def test_zero_charge_raises(self):
        with pytest.raises(ValueError, match="Charge cannot be zero"):
            parse_mzpaf_annotation("b2^0")
    
    def test_empty_string(self):
        result = parse_mzpaf_annotation("")
        assert result == []


class TestJSONSerialization:
    """Test JSON serialization."""
    
    def test_to_json(self):
        result = parse_mzpaf_annotation("b2{AC}-H2O^2")
        json_data = result[0].to_json()
        
        assert json_data["charge"] == 2
        assert json_data["molecule_description"]["series_label"] == "peptide"
        assert json_data["molecule_description"]["series"] == "b"
        assert json_data["molecule_description"]["position"] == 2
        assert len(json_data["neutral_losses"]) == 1
    
    def test_mass_error_json(self):
        result = parse_mzpaf_annotation("y3/0.5ppm")
        json_data = result[0].to_json()
        
        assert json_data["mass_error"]["value"] == 0.5
        assert json_data["mass_error"]["unit"] == "ppm"


class TestNeutralLoss:
    """Test NeutralLoss parsing and formatting."""
    
    def test_parse_simple_formula(self):
        losses = NeutralLoss.parse_many("+H2O")
        assert len(losses) == 1
        assert losses[0].name == "H2O"
        assert losses[0].coefficient == 1
    
    def test_parse_negative_loss(self):
        losses = NeutralLoss.parse_many("-NH3")
        assert losses[0].coefficient == -1
    
    def test_parse_multiple_losses(self):
        losses = NeutralLoss.parse_many("+H2O-NH3")
        assert len(losses) == 2
        assert losses[0].coefficient == 1
        assert losses[1].coefficient == -1


class TestIsotopeLabel:
    """Test IsotopeLabel parsing and formatting."""
    
    def test_parse_simple_isotope(self):
        iso = IsotopeLabel.from_string("+i")
        assert iso.offset == 1
        assert iso.element is None
    
    def test_parse_element_specific(self):
        iso = IsotopeLabel.from_string("+i13C")
        assert iso.offset == 1
        assert iso.element == "C"
        assert iso.nucleon_count == 13
    
    def test_format_isotope(self):
        iso = IsotopeLabel(offset=1, element="C", nucleon_count=13)
        assert str(iso) == "+i13C"


class TestConvenienceFunctions:
    """Test convenience functions."""
    
    def test_parse_mzpaf_annotation(self):
        result = parse_mzpaf_annotation("b2")
        assert len(result) == 1
        assert isinstance(result[0], FragmentAnnotation)
    
    def test_format_mzpaf_annotation(self):
        result = parse_mzpaf_annotation("y5-H2O")
        formatted = format_mzpaf_annotation(result[0])
        assert formatted == "y5-H2O"
