import json
from pathlib import Path
import warnings

import pytest

from peptacular.mzpaf.official_mzpaf import (
    parse_annotation, Unannotated, MassError, SMILESAnnotation,
    NamedCompoundIonAnnotation, FormulaAnnotation, PeptideFragmentIonAnnotation,
    ReferenceIonAnnotation, IonAnnotationBase)


SPEC_PATH = Path(__file__).parent / "../../../specification/"
SCHEMA_PATH = SPEC_PATH / "annotation-schema.json"

HAS_SCHEMA = SCHEMA_PATH.exists()

if HAS_SCHEMA:
    SCHEMA = json.load(SCHEMA_PATH.open('rt'))
else:
    SCHEMA = None


def validate_schema(annotation: IonAnnotationBase):
    """Helper function to validate annotation against JSONSchema."""
    PASS


def test_parse_unannotated():
    x, = parse_annotation("?^2")
    assert isinstance(x, Unannotated)
    validate_schema(x)


def test_satellite_ion_series():
    base = "da32"
    parsed = parse_annotation(base)[0]
    assert parsed.series == "da"
    assert parsed.position == 32


def test_parse_basic_fragment():
    base = "b14"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'b'
    assert parsed.position == 14
    validate_schema(parsed)


def test_parse_with_neutral_losses():
    base = "b14-H2O-NH3+[Foo]"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'b'
    assert parsed.position == 14
    assert parsed.neutral_losses == ['-H2O', '-NH3', '[Foo]']
    validate_schema(parsed)


def test_parse_with_isotope():
    base = "b14+2i"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'b'
    assert parsed.position == 14
    assert parsed.isotope[0] == 2
    validate_schema(parsed)


def test_parse_with_charge_adducts():
    base = "b14[M+NH4]^2"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'b'
    assert parsed.position == 14
    assert parsed.charge == 2
    assert parsed.adducts == ["M", "NH4"]
    validate_schema(parsed)


def test_parse_with_mass_error():
    base = "b14/0.5ppm"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'b'
    assert parsed.position == 14
    assert parsed.mass_error == MassError(0.5, 'ppm')
    validate_schema(parsed)


def test_parse_with_analyte_reference():
    base = "2@b14"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'b'
    assert parsed.position == 14
    assert parsed.analyte_reference == 2
    validate_schema(parsed)


def test_parse_with_confidence():
    base = "b14*0.05"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'b'
    assert parsed.position == 14
    assert parsed.confidence == 0.05
    validate_schema(parsed)


def test_parse_auxiliary_flag():
    base = "&b14"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'b'
    assert parsed.position == 14
    assert parsed.is_auxiliary
    validate_schema(parsed)


def test_parse_unannotated_labeled():
    base = "?17"
    parsed = parse_annotation(base)[0]
    assert isinstance(parsed, Unannotated)
    assert parsed.unannotated_label == '17'
    validate_schema(parsed)


def test_parse_smiles():
    base = "s{CCC(=O)O}"
    parsed = parse_annotation(base)[0]
    assert isinstance(parsed, SMILESAnnotation)
    assert parsed.smiles == "CCC(=O)O"
    validate_schema(parsed)


def test_parse_external():
    base = "_{foobar}"
    parsed = parse_annotation(base)[0]
    assert isinstance(parsed, NamedCompoundIonAnnotation)
    assert parsed.compound_name == 'foobar'
    validate_schema(parsed)


def test_parse_formula():
    base = "f{C34H53N7O15}"
    parsed = parse_annotation(base)[0]
    assert isinstance(parsed, FormulaAnnotation)
    assert parsed.formula == "C34H53N7O15"
    validate_schema(parsed)


def test_parse_ordinal_with_sequence():
    x, = parse_annotation("y1{{Glycan:Hex1}PEPTIDE}-H2O")
    validate_schema(x)
    assert x.sequence == "{Glycan:Hex1}PEPTIDE"
    assert x.neutral_losses == ["-H2O"]
    assert isinstance(x, PeptideFragmentIonAnnotation)


def test_parse_reference():
    x, = parse_annotation("r[TMT126]")
    assert isinstance(x, ReferenceIonAnnotation)
    assert x.reference == 'TMT126'
    validate_schema(x)


def test_isotope_variants_average():
    average = "b14+2iA"
    parsed = parse_annotation(average)[0]
    assert parsed.isotope[0].average
    assert parsed.isotope[0].isotope == 2


def test_isotope_variants_element():
    average = "b14+2i13C"
    parsed = parse_annotation(average)[0]
    assert parsed.isotope[0].element == 'C'
    assert parsed.isotope[0].isotope == 2
    assert parsed.isotope[0].nucleon_count == 13


def test_roundtrip_json():
    base = "b14-H2O+2i[M+NH4]^2/0.5ppm"
    parsed = parse_annotation(base)[0]
    dup = parsed.from_json(parsed.to_json())
    assert parsed == dup


def test_string_equality():
    base = "2@b14-H2O-NH3+[Foo]+2i[M+NH4]^2/0.5ppm*0.05"
    parsed = parse_annotation(base)[0]
    assert parsed == base


def test_y_ion_series():
    base = "y7"
    parsed = parse_annotation(base)[0]
    assert parsed.series == 'y'
    assert parsed.position == 7
    validate_schema(parsed)


def test_multiple_neutral_losses():
    base = "b5-H2O-CO-NH3"
    parsed = parse_annotation(base)[0]
    assert parsed.neutral_losses == ['-H2O', '-CO', '-NH3']
    validate_schema(parsed)