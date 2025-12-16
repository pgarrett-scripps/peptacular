"""
mzPAF Parser Module - Enhanced with backbone cleavage notation for internal fragments
Parses and serializes mzPAF (Peak Annotation Format) strings for peptide mass spectra.
Based on mzPAF specification v1.0 with extended internal fragment support
"""

import re
from dataclasses import dataclass
from typing import Literal
from enum import Enum


# Table from the specification showing differences from yb
INTERNAL_MASS_DIFFS: dict[tuple[str, str], None | str] = {
    ("a", "x"): None,  #  Default, no difference
    ("b", "x"): "+CO",
    ("c", "x"): "+CHNO",
    ("a", "y"): "-CO",
    ("b", "y"): None,  # Default, no difference
    ("c", "y"): "+NH",
    ("a", "z"): "-CHNO",
    ("b", "z"): "-NH",
    ("c", "z"): None,  # No difference
}


class IonSeries(Enum):
    """Enumeration of ion series types"""

    A = "a"
    B = "b"
    C = "c"
    D = "d"
    V = "v"
    W = "w"
    X = "x"
    Y = "y"
    Z = "z"
    DA = "da"
    DB = "db"
    WA = "wa"
    WB = "wb"


class BackboneCleavageType(Enum):
    """Types of backbone cleavages for internal fragments"""

    A = "a"  # C-CO bond cleavage
    B = "b"  # CO-NH bond cleavage
    C = "c"  # NH-CH bond cleavage
    X = "x"  # CH-CO bond cleavage
    Y = "y"  # CO-NH bond cleavage
    Z = "z"  # NH-CH bond cleavage


@dataclass
class MassError:
    """Represents mass error with value and unit"""

    value: float
    unit: Literal["da", "ppm"] = "da"  # "da" or "ppm"

    def serialize(self) -> str:
        if self.unit == "ppm":
            return f"{self.value:g}ppm"
        elif self.unit == "da":
            return f"{self.value:g}"
        else:
            raise ValueError(f"Unknown mass error unit: {self.unit}")

    def __str__(self) -> str:
        return self.serialize()


@dataclass
class IsotopeSpecification:
    """Represents isotope information"""

    count: int = 0  # number of isotopes above/below monoisotope
    element: str | None = None  # e.g., "13C", "15N"
    is_average: bool = False  # True for averaged isotopomers

    def serialize(self) -> str:
        if self.count == 0:
            return ""

        sign = "+" if self.count > 0 else ""
        count_str = "" if abs(self.count) == 1 else str(abs(self.count))

        if self.is_average is True:
            return f"{sign}{count_str}iA"
        elif self.element is not None:
            return f"{sign}{count_str}i{self.element}"
        else:
            return f"{sign}{count_str}i"

    def __str__(self) -> str:
        return self.serialize()


@dataclass
class PeptideIon:
    """Represents a primary peptide fragment ion"""

    series: str  # b, y, a, x, c, z, etc.
    position: int
    sequence: str | None = None  # ProForma sequence

    def serialize(self) -> str:
        result = f"{self.series}{self.position}"
        if self.sequence:
            result += f"{{{self.sequence}}}"
        return result

    def __str__(self):
        return self.serialize()


@dataclass
class InternalFragment:
    """Represents an internal fragment ion with optional backbone cleavage specification"""

    start_position: int
    end_position: int
    sequence: str | None = None

    def serialize(self) -> str:
        # If using default yb cleavage, just use 'm'
        result = f"m{self.start_position}:{self.end_position}"

        if self.sequence:
            result += f"{{{self.sequence}}}"
        return result

    def __str__(self):
        return self.serialize()


@dataclass
class ImmoniumIon:
    """Represents an immonium ion"""

    amino_acid: str
    modification: str | None = None

    def serialize(self) -> str:
        result = f"I{self.amino_acid}"
        if self.modification:
            result += f"[{self.modification}]"
        return result

    def __str__(self):
        return self.serialize()


@dataclass
class ReferenceIon:
    """Represents a reference ion"""

    name: str

    def serialize(self) -> str:
        return f"r[{self.name}]"

    def __str__(self):
        return self.serialize()


@dataclass
class NamedCompound:
    """Represents a named compound"""

    name: str

    def serialize(self) -> str:
        return f"_{{{self.name}}}"

    def __str__(self):
        return self.serialize()


@dataclass
class ChemicalFormula:
    """Represents a chemical formula"""

    formula: str

    def serialize(self) -> str:
        return f"f{{{self.formula}}}"

    def __str__(self):
        return self.serialize()


@dataclass
class SMILESCompound:
    """Represents a SMILES string"""

    smiles: str

    def serialize(self) -> str:
        return f"s{{{self.smiles}}}"

    def __str__(self):
        return self.serialize()


@dataclass
class UnknownIon:
    """Represents an unknown/unannotated ion"""

    label: int | None = None

    def serialize(self) -> str:
        if self.label is not None:
            return f"?{self.label}"
        return "?"

    def __str__(self):
        return self.serialize()


@dataclass
class PrecursorIon:
    """Represents a precursor ion"""

    def serialize(self) -> str:
        return "p"

    def __str__(self):
        return self.serialize()


@dataclass
class FragmentAnnotation:
    """Complete fragment ion annotation"""

    # Core ion description
    ion_type: (
        PeptideIon
        | InternalFragment
        | ImmoniumIon
        | ReferenceIon
        | NamedCompound
        | ChemicalFormula
        | SMILESCompound
        | UnknownIon
        | PrecursorIon
    )

    # Optional components
    analyte_reference: int | None = None
    is_auxiliary: bool = False
    neutral_losses: list[str] | None = None  # Optional list of neutral losses
    isotope: list[IsotopeSpecification] | None = (
        None  # Optional list of isotope specifications
    )
    adducts: list[str] | None = None  # Optional list of adducts, e.g., ["M+H+Na"]
    charge: int = 1
    mass_error: MassError | None = None
    confidence: float | None = None

    def serialize(self) -> str:
        """Serialize the annotation back to mzPAF string format"""
        parts: list[str] = []

        # Auxiliary marker
        if self.is_auxiliary:
            parts.append("&")

        # Analyte reference
        if self.analyte_reference is not None:
            parts.append(f"{self.analyte_reference}@")

        # Ion type
        parts.append(str(self.ion_type))

        # Neutral losses
        if self.neutral_losses:
            for loss in self.neutral_losses:
                parts.append(loss)

        # Isotope
        if self.isotope:
            for iso in self.isotope:
                if iso.count != 0:
                    parts.append(str(iso))

        # Adduct type
        if self.adducts:
            for adduct in self.adducts:
                parts.append(f"[{adduct}]")

        # Charge state
        if self.charge > 1:
            parts.append(f"^{self.charge}")

        # Mass error
        if self.mass_error:
            parts.append(f"/{self.mass_error}")

        # Confidence
        if self.confidence is not None:
            parts.append(f"*{self.confidence:g}")

        return "".join(parts)

    def __str__(self) -> str:
        return self.serialize()

    @staticmethod
    def parse(annotation_str: str) -> "FragmentAnnotation":
        """Parse a single mzPAF annotation string into a FragmentAnnotation object"""
        parser = mzPAFParser()
        return parser.parse_single(annotation_str)


class mzPAFParser:
    """Parser for mzPAF annotation strings with enhanced internal fragment support"""

    # Updated regex pattern to support backbone cleavage notation for internal fragments
    PATTERN = re.compile(
        r"^(?P<is_auxiliary>&)?"
        r"(?:(?P<analyte_reference>\d+)@)?"
        r"(?:"
        r"(?:(?P<series>(?:da|db|wa|wb)|[axbyczdwv]\.?)(?P<ordinal>\d+)(?:\{(?P<sequence_ordinal>.+)\})?)|"
        r"(?P<series_internal>m(?P<internal_start>\d+):(?P<internal_end>\d+)(?:\{(?P<sequence_internal>.+)\})?)|"
        r"(?P<precursor>p)|"
        r"(?:I(?P<immonium>[A-Z])(?:\[(?P<immonium_modification>(?:[^\]]+))\])?)|"
        r"(?P<reference>r(?:(?:\[(?P<reference_label>[^\]]+)\])))|"
        r"(?:f\{(?P<formula>[A-Za-z0-9\[\]]+)\})|"
        r"(?:_\{(?P<named_compound>[^\{\}\s,/]+)\})|"
        r"(?:s\{(?P<smiles>[^\}]+)\})|"
        r"(?:(?P<unannotated>\?)(?P<unannotated_label>\d+)?)"
        r")"
        r"(?P<neutral_losses>(?:[+-](?:\d+(?:\.\d+)?|\d*(?:(?:(?:\[[0-9]+[A-Z][A-Za-z0-9]*\])|(?:[A-Z][A-Za-z0-9]*))+)|(?:\[(?:(?:[A-Za-z0-9:\.]+)(?:\[(?:[A-Za-z0-9\.:\-]+)\])?)\])))+)?"
        r"(?P<isotope>(?:(?:[+-]\d*)i(?:(?:\d+)?(?:[A-Z][a-z]*)?|A)?)+)?"
        r"(?:\[(?P<adducts>M(?:[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?"
        r"(?:\^(?P<charge>[+-]?\d+))?"
        r"(?:/(?P<mass_error>-?\d+(?:\.\d+)?)(?P<mass_error_unit>ppm)?)?"
        r"(?:\*(?P<confidence>\d*(?:\.\d+)?))?"
        r"$"
    )

    def parse_single(self, annotation_str: str) -> FragmentAnnotation:
        """Parse a single annotation string"""
        match = self.PATTERN.match(annotation_str)
        if not match:
            raise ValueError(f"Invalid mzPAF annotation: {annotation_str}")

        groups = match.groupdict()

        # Parse ion type
        ion_type = self._parse_ion_type(groups)

        # Parse isotope (returns list or None)
        isotope = self._parse_isotope(groups.get("isotope"))

        # Parse neutral losses (returns list or None)
        neutral_losses = self._parse_neutral_losses(groups.get("neutral_losses"))

        # Parse mass error
        mass_error = None
        if groups.get("mass_error"):
            unit = groups.get("mass_error_unit") or "da"
            if unit not in ("da", "ppm"):
                raise ValueError(f"Unknown mass error unit: {unit}")
            mass_error = MassError(float(groups["mass_error"]), unit)

        # Parse charge
        charge = int(groups["charge"]) if groups.get("charge") else 1

        # Parse confidence
        confidence = float(groups["confidence"]) if groups.get("confidence") else None

        # Parse analyte reference
        analyte_ref = (
            int(groups["analyte_reference"])
            if groups.get("analyte_reference")
            else None
        )

        # Parse adducts (returns list or None)
        adducts = self._parse_adducts(groups.get("adducts"))

        return FragmentAnnotation(
            ion_type=ion_type,
            analyte_reference=analyte_ref,
            is_auxiliary=bool(groups.get("is_auxiliary")),
            neutral_losses=neutral_losses,
            isotope=isotope,
            adducts=adducts,
            charge=charge,
            mass_error=mass_error,
            confidence=confidence,
        )

    def _parse_ion_type(
        self, groups: dict[str, str | None]
    ) -> (
        PeptideIon
        | InternalFragment
        | ImmoniumIon
        | ReferenceIon
        | NamedCompound
        | ChemicalFormula
        | SMILESCompound
        | UnknownIon
        | PrecursorIon
    ):
        """Parse the ion type from regex groups."""

        # Peptide ion (b, y, c, z series)
        if groups.get("series"):
            ordinal = groups["ordinal"]
            if ordinal is None:
                raise ValueError("Peptide ion missing ordinal position")
            if "series" not in groups or not isinstance(groups["series"], str):
                raise ValueError("Peptide ion missing series type")
            return PeptideIon(
                series=groups["series"],
                position=int(ordinal),
                sequence=groups.get("sequence_ordinal"),
            )

        # Internal fragment
        if groups.get("internal_start"):
            start = groups["internal_start"]
            end = groups["internal_end"]
            if start is None or end is None:
                raise ValueError("Internal fragment missing start or end position")

            return InternalFragment(
                start_position=int(start),
                end_position=int(end),
                sequence=groups.get("sequence_internal"),
            )

        # Precursor ion
        if groups.get("precursor"):
            return PrecursorIon()

        # Immonium ion
        if groups.get("immonium"):
            if "immonium" not in groups or not isinstance(groups["immonium"], str):
                raise ValueError("Immonium ion missing amino acid")
            return ImmoniumIon(
                amino_acid=groups["immonium"],
                modification=groups.get("immonium_modification"),
            )

        # Reference ion
        if groups.get("reference"):
            if "reference" not in groups or not isinstance(groups["reference"], str):
                raise ValueError("Reference ion missing label")
            return ReferenceIon(name=groups["reference"])

        # Chemical formula
        if groups.get("formula"):
            if "formula" not in groups or not isinstance(groups["formula"], str):
                raise ValueError("Chemical formula missing")
            return ChemicalFormula(formula=groups["formula"])

        # Named compound
        if groups.get("named_compound"):
            if "named_compound" not in groups or not isinstance(
                groups["named_compound"], str
            ):
                raise ValueError("Named compound missing")
            return NamedCompound(name=groups["named_compound"])

        # SMILES compound
        if groups.get("smiles"):
            if "smiles" not in groups or not isinstance(groups["smiles"], str):
                raise ValueError("SMILES compound missing")
            return SMILESCompound(smiles=groups["smiles"])

        # Unknown/unannotated ion
        if groups.get("unannotated"):
            label_str = groups.get("unannotated_label")
            label = int(label_str) if label_str else None
            return UnknownIon(label=label)

        # Nothing matched
        available_keys = [k for k, v in groups.items() if v is not None]
        raise ValueError(
            f"Unable to parse ion type from annotation. "
            f"Non-null keys found: {available_keys}"
        )

    def _parse_isotope(
        self, isotope_str: str | None
    ) -> list[IsotopeSpecification] | None:
        """Parse isotope notation - returns list or None"""
        if not isotope_str:
            return None

        isotopes: list[IsotopeSpecification] = []
        # Pattern to match individual isotope specifications
        pattern = r"([+-]?)(\d*)i((?:\d+)?(?:[A-Z][a-z]*)?|A)?"

        for match in re.finditer(pattern, isotope_str):
            sign_str = match.group(1)
            count_str = match.group(2)
            element_or_avg = match.group(3)

            sign = -1 if sign_str == "-" else 1
            count = int(count_str) if count_str else 1
            count *= sign

            # Check if it's an averaged isotopomer
            if element_or_avg == "A":
                isotopes.append(IsotopeSpecification(count=count, is_average=True))
            else:
                element = element_or_avg if element_or_avg else None
                isotopes.append(IsotopeSpecification(count=count, element=element))

        return isotopes if isotopes else None

    def _parse_neutral_losses(self, losses_str: str | None) -> list[str] | None:
        """Parse neutral losses/gains (mass, formula, or named) - returns list or None"""
        if not losses_str:
            return None

        # Match: +/- followed by either:
        # - A decimal number (mass): +17.03, -18.01
        # - A count + formula: -H2O, -2H2O, +NH3
        # - A bracketed modification: -[Phospho]
        pattern = r"([+-](?:\d+(?:\.\d+)?|\d*(?:(?:\[[^\]]+\])|(?:[A-Z][A-Za-z0-9]*)+|\[(?:[^\]]+)\])))"
        losses = re.findall(pattern, losses_str)
        return losses if losses else None

    def _parse_adducts(self, adduct_str: str | None) -> list[str] | None:
        """Parse adduct notation - returns list or None"""
        if not adduct_str:
            return None

        # For now, return as a single-item list
        # Future enhancement could parse multiple adducts if the spec supports it
        return [adduct_str]

    def parse(self, annotation_str: str) -> list[FragmentAnnotation]:
        """Parse potentially multiple comma-separated annotations"""
        if not annotation_str:
            return []

        # Split by comma to handle multiple annotations
        parts = annotation_str.split(",")
        annotations: list[FragmentAnnotation] = []

        for part in parts:
            part = part.strip()
            if part:
                annotations.append(self.parse_single(part))

        return annotations


# Example usage and testing
if __name__ == "__main__":
    parser = mzPAFParser()

    # Test cases including enhanced internal fragment notation
    test_cases = [
        # Standard test cases
        "b2-H2O/3.2ppm",
        "b2-18/3.2ppm",
        "b4-H2O^2/3.2ppm",
        "1@y12/0.13,2@b9-NH3/0.23",
        "y4-H2O+H2O+2i+3i[M+H+Na]^2",
        "IY[Phospho]",
        "r[TMT127N]",
        "?^3",
        # Mass error test cases (from spec)
        "y1/-1.4ppm",
        "y1/-0.0002",
        "y1/0.0002",  # positive, no sign
        # Internal fragment test cases
        "m3:6",  # Default yb cleavage
        "m3:6-CO-H20",  # yb with CO loss (equivalent to ya)
        # Neutral loss test cases with mass values
        "b5-18.01",  # Mass-based loss (water)
        "y3+17.03",  # Mass-based gain (ammonia)
        "b7-97.98",  # Mass-based loss (phosphoric acid)
        "y5-18.01-17.03",  # Multiple mass losses
        "b4-H2O-18.01",  # Mixed formula and mass
    ]

    for test in test_cases:
        print(f"\nOriginal: {test}")
        try:
            annotations = parser.parse(test)
            for ann in annotations:
                reconstructed = ann.serialize()
                print(f"  Parsed and reconstructed: {reconstructed}")

                # If it's an internal fragment, show the cleavage types
                if isinstance(ann.ion_type, InternalFragment):
                    frag = ann.ion_type
                    print(frag.serialize())
                print(ann.neutral_losses)

        except Exception as e:
            print(f"  Error: {e}")

    print("\n✅ Testing complete!")
