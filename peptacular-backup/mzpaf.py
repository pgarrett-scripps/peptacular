"""
mzPAF Parser Module - Enhanced with backbone cleavage notation for internal fragments
Parses and serializes mzPAF (Peak Annotation Format) strings for peptide mass spectra.
Based on mzPAF specification v1.0 with extended internal fragment support
"""

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Union


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
    unit: str = "Da"  # "Da" or "ppm"

    def __str__(self):
        if self.unit == "ppm":
            return f"{self.value:g}ppm"
        return f"{self.value:g}"


@dataclass
class IsotopeSpecification:
    """Represents isotope information"""

    count: int = 0  # number of isotopes above/below monoisotope
    element: Optional[str] = None  # e.g., "13C", "15N"
    is_average: bool = False  # True for averaged isotopomers

    def __str__(self):
        if self.count == 0:
            return ""

        sign = "+" if self.count > 0 else ""
        count_str = "" if abs(self.count) == 1 else str(abs(self.count))

        if self.is_average:
            return f"{sign}{count_str}iA"
        elif self.element:
            return f"{sign}{count_str}i{self.element}"
        else:
            return f"{sign}{count_str}i"


@dataclass
class PeptideIon:
    """Represents a primary peptide fragment ion"""

    series: str  # b, y, a, x, c, z, etc.
    position: int
    sequence: Optional[str] = None  # ProForma sequence

    def __str__(self):
        result = f"{self.series}{self.position}"
        if self.sequence:
            result += f"{{{self.sequence}}}"
        return result


@dataclass
class InternalFragment:
    """Represents an internal fragment ion with optional backbone cleavage specification"""

    start_position: int
    end_position: int
    n_terminal_cleavage: str = "y"  # Default to y-type cleavage at N-terminal side
    c_terminal_cleavage: str = "b"  # Default to b-type cleavage at C-terminal side
    sequence: Optional[str] = None

    def __str__(self):
        # If using default yb cleavage, just use 'm'
        if self.n_terminal_cleavage == "y" and self.c_terminal_cleavage == "b":
            result = f"m{self.start_position}:{self.end_position}"
        else:
            # Use extended notation with cleavage types
            result = f"m{self.n_terminal_cleavage}{self.c_terminal_cleavage}{self.start_position}:{self.end_position}"

        if self.sequence:
            result += f"{{{self.sequence}}}"
        return result


@dataclass
class ImmoniumIon:
    """Represents an immonium ion"""

    amino_acid: str
    modification: Optional[str] = None

    def __str__(self):
        result = f"I{self.amino_acid}"
        if self.modification:
            result += f"[{self.modification}]"
        return result


@dataclass
class ReferenceIon:
    """Represents a reference ion"""

    name: str

    def __str__(self):
        return f"r[{self.name}]"


@dataclass
class NamedCompound:
    """Represents a named compound"""

    name: str

    def __str__(self):
        return f"_{{{self.name}}}"


@dataclass
class ChemicalFormula:
    """Represents a chemical formula"""

    formula: str

    def __str__(self):
        return f"f{{{self.formula}}}"


@dataclass
class SMILESCompound:
    """Represents a SMILES string"""

    smiles: str

    def __str__(self):
        return f"s{{{self.smiles}}}"


@dataclass
class UnknownIon:
    """Represents an unknown/unannotated ion"""

    label: Optional[int] = None

    def __str__(self):
        if self.label is not None:
            return f"?{self.label}"
        return "?"


@dataclass
class PrecursorIon:
    """Represents a precursor ion"""

    def __str__(self):
        return "p"


@dataclass
class FragmentAnnotation:
    """Complete fragment ion annotation"""

    # Core ion description
    ion_type: Union[
        PeptideIon,
        InternalFragment,
        ImmoniumIon,
        ReferenceIon,
        NamedCompound,
        ChemicalFormula,
        SMILESCompound,
        UnknownIon,
        PrecursorIon,
    ]

    # Optional components
    analyte_reference: Optional[int] = None
    is_auxiliary: bool = False
    neutral_losses: List[str] = field(default_factory=list)
    isotope: Optional[IsotopeSpecification] = None
    adducts: Optional[str] = None  # e.g., "M+H+Na"
    charge: int = 1
    mass_error: Optional[MassError] = None
    confidence: Optional[float] = None

    def to_string(self) -> str:
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
        for loss in self.neutral_losses:
            parts.append(loss)

        # Isotope
        if self.isotope and self.isotope.count != 0:
            parts.append(str(self.isotope))

        # Adduct type
        if self.adducts:
            parts.append(f"[{self.adducts}]")

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


class mzPAFParser:
    """Parser for mzPAF annotation strings with enhanced internal fragment support"""

    # Updated regex pattern to support backbone cleavage notation for internal fragments
    PATTERN = re.compile(
        r"^(?P<is_auxiliary>&)?"
        r"(?:(?P<analyte_reference>\d+)@)?"
        r"(?:"
        r"(?:(?P<series>(?:da|db|wa|wb)|[axbyczdwv]\.?)(?P<ordinal>\d+)(?:\{(?P<sequence_ordinal>.+)\})?)|"
        # Enhanced internal fragment pattern - fixed parenthesis issue
        r"(?:m(?P<backbone_cleavage>[xyzabc]{2})?(?P<internal_start>\d+):(?P<internal_end>\d+)(?:\{(?P<sequence_internal>.+)\})?)|"
        r"(?P<precursor>p)|"
        r"(?:I(?P<immonium>[A-Z])(?:\[(?P<immonium_modification>(?:[^\]]+))\])?)|"
        r"(?P<reference>r(?:(?:\[(?P<reference_label>[^\]]+)\])))|"
        r"(?:f\{(?P<formula>[A-Za-z0-9\[\]]+)\})|"
        r"(?:_\{(?P<named_compound>[^\{\}\s,/]+)\})|"
        r"(?:s\{(?P<smiles>[^\}]+)\})|"
        r"(?:(?P<unannotated>\?)(?P<unannotated_label>\d+)?)"
        r")"
        r"(?P<neutral_losses>(?:[+-]\d*(?:(?:(?:(?:\[[0-9]+[A-Z][A-Za-z0-9]*\])|(?:[A-Z][A-Za-z0-9]*))+)|(?:\[(?:(?:[A-Za-z0-9:\.]+)(?:\[(?:[A-Za-z0-9\.:\-]+)\])?)\])))+)?"
        r"(?P<isotope>(?:[+-]\d*)i(?:(?:\d+)?(?:[A-Z][a-z]*)?|A)?)?"
        r"(?:\[(?P<adducts>M(?:[+-]\d*[A-Z][A-Za-z0-9]*)+)\])?"
        r"(?:\^(?P<charge>[+-]?\d+))?"
        r"(?:/(?P<mass_error>[+-]?\d+(?:\.\d+)?)(?P<mass_error_unit>ppm)?)?"
        r"(?:\*(?P<confidence>\d*(?:\.\d+)?))?"
    )

    def parse_single(self, annotation_str: str) -> FragmentAnnotation:
        """Parse a single annotation string"""
        match = self.PATTERN.match(annotation_str)
        if not match:
            raise ValueError(f"Invalid mzPAF annotation: {annotation_str}")

        groups = match.groupdict()

        # Parse ion type
        ion_type = self._parse_ion_type(groups)

        # Parse isotope
        isotope = self._parse_isotope(groups.get("isotope"))

        # Parse neutral losses
        neutral_losses = self._parse_neutral_losses(groups.get("neutral_losses"))

        # Parse mass error
        mass_error = None
        if groups.get("mass_error"):
            mass_error = MassError(
                float(groups["mass_error"]), groups.get("mass_error_unit", "Da")
            )

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

        return FragmentAnnotation(
            ion_type=ion_type,
            analyte_reference=analyte_ref,
            is_auxiliary=bool(groups.get("is_auxiliary")),
            neutral_losses=neutral_losses,
            isotope=isotope,
            adducts=groups.get("adducts"),
            charge=charge,
            mass_error=mass_error,
            confidence=confidence,
        )

    def _parse_ion_type(self, groups: dict):
        """Parse the ion type from regex groups"""
        if groups.get("series"):
            return PeptideIon(
                series=groups["series"],
                position=int(groups["ordinal"]),
                sequence=groups.get("sequence_ordinal"),
            )
        elif (
            groups.get("internal_start") is not None
        ):  # Changed: check if not None instead of truthiness
            # Parse internal fragment with optional backbone cleavage specification
            backbone = groups.get("backbone_cleavage", "")

            if len(backbone) == 2:
                # Extended notation with specified cleavage types
                n_terminal_cleavage = backbone[0]
                c_terminal_cleavage = backbone[1]
            else:
                # Default yb cleavage
                n_terminal_cleavage = "y"
                c_terminal_cleavage = "b"

            return InternalFragment(
                start_position=int(groups["internal_start"]),
                end_position=int(groups["internal_end"]),
                n_terminal_cleavage=n_terminal_cleavage,
                c_terminal_cleavage=c_terminal_cleavage,
                sequence=groups.get("sequence_internal"),
            )
        elif groups.get("precursor"):
            return PrecursorIon()
        elif groups.get("immonium"):
            return ImmoniumIon(
                amino_acid=groups["immonium"],
                modification=groups.get("immonium_modification"),
            )
        elif groups.get("reference"):
            return ReferenceIon(name=groups["reference_label"])
        elif groups.get("formula"):
            return ChemicalFormula(formula=groups["formula"])
        elif groups.get("named_compound"):
            return NamedCompound(name=groups["named_compound"])
        elif groups.get("smiles"):
            return SMILESCompound(smiles=groups["smiles"])
        elif groups.get("unannotated"):
            label = (
                int(groups["unannotated_label"])
                if groups.get("unannotated_label")
                else None
            )
            return UnknownIon(label=label)
        else:
            raise ValueError("Could not determine ion type")

    def _parse_isotope(
        self, isotope_str: Optional[str]
    ) -> Optional[IsotopeSpecification]:
        """Parse isotope notation"""
        if not isotope_str:
            return None

        # Check for averaged isotopomer
        if isotope_str.endswith("A"):
            count = 1
            if len(isotope_str) > 2:
                # Extract count
                count_match = re.match(r"([+-]?)(\d*)iA", isotope_str)
                if count_match:
                    sign = -1 if count_match.group(1) == "-" else 1
                    count = int(count_match.group(2)) if count_match.group(2) else 1
                    count *= sign
            return IsotopeSpecification(count=count, is_average=True)

        # Parse regular isotope
        pattern = r"([+-]?)(\d*)i(\d*[A-Z][a-z]*)?"
        match = re.match(pattern, isotope_str)
        if match:
            sign = -1 if match.group(1) == "-" else 1
            count = int(match.group(2)) if match.group(2) else 1
            count *= sign
            element = match.group(3) if match.group(3) else None
            return IsotopeSpecification(count=count, element=element)

        return None

    def _parse_neutral_losses(self, losses_str: Optional[str]) -> List[str]:
        """Parse neutral losses/gains"""
        if not losses_str:
            return []

        # Split on + or - (keeping the sign)
        pattern = r"([+-]\d*(?:(?:\[[^\]]+\])|(?:[A-Z][A-Za-z0-9]*)+|\[(?:[^\]]+)\]))"
        losses = re.findall(pattern, losses_str)
        return losses

    def parse(self, annotation_str: str) -> List[FragmentAnnotation]:
        """Parse potentially multiple comma-separated annotations"""
        if not annotation_str:
            return []

        # Split by comma to handle multiple annotations
        parts = annotation_str.split(",")
        annotations = []

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
        "b4-H2O^2/3.2ppm",
        "1@y12/0.13,2@b9-NH3/0.23",
        "y4-H2O+2i[M+H+Na]^2",
        "IY[Phospho]",
        "r[TMT127N]",
        "p-H3PO4^2",
        "?^3",
        # Internal fragment test cases
        "m3:6",  # Default yb cleavage
        "m3:6-CO",  # yb with CO loss (equivalent to ya)
        "mxa3:6",  # xa cleavage specified
        "mzc3:6",  # zc cleavage specified
        "myc3:6-H2O^2",  # yc cleavage with water loss, doubly charged
        "mab5:8{PEPT}",  # ab cleavage with sequence
        "1@mxb2:4/1.2ppm*0.85",  # xb cleavage with analyte ref, mass error, and confidence
    ]

    for test in test_cases:
        print(f"\nOriginal: {test}")
        try:
            annotations = parser.parse(test)
            for ann in annotations:
                reconstructed = ann.to_string()
                print(f"  Parsed and reconstructed: {reconstructed}")

                # If it's an internal fragment, show the cleavage types
                if isinstance(ann.ion_type, InternalFragment):
                    frag = ann.ion_type
                    print(f"    N-terminal cleavage: {frag.n_terminal_cleavage}")
                    print(f"    C-terminal cleavage: {frag.c_terminal_cleavage}")

        except Exception as e:
            print(f"  Error: {e}")

    print("\n✅ Testing complete!")
