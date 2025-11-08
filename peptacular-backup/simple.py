"""
Modern implementation of the mzPAF product ion annotation format.

This implementation includes parsing the compact text format and serializing,
with full type annotations and simplified architecture.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Dict, List, Optional, Protocol, Union

# Import peptacular utilities for mass/composition calculations
from ..chem import mod_comp
from ..chem.chem_util import parse_chem_formula
from ..constants import ISOTOPIC_ATOMIC_MASSES, NEUTRON_MASS, PROTON_MASS
from ..mass_calc import mod_mass
from ..proforma import parse
from ..proforma.mass_comp import chem_mass
from .reference_molecule import ReferenceMolecule

# Type aliases
JSONDict = Dict[str, Any]
Composition = Dict[str, int]

# Core regex pattern for mzPAF - preserved for standard compatibility
ANNOTATION_PATTERN = re.compile(
    r"""
^(?P<is_auxiliary>&)?
   (?:(?P<analyte_reference>\d+)@)?
   (?:(?:(?P<series>(?:da|db|wa|wb)|[axbyczdwv]\.?)(?P<ordinal>\d+)(?:\{(?P<sequence_ordinal>.+)\})?)|
   (?P<series_internal>m(?P<backbone_cleavage>[a-z]{0,2})(?P<internal_start>\d+):(?P<internal_end>\d+)(?:\{(?P<sequence_internal>.+)\})?)|
   (?P<precursor>p(?:\{(?P<sequence_precursor>.+)\})?)|
   (:?I(?P<immonium>[A-Z])(?:\[(?P<immonium_modification>(?:[^\]]+))\])?)|
   (?P<reference>r(?:
    (?:\[
        (?P<reference_label>[^\]]+)
    \])
   ))|
   (?:f\{(?P<formula>[A-Za-z0-9\[\]]+)\})|
   (?:_\{
    (?P<named_compound>[^\{\}\s,/]+)
    \})|
   (?:s\{(?P<smiles>[^\}]+)\})|
   (?:(?P<unannotated>\?)(?P<unannotated_label>\d+)?)
)
(?P<neutral_losses>(?:[+-]\d*
    (?:(?:(?:(?:\[[0-9]+[A-Z][A-Za-z0-9]*\])|(?:[A-Z][A-Za-z0-9]*))+)|
        (?:\[
            (?:
                (?:[A-Za-z0-9:\.]+)(?:\[(?:[A-Za-z0-9\.:\-\ ]+)\])?
            )
            \])
    )
)+)?
(?P<isotope>
    (?:[+-]\d*)i(?:(?:\d+)(?:[A-Z][a-z]?)|(?:A))?
)?
(?:\[(?P<adducts>M(:?[+-]\d*(?:\[\d+\])?[A-Z][A-Za-z0-9]*)+)\])?
(?:\^(?P<charge>[+-]?\d+))?
(?:/(?P<mass_error>[+-]?\d+(?:\.\d+)?)(?P<mass_error_unit>ppm)?)?
(?:\*(?P<confidence>\d*(?:\.\d+)?))?
""",
    re.X,
)

ISOTOPE_PATTERN = re.compile(
    r"(?P<isotope>[+-]\d*)i(?:(?P<nucleon_count>\d+)(?P<element>[A-Z][a-z]?)|(?P<average_isotopologue>A))?"
)

NEUTRAL_LOSS_PATTERN = re.compile(
    r"""\s*(?P<neutral_loss>(?:(?P<sign>[+-])?\s*(?P<coefficient>\d*)\s*
    (?:(?P<formula>(?:(?:[A-Z][A-Za-z0-9]*)|(?:\[\d+[A-Z][A-Za-z0-9]\]*))+)|
        (?P<braced_name>\[
            (?:
                (?:[A-Za-z0-9:\.]+)(?:\[(?:[A-Za-z0-9\.:\-\ ]+)\])?
            )
            \])
    )
))""",
    re.X,
)


def count_elements_in_smiles(smiles_string, include_hydrogen=True):
    """
    Parse a SMILES string and return a counter of elements.

    Parameters
    ----------
    smiles_string : str
        The SMILES string to parse
    include_hydrogen : bool
        Whether to include implicit hydrogens in the count

    Returns
    -------
    Counter
        A Counter object with element symbols as keys and counts as values

    Examples
    --------
    >>> count_elements_in_smiles('CCO')
    Counter({'C': 2, 'O': 1, 'H': 6})
    >>> count_elements_in_smiles('CCO', include_hydrogen=False)
    Counter({'C': 2, 'O': 1})
    >>> count_elements_in_smiles('c1ccccc1')  # benzene
    Counter({'C': 6, 'H': 6})
    """
    from pysmiles import read_smiles  # Import here to avoid hard dependency if not used

    # Parse the SMILES string
    mol = read_smiles(smiles_string, explicit_hydrogen=False)

    element_counter = Counter()

    # Count each atom in the molecule
    for node in mol.nodes():
        node_data = mol.nodes[node]

        # Get the element, defaulting to '*' for wildcards
        element = node_data.get("element", "*")
        element_counter[element] += 1

        # Add implicit hydrogens if requested
        if include_hydrogen and "hcount" in node_data:
            h_count = node_data["hcount"]
            if h_count > 0:
                element_counter["H"] += h_count

    return element_counter


class NeutralLossType(Enum):
    """Classification for neutral loss types."""

    FORMULA = auto()
    REFERENCE = auto()
    BRACED_NAME = auto()


class IonSeriesType(Enum):
    """Types of ion series."""

    PEPTIDE = auto()
    INTERNAL = auto()
    PRECURSOR = auto()
    IMMONIUM = auto()
    REFERENCE = auto()
    FORMULA = auto()
    SMILES = auto()
    NAMED_COMPOUND = auto()
    UNANNOTATED = auto()


@dataclass(frozen=True)
class MassError:
    """Represents mass error with units."""

    value: float
    unit: str = "Da"

    def __str__(self) -> str:
        unit = "" if self.unit == "Da" else self.unit
        return f"{self.value}{unit}"

    def to_json(self) -> JSONDict:
        return {"value": self.value, "unit": self.unit}

    @classmethod
    def from_json(cls, data: JSONDict) -> MassError:
        return cls(data["value"], data.get("unit", "Da"))


@dataclass(frozen=True)
class IsotopeLabel:
    """Represents an isotope specification."""

    offset: int
    element: Optional[str] = None
    nucleon_count: Optional[int] = None
    is_average: bool = False

    def __str__(self) -> str:
        suffix = ""
        if self.element and self.nucleon_count:
            suffix = f"{self.nucleon_count}{self.element}"
        elif self.is_average:
            suffix = "A"

        prefix = ""
        if self.offset >= 0:
            prefix = f"+{self.offset}" if self.offset > 1 else "+"
        else:
            prefix = f"{self.offset}" if self.offset < -1 else "-"

        return f"{prefix}i{suffix}"

    def to_json(self) -> JSONDict:
        data = {"isotope": self.offset}
        if self.is_average:
            data["averaged"] = True
        elif self.element:
            data["element"] = self.element
            data["nucleon_count"] = self.nucleon_count
        return data

    @classmethod
    def from_string(cls, text: str) -> Optional[IsotopeLabel]:
        match = ISOTOPE_PATTERN.search(text)
        if not match:
            return None

        groups = match.groupdict()
        offset = cls._parse_sign_or_int(groups.get("isotope", ""))

        return cls(
            offset=offset,
            element=groups.get("element"),
            nucleon_count=int(groups["nucleon_count"])
            if groups.get("nucleon_count")
            else None,
            is_average=bool(groups.get("average_isotopologue")),
        )

    @staticmethod
    def _parse_sign_or_int(text: str) -> int:
        if text == "+":
            return 1
        elif text == "-":
            return -1
        return int(text) if text else 0


@dataclass(frozen=True)
class NeutralLoss:
    """Represents a neutral loss or gain."""

    name: str
    coefficient: int = 1
    loss_type: NeutralLossType = NeutralLossType.FORMULA

    def __str__(self) -> str:
        name = (
            f"[{self.name}]" if self.loss_type != NeutralLossType.FORMULA else self.name
        )

        if self.coefficient >= 0:
            prefix = f"+{self.coefficient}" if self.coefficient > 1 else "+"
        else:
            prefix = f"{self.coefficient}" if self.coefficient < -1 else "-"

        return f"{prefix}{name}"

    def to_json(self) -> str:
        return str(self)

    @classmethod
    def parse_many(cls, text: str) -> List[NeutralLoss]:
        if not text:
            return []

        losses = []
        for match in NEUTRAL_LOSS_PATTERN.finditer(text):
            groups = match.groupdict()
            coef = int(groups.get("coefficient") or 1)
            if groups.get("sign") == "-":
                coef = -coef

            name = groups.get("formula") or groups.get("braced_name", "")
            loss_type = (
                NeutralLossType.FORMULA
                if groups.get("formula")
                else NeutralLossType.BRACED_NAME
            )

            if name.startswith("[") and name.endswith("]"):
                name = name[1:-1]
                loss_type = NeutralLossType.BRACED_NAME

            losses.append(cls(name, coef, loss_type))

        return losses


# Base classes for ion types
@dataclass
class IonType(Protocol):
    """Protocol for ion types."""

    def format(self) -> str: ...
    def to_json(self) -> JSONDict: ...
    def mass(self) -> float: ...
    def comp(self) -> dict: ...


@dataclass(frozen=True)
class PeptideIon:
    """Standard peptide fragment ion (a, b, c, x, y, z)."""

    series: str
    position: int
    sequence: Optional[str] = None

    def format(self) -> str:
        base = f"{self.series}{self.position}"
        return f"{base}{{{self.sequence}}}" if self.sequence else base

    def to_json(self) -> JSONDict:
        data = {"series": self.series, "position": self.position}
        if self.sequence:
            data["sequence"] = self.sequence
        return data

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass of the peptide fragment."""
        if not self.sequence:
            raise ValueError(
                "Cannot calculate mass without sequence. See the referenced analyte."
            )
        ann = parse(self.sequence)
        # Get fragment masses for this series
        return ann.mass(ion_type=self.series, monoisotopic=monoisotopic)

    def comp(self) -> dict:
        """Get chemical composition of the peptide fragment."""
        if not self.sequence:
            raise ValueError(
                "Cannot calculate composition without sequence. See the referenced analyte."
            )
        ann = parse(self.sequence)
        # Get fragment compositions for this series
        return ann.comp(ion_type=self.series)


@dataclass(frozen=True)
class InternalFragment:
    """Internal peptide fragment."""

    start: int
    end: int
    n_terminal_cleavage: str = "y"
    c_terminal_cleavage: str = "b"
    sequence: Optional[str] = None

    def format(self) -> str:
        # Only include backbone specification if not default yb
        backbone = ""
        if not (self.n_terminal_cleavage == "y" and self.c_terminal_cleavage == "b"):
            backbone = f"{self.n_terminal_cleavage}{self.c_terminal_cleavage}"

        base = f"m{backbone}{self.start}:{self.end}"
        return f"{base}{{{self.sequence}}}" if self.sequence else base

    def to_json(self) -> JSONDict:
        data = {
            "start_position": self.start,
            "end_position": self.end,
            "n_terminal_cleavage": self.n_terminal_cleavage,
            "c_terminal_cleavage": self.c_terminal_cleavage,
        }
        if self.sequence:
            data["sequence"] = self.sequence
        return data

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass of the internal fragment."""
        if not self.sequence:
            raise ValueError(
                "Cannot calculate mass without sequence. See the referenced analyte."
            )
        # Internal fragments require extracting the subsequence
        ann = parse(self.sequence)
        return ann.mass(
            ion_type=f"{self.c_terminal_cleavage}{self.n_terminal_cleavage}",
            monoisotopic=monoisotopic,
        )

    def comp(self) -> dict:
        if not self.sequence:
            raise ValueError(
                "Cannot calculate composition without sequence. See the referenced analyte."
            )
        # Internal fragments require extracting the subsequence
        ann = parse(self.sequence)
        # Extract residues from start to end (1-indexed positions)
        return ann.comp(
            ion_type=f"{self.c_terminal_cleavage}{self.n_terminal_cleavage}"
        )


@dataclass(frozen=True)
class PrecursorIon:
    """Precursor ion."""

    sequence: Optional[str] = None

    def format(self) -> str:
        if self.sequence:
            return f"p{{{self.sequence}}}"
        return "p"

    def to_json(self) -> JSONDict:
        data = {}
        if self.sequence:
            data["sequence"] = self.sequence
        return data

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate precursor mass."""
        if not self.sequence:
            raise ValueError(
                "Cannot calculate precursor mass without sequence or analyte reference."
            )
        ann = parse(self.sequence)
        return ann.mass(monoisotopic=monoisotopic)

    def comp(self) -> dict:
        """Get chemical composition of the precursor."""
        if not self.sequence:
            raise ValueError(
                "Cannot calculate composition without sequence or analyte reference."
            )
        ann = parse(self.sequence)
        return ann.comp()


@dataclass(frozen=True)
class ImmoniumIon:
    """Immonium ion from a single amino acid."""

    amino_acid: str
    modification: Optional[str] = None

    def format(self) -> str:
        mod = f"[{self.modification}]" if self.modification else ""
        return f"I{self.amino_acid}{mod}"

    def to_json(self) -> JSONDict:
        data = {"amino_acid": self.amino_acid}
        if self.modification:
            data["modification"] = self.modification
        return data

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate immonium ion mass."""
        # Base immonium formula: amino acid - CO - NH

        annot = parse(self.amino_acid)
        return annot.mass(ion_type="I", monoisotopic=monoisotopic)

    def comp(self) -> dict:
        """Get chemical composition of the immonium ion."""
        annot = parse(self.amino_acid)
        return annot.comp(ion_type="I")


@dataclass(frozen=True)
class ReferenceIon:
    """Reference ion (e.g., reporter ions)."""

    name: str

    def format(self) -> str:
        return f"r[{self.name}]"

    def to_json(self) -> JSONDict:
        return {"reference": self.name}

    def mass(self, monoisotopic: bool = True) -> float:
        """Get mass from reference molecule database."""

        try:
            ref = ReferenceMolecule.get(self.name)
            return ref.mass(monoisotopic=monoisotopic)
        except KeyError:
            raise ValueError(f"Unknown reference molecule: {self.name}")

    def comp(self) -> dict:
        """Get composition from reference molecule database."""
        try:
            ref = ReferenceMolecule.get(self.name)
            return ref.comp()
        except KeyError:
            raise ValueError(f"Unknown reference molecule: {self.name}")


@dataclass(frozen=True)
class FormulaIon:
    """Ion specified by chemical formula."""

    formula: str

    def format(self) -> str:
        return f"f{{{self.formula}}}"

    def to_json(self) -> JSONDict:
        return {"formula": self.formula}

    def mass(self, monoisotopic: bool = True) -> float:
        """Calculate mass from chemical formula."""
        return chem_mass(self.formula)

    def composition(self) -> Composition:
        """Get chemical composition."""
        return parse_chem_formula(self.formula)


@dataclass(frozen=True)
class SMILESIon:
    """Ion specified by SMILES notation."""

    smiles: str

    def format(self) -> str:
        return f"s{{{self.smiles}}}"

    def to_json(self) -> JSONDict:
        return {"smiles": self.smiles}

    def mass(self, monoisotopic: bool = True) -> float:
        """SMILES mass calculation requires external library."""
        comp = count_elements_in_smiles(self.smiles, include_hydrogen=True)
        return chem_mass(comp)

    def comp(self) -> dict:
        """SMILES composition calculation requires external library."""
        comp = count_elements_in_smiles(self.smiles, include_hydrogen=True)
        return dict(comp)


@dataclass(frozen=True)
class NamedCompound:
    """Named compound ion."""

    name: str

    def format(self) -> str:
        return f"_{{{self.name}}}"

    def to_json(self) -> JSONDict:
        return {"compound_name": self.name}

    def mass(self, monoisotopic: bool = True) -> float:
        """Named compound mass requires lookup."""
        raise ValueError(f"Cannot determine mass for named compound: {self.name}")

    def comp(self) -> dict:
        """Named compound composition requires lookup."""
        raise ValueError(
            f"Cannot determine composition for named compound: {self.name}"
        )


@dataclass(frozen=True)
class UnknownIon:
    """Unannotated/unknown ion."""

    label: Optional[int] = None

    def format(self) -> str:
        return f"?{self.label}" if self.label is not None else "?"

    def to_json(self) -> JSONDict:
        return {"unannotated_label": self.label}

    def mass(self, monoisotopic: bool = True) -> float:
        """Unknown ions have no calculable mass."""
        raise ValueError("Cannot calculate mass for unknown ion")

    def comp(self) -> dict:
        """Unknown ions have no calculable composition."""
        raise ValueError("Cannot calculate composition for unknown ion")


# Main annotation class
@dataclass
class FragmentAnnotation:
    """Complete fragment ion annotation."""

    ion_type: Union[
        PeptideIon,
        InternalFragment,
        PrecursorIon,
        ImmoniumIon,
        ReferenceIon,
        FormulaIon,
        SMILESIon,
        NamedCompound,
        UnknownIon,
    ]
    neutral_losses: List[NeutralLoss] = field(default_factory=list)
    isotopes: List[IsotopeLabel] = field(default_factory=list)
    adducts: List[str] = field(default_factory=list)
    charge: int = 1
    analyte_reference: Optional[int] = None
    mass_error: Optional[MassError] = None
    confidence: Optional[float] = None
    is_auxiliary: bool = False

    def __str__(self) -> str:
        parts = []

        # Analyte reference
        if self.analyte_reference is not None:
            parts.append(f"{self.analyte_reference}@")

        # Ion type
        parts.append(self.ion_type.format())

        # Neutral losses
        if self.neutral_losses:
            parts.extend(str(loss) for loss in self.neutral_losses)

        # Isotopes
        if self.isotopes:
            parts.extend(str(iso) for iso in self.isotopes)

        # Adducts
        if self.adducts:
            adduct_str = "".join(
                f"+{a}" if not a.startswith(("+", "-")) else a for a in self.adducts
            )
            if adduct_str.startswith("+"):
                adduct_str = adduct_str[1:]
            parts.append(f"[M{adduct_str}]")

        # Charge
        if self.charge not in (0, 1):
            parts.append(f"^{abs(self.charge)}")

        # Mass error
        if self.mass_error:
            parts.append(f"/{self.mass_error}")

        # Confidence
        if self.confidence is not None:
            parts.append(f"*{self.confidence}")

        result = "".join(parts)
        return f"&{result}" if self.is_auxiliary else result

    def to_json(self) -> JSONDict:
        data = {
            "molecule_description": {
                **self.ion_type.to_json(),
                "series_label": self._get_series_label(),
            },
            "charge": self.charge,
            "neutral_losses": [str(loss) for loss in self.neutral_losses],
            "isotope": [iso.to_json() for iso in self.isotopes],
            "adducts": self.adducts,
        }

        if self.analyte_reference is not None:
            data["analyte_reference"] = self.analyte_reference
        if self.mass_error:
            data["mass_error"] = self.mass_error.to_json()
        if self.confidence is not None:
            data["confidence"] = self.confidence

        return data

    def _get_series_label(self) -> str:
        type_map = {
            PeptideIon: "peptide",
            InternalFragment: "internal",
            PrecursorIon: "precursor",
            ImmoniumIon: "immonium",
            ReferenceIon: "reference",
            FormulaIon: "formula",
            SMILESIon: "smiles",
            NamedCompound: "named_compound",
            UnknownIon: "unannotated",
        }
        return type_map.get(type(self.ion_type), "unknown")

    def neutral_mass(self) -> float:
        """
        Calculate the neutral mass of the annotated ion.

        Returns
        -------
        float
            The neutral mass including modifications, losses, and isotopes
        """
        # Base ion mass
        mass = self.ion_type.mass()

        # Apply neutral losses/gains
        for loss in self.neutral_losses:
            mass += self._neutral_loss_mass(loss)

        # Apply isotope shifts (simplified - actual isotope mass deltas vary)
        for iso in self.isotopes:
            mass += self._isotope_mass_shift(iso)

        # Apply adducts
        for adduct in self.adducts:
            mass += self._adduct_mass(adduct)

        return mass

    def mz(self) -> float:
        """
        Calculate the m/z value for the charged ion.

        Returns
        -------
        float
            The mass-to-charge ratio
        """
        if self.charge == 0:
            raise ValueError("Cannot calculate m/z for charge state 0")

        neutral = self.neutral_mass()
        # Add/remove proton mass based on charge
        return (neutral + (self.charge * PROTON_MASS)) / abs(self.charge)

    def _neutral_loss_mass(self, loss: NeutralLoss, monoisotopic: bool = True) -> float:
        """Calculate mass contribution of a neutral loss."""
        if loss.loss_type == NeutralLossType.FORMULA:
            result = chem_mass(loss.name)
            base_mass = result.mass if hasattr(result, "mass") else float(result)
        elif loss.loss_type == NeutralLossType.BRACED_NAME:
            base_mass = mod_mass(loss.name)
        elif loss.loss_type == NeutralLossType.REFERENCE:
            base_mass = ReferenceMolecule.get(loss.name).mass(monoisotopic=monoisotopic)
        else:
            raise ValueError(f"Unknown neutral loss type: {loss.loss_type}")

        return loss.coefficient * base_mass

    def _isotope_mass_shift(self, iso: IsotopeLabel) -> float:
        """Calculate mass shift from isotope specification."""
        if iso.element and iso.nucleon_count:
            # Element-specific isotope (e.g., 13C, 15N)
            # Use the actual mass difference for the specific isotope
            if iso.element in ISOTOPIC_ATOMIC_MASSES:
                isotopes = ISOTOPIC_ATOMIC_MASSES[iso.element]
                if iso.nucleon_count in isotopes:
                    # Get the mass of the specific isotope
                    isotope_mass = isotopes[iso.nucleon_count]
                    # Get the mass of the most abundant isotope (typically lowest A)
                    common_mass = min(isotopes.values())
                    # Mass shift per isotope substitution
                    delta = isotope_mass - common_mass
                    return iso.offset * delta

            # Fallback to neutron mass
            return iso.offset * NEUTRON_MASS
        elif iso.is_average:
            # Average isotopologue - no simple mass shift
            return 0.0
        else:
            # Generic isotope offset - use neutron mass
            return iso.offset * NEUTRON_MASS

    def _adduct_mass(self, adduct: str) -> float:
        """Calculate mass contribution of an adduct."""
        # Parse adduct string (e.g., "+Na", "-H2O", "+2H")

        # Handle multiplier
        match = re.match(r"([+-]?)(\d*)(.+)", adduct)
        if not match:
            raise ValueError(f"Invalid adduct format: {adduct}")

        sign_str, mult_str, formula = match.groups()
        sign = -1 if sign_str == "-" else 1
        mult = int(mult_str) if mult_str else 1

        # Calculate formula mass
        result = chem_mass(formula)
        mass = result.mass if hasattr(result, "mass") else float(result)

        return sign * mult * mass


class MzpafAnnotationParser:
    """Parser for mzPAF annotation strings."""

    def __init__(self, pattern: re.Pattern = ANNOTATION_PATTERN):
        self.pattern = pattern

    def parse(
        self, annotation_str: str, wrap_errors: bool = False
    ) -> List[FragmentAnnotation]:
        """Parse annotation string, potentially containing multiple annotations."""
        if not annotation_str:
            return []

        try:
            return self._parse_annotations(annotation_str)
        except Exception as e:
            if wrap_errors:
                # Return a special error annotation
                return [self._create_error_annotation(annotation_str, str(e))]
            raise

    def _parse_annotations(self, annotation_str: str) -> List[FragmentAnnotation]:
        """Parse potentially multiple comma-separated annotations."""
        annotations = []
        remaining = annotation_str

        while remaining:
            annotation, remaining = self._parse_single(remaining)
            annotations.append(annotation)

            if remaining:
                if not remaining.startswith(","):
                    raise ValueError(f"Expected ',' but found: {remaining}")
                remaining = remaining[1:]  # Skip comma

        # Validate total confidence
        total_conf = sum(a.confidence or 0 for a in annotations)
        if total_conf > 1.001:  # Small tolerance for floating point
            raise ValueError(f"Total confidence {total_conf} exceeds 1.0")

        return annotations

    def _parse_single(self, annotation_str: str) -> tuple[FragmentAnnotation, str]:
        """Parse a single annotation and return it with remaining string."""
        match = self.pattern.search(annotation_str)
        if not match:
            raise ValueError(f"Invalid mzPAF annotation: {annotation_str}")

        groups = match.groupdict()

        # Parse components
        ion_type = self._parse_ion_type(groups)
        neutral_losses = NeutralLoss.parse_many(groups.get("neutral_losses", ""))
        isotopes = self._parse_isotopes(groups.get("isotope", ""))
        adducts = self._parse_adducts(groups.get("adducts", ""))
        charge = self._parse_charge(groups.get("charge"))
        analyte_ref = (
            int(groups["analyte_reference"])
            if groups.get("analyte_reference")
            else None
        )
        mass_error = self._parse_mass_error(groups)
        confidence = float(groups["confidence"]) if groups.get("confidence") else None
        is_auxiliary = bool(groups.get("is_auxiliary"))

        annotation = FragmentAnnotation(
            ion_type=ion_type,
            neutral_losses=neutral_losses,
            isotopes=isotopes,
            adducts=adducts,
            charge=charge,
            analyte_reference=analyte_ref,
            mass_error=mass_error,
            confidence=confidence,
            is_auxiliary=is_auxiliary,
        )

        remaining = annotation_str[match.end() :]
        return annotation, remaining

    def _parse_ion_type(self, groups: dict) -> Any:
        """Parse the ion type from regex groups."""
        if groups.get("series"):
            return PeptideIon(
                series=groups["series"],
                position=int(groups["ordinal"]),
                sequence=groups.get("sequence_ordinal"),
            )

        elif groups.get("series_internal"):
            backbone = groups.get("backbone_cleavage", "")
            n_term = backbone[0] if len(backbone) >= 1 else "y"
            c_term = backbone[1] if len(backbone) >= 2 else "b"

            return InternalFragment(
                start=int(groups["internal_start"]),
                end=int(groups["internal_end"]),
                n_terminal_cleavage=n_term,
                c_terminal_cleavage=c_term,
                sequence=groups.get("sequence_internal"),
            )

        elif groups.get("precursor"):
            # Precursor can have an inline sequence specification
            # The regex pattern would need to capture this, but for now
            # we'll check for sequence in a similar way to other ion types
            return PrecursorIon(sequence=groups.get("sequence_precursor"))

        elif groups.get("immonium"):
            return ImmoniumIon(
                amino_acid=groups["immonium"],
                modification=groups.get("immonium_modification"),
            )

        elif groups.get("reference"):
            return ReferenceIon(name=groups["reference_label"])

        elif groups.get("formula"):
            return FormulaIon(formula=groups["formula"])

        elif groups.get("smiles"):
            return SMILESIon(smiles=groups["smiles"])

        elif groups.get("named_compound"):
            return NamedCompound(name=groups["named_compound"])

        elif groups.get("unannotated"):
            label = (
                int(groups["unannotated_label"])
                if groups.get("unannotated_label")
                else None
            )
            return UnknownIon(label=label)

        raise ValueError("Could not determine ion type")

    def _parse_isotopes(self, isotope_str: str) -> List[IsotopeLabel]:
        """Parse isotope specifications."""
        if not isotope_str:
            return []

        isotopes = []
        for match in ISOTOPE_PATTERN.finditer(isotope_str):
            iso = IsotopeLabel.from_string(match.group())
            if iso:
                isotopes.append(iso)
        return isotopes

    def _parse_adducts(self, adduct_str: str) -> List[str]:
        """Parse adduct specifications."""
        if not adduct_str:
            return []

        # Remove 'M' prefix if present
        if adduct_str.startswith("M"):
            adduct_str = adduct_str[1:]

        # Split on + or - (keeping the signs)
        parts = re.split(r"(?=[+-])", adduct_str)
        return [p for p in parts if p]

    def _parse_charge(self, charge_str: Optional[str]) -> int:
        """Parse charge state."""
        if not charge_str:
            return 1
        charge = int(charge_str)
        if charge == 0:
            raise ValueError("Charge cannot be zero")
        return charge

    def _parse_mass_error(self, groups: dict) -> Optional[MassError]:
        """Parse mass error specification."""
        if not groups.get("mass_error"):
            return None

        unit = groups.get("mass_error_unit", None)
        if unit is None:
            unit = "Da"
        return MassError(value=float(groups["mass_error"]), unit=unit)

    def _create_error_annotation(self, content: str, error: str) -> FragmentAnnotation:
        """Create a special annotation for parsing errors."""
        return FragmentAnnotation(
            ion_type=UnknownIon(), mass_error=MassError(0.0), is_auxiliary=True
        )


# Convenience functions
def parse_mzpaf_annotation(
    annotation_str: str, wrap_errors: bool = False
) -> List[FragmentAnnotation]:
    """Parse an mzPAF annotation string."""
    parser = MzpafAnnotationParser()
    return parser.parse(annotation_str, wrap_errors)


def format_mzpaf_annotation(annotation: FragmentAnnotation) -> str:
    """Format an annotation to mzPAF string."""
    return str(annotation)
