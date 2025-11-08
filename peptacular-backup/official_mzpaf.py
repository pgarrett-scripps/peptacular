"""
An implementation of the mzPAF product ion annotation format.

This reference implementation includes parsing the compact text format
and serializing.
"""

import re
import warnings
from dataclasses import dataclass
from enum import Flag
from enum import auto as enauto
from sys import intern
from typing import Any, Dict, List, Optional, Pattern, Tuple, Union

from ..chem import mod_comp
from ..mass_calc import mod_mass
from ..proforma import parse
from ..proforma.mass_comp import chem_mass
from .reference import ReferenceMolecule

JSONDict = Dict[str, Union[List, Dict, int, float, str, bool, None]]

annotation_pattern = re.compile(
    r"""
^(?P<is_auxiliary>&)?
   (?:(?P<analyte_reference>\d+)@)?
   (?:(?:(?P<series>(?:da|db|wa|wb)|[axbyczdwv]\.?)(?P<ordinal>\d+)(?:\{(?P<sequence_ordinal>.+)\})?)|
   (?P<series_internal>[m](?P<internal_start>\d+):(?P<internal_end>\d+)(?:\{(?P<sequence_internal>.+)\})?)|
   (?P<precursor>p)|
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


# (?:([+-]\d*)i(:?\d+(:?[A-Z][a-z]*))?)*
isotope_pattern = re.compile(
    r"(?P<isotope>[+-]\d*)i(?:(?P<nucleon_count>\d+)(?P<element>[A-Z][a-z]?)|(?P<average_isotopologue>A))?"
)


neutral_loss_pattern = re.compile(
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


def _sre_to_ecma(pattern):
    # Assumes that expected whitespace matches are denoted with \s
    return pattern.replace("?P<", "?<").replace("\n", "").replace(" ", "")


def tokenize_signed_symbol_list(string: str) -> List[str]:
    """
    Parse a string containing signed lists of symbols into
    a signed token list

    Parameters
    ----------
    string: str
        The symbol sequence to parse

    Returns
    -------
    list
    """
    if not string:
        return []
    tokens = re.split(r"\s*(\+|-)\s*", string)
    if not tokens:
        return []

    result = []
    sign = "+"
    symbol = None
    for token in tokens:
        if not token.strip():
            continue
        if token in ("+", "-"):
            if symbol is not None:
                result.append(sign + symbol if sign == "-" else symbol)
            sign = token
            symbol = None
        else:
            symbol = token
    if symbol is not None:
        result.append(sign + symbol if sign == "-" else symbol)
    return result


def combine_formula(tokens: List[str], leading_sign: bool = False) -> str:
    """
    Combine the tokens of one or more formulae into a parseable string.

    Parameters
    ----------
    tokens : list[str]
        The formula tokens to merge, should be composed of valid chemical group
        notation, optionally with leading signs (+/-).
    leading_sign : bool
        Whether or not to enforce that the formula must have a leading sign, adding
        the default ``+`` sign if no leading sign is present.

    Returns
    -------
    str
    """
    if not tokens:
        return ""
    if not isinstance(tokens[0], str):
        tokens = [str(t) for t in tokens]
    if not tokens[0].startswith("-") and leading_sign:
        out = ["+" + tokens[0]]
    else:
        out = [tokens[0]]
    for token in tokens[1:]:
        if token.startswith("-") or token.startswith("+"):
            out.append(token)
        else:
            out.append("+" + token)
    return "".join(out)


class NeutralNameType(Flag):
    """:class:`NeutralName` classification"""

    Reference = enauto()
    Formula = enauto()
    BracedName = enauto()
    Unknown = enauto()


@dataclass
class NeutralName(object):
    """
    Describe an entity that may appear in the neutral loss list.

    Attributes
    ----------
    name :  str
        The name of the neutral loss which may be a formula or braced name
        denoting either a "reference molecule" or a UNIMOD name
    delta_type : :class:`NeutralNameType`
        The kind of thing :attr:`name` denotes
    coefficient : int
        A signed multiplicative coefficient of the entity
    """

    name: str
    delta_type: NeutralNameType = NeutralNameType.Unknown
    coefficient: int = 1

    def __post_init__(self):
        self.delta_type = self._infer_type()

    def format_name(self, leading_sign: bool = True) -> str:
        name = self.name
        if (
            self.delta_type == NeutralNameType.Reference
            or self.delta_type == NeutralNameType.BracedName
        ):
            name = f"[{name}]"
        if self.coefficient >= 0 and leading_sign:
            if self.coefficient > 1:
                return f"+{self.coefficient}{name}"
            else:
                return f"+{name}"
        elif self.coefficient < 0:
            if self.coefficient < -1:
                return f"{self.coefficient}{name}"
            else:
                return f"-{name}"
        else:
            if self.coefficient > 1:
                return f"{self.coefficient}{name}"
            else:
                return f"{name}"

    def __str__(self):
        return self.format_name()

    def _infer_type(self):
        if self.name.startswith("[") and self.name.endswith("]"):
            inner_name = self.name = self.name[1:-1]
            if ReferenceMolecule.is_reference(inner_name):
                tp = NeutralNameType.Reference
            else:
                tp = NeutralNameType.BracedName
            self.name = inner_name
            return tp
        else:
            return NeutralNameType.Formula

    def mass(self) -> float:
        """
        Calculcate the neutral mass.

        This may involve loading the UNIMOD controlled vocabulary from a remote
        source. See :mod:`pyteomics.proforma` for more information.
        """
        mass: float = 0.0
        if self.delta_type == NeutralNameType.Formula:
            mass = chem_mass(self.name).mass
        elif self.delta_type == NeutralNameType.Reference:
            mass = ReferenceMolecule.get(self.name[1:-1]).neutral_mass
        elif self.delta_type == NeutralNameType.BracedName:
            mass = mod_mass(self.name[1:-1])
        else:
            raise ValueError(
                f"Cannot interpret {self.name} with type {self.delta_type}"
            )
        return self.coefficient * mass

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, str):
            if other.startswith("+"):
                return self.format_name(True)
            else:
                return self.format_name(False) == other
        return self.name == other.name and self.coefficient == other.coefficient

    @classmethod
    def parse(cls, string: str) -> List["NeutralName"]:
        """
        Parse a string of expressions into a list of :class:`NeutralName` using the tokenizer
        from the main parser.
        """
        if not string:
            return []
        names = []
        for match in neutral_loss_pattern.finditer(string):
            groups = match.groupdict()
            coef = int(groups["coefficient"] or 1)
            sign = groups["sign"] or "+"
            if sign == "-":
                coef = -1 * coef
            names.append(
                NeutralName(
                    groups["formula"] or groups["braced_name"], coefficient=coef
                )
            )
        return names

    @classmethod
    def combine(cls, tokens: List["NeutralName"], leading_sign: bool) -> str:
        """
        Combine a list of :class:`NeutralName` into an expression string that is compatible
        with :meth:`parse`.
        """
        if not tokens:
            return ""
        if tokens[0].coefficient >= 0 and leading_sign:
            out = [tokens[0].format_name(leading_sign=leading_sign)]
        else:
            out = [tokens[0].format_name(leading_sign=False)]
        for token in tokens[1:]:
            out.append(token.format_name(leading_sign=True))
        return "".join(out)


@dataclass
class _IsotopeVariant:
    isotope: int
    element: Optional[str] = None
    nucleon_count: Optional[int] = None
    average: bool = False

    def to_dict(self) -> Dict[str, Any]:
        state = {
            "isotope": self.isotope,
        }

        if self.average:
            state["variant"] = {"averaged": True}
        elif self.element:
            state["variant"] = {
                "element": self.element,
                "nucleon_count": self.nucleon_count,
            }
        return state

    @classmethod
    def from_dict(cls, state: Union[int, Dict[str, Any]]):
        if isinstance(state, int):
            return cls(isotope=state)
        isotope = state["isotope"]
        if state.get("averaged"):
            return cls(isotope=int(isotope), average=True)
        element = state.get("element")
        nucleon_count = state.get("nucleon_count")
        if element and nucleon_count:
            return cls(
                isotope=int(isotope), element=element, nucleon_count=int(nucleon_count)
            )
        return cls(isotope=int(isotope))

    @classmethod
    def parse(cls, string: str) -> Optional["IsotopeVariant"]:
        hit = isotope_pattern.search(string)
        if hit:
            return cls.from_parsed(hit.groupdict())
        return None

    @classmethod
    def from_parsed(cls, state: dict) -> "IsotopeVariant":
        isotope = int_or_sign(state.get("isotope", ""))
        nucleon_count = state.get("nucleon_count")
        if nucleon_count:
            nucleon_count = int(nucleon_count)
        else:
            nucleon_count = None
        element = state.get("element")
        average = bool(state.get("average_isotopologue"))
        return cls(isotope, element, nucleon_count, average)

    def format_name(self, leading_sign: bool = True) -> str:
        if self.element:
            element_suffix = f"{self.nucleon_count}{self.element}"
        elif self.average:
            element_suffix = "A"
        else:
            element_suffix = ""
        if self.isotope >= 0 and leading_sign:
            if self.isotope > 1:
                return f"+{self.isotope}i{element_suffix}"
            else:
                return f"+i{element_suffix}"
        elif self.isotope < 0:
            if self.isotope < -1:
                return f"{self.isotope}i{element_suffix}"
            else:
                return f"-i{element_suffix}"
        else:
            if self.isotope > 1:
                return f"{self.isotope}i{element_suffix}"
            else:
                return f"i{element_suffix}"

    def __str__(self):
        return self.format_name()


class IsotopeVariant(_IsotopeVariant):
    def __eq__(self, other: "IsotopeVariant"):
        if other is None:
            return False
        if isinstance(other, int):
            other = IsotopeVariant(other)
        elif isinstance(other, str):
            other = IsotopeVariant.parse(other)
        return super().__eq__(other)

    def __ne__(self, other):
        return not self == other


class MassError(object):
    """
    Represent the mass error of a peak annotation.

    Attributes
    ----------
    unit : str
        The unit of the mass error, may be Da (daltons) or ppm (parts-per-million)
    mass_error : float
        The magnitude of the mass error
    """

    _DEFAULT_UNIT = "Da"

    unit: str
    mass_error: float

    def __init__(self, mass_error=0, unit=None):
        if unit is None:
            unit = self._DEFAULT_UNIT
        self.mass_error = float(mass_error)
        self.unit = unit

    def serialize(self) -> str:
        """
        Serialize the mass error to the notation used in mzPAF.

        Returns
        -------
        str
        """
        unit = self.unit
        if unit == self._DEFAULT_UNIT:
            unit = ""
        return f"{self.mass_error}{unit}"

    def __eq__(self, other):
        return self.mass_error == other.mass_error and self.unit == other.unit

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return self.serialize()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass_error}, {self.unit})"

    def to_json(self) -> JSONDict:
        """
        Serialize the mass error to a JSON-compatible :class:`dict`

        Returns
        -------
        dict
        """
        return {"value": self.mass_error, "unit": self.unit}


class _HasSequenceMixin:
    __slots__ = ()
    sequence: str

    def has_sequence(self) -> bool:
        return bool(self.sequence)

    def parse_sequence(self):
        if not self.has_sequence():
            raise ValueError(
                "This annotation has no specified sequence. See the referenced analyte to get the sequence."
            )
        return parse(self.sequence)


class _SeriesLabelSubclassRegisteringMeta(type):
    """
    A metaclass which does some custom registration whenever a subclass is defined.

    Attributes
    ----------
    _label_registry : :class:`dict` mapping :class:`str` to implementing :class:`_SeriesLabelSubclassRegisteringMeta` types.
    """

    _label_registry: Dict[str, "_SeriesLabelSubclassRegisteringMeta"]

    def __new__(mcls, name, bases, attrs):
        label = attrs.get("series_label")
        if label and isinstance(label, str):
            label = intern(label)
        override_label = attrs.get("override_label", False)
        cls = super(_SeriesLabelSubclassRegisteringMeta, mcls).__new__(
            mcls, name, bases, attrs
        )
        if not hasattr(cls, "_label_registry"):
            registry = cls._label_registry = {}
        else:
            registry = cls._label_registry
        if label:
            if label not in registry or override_label:
                registry[label] = cls
        return cls


class IonAnnotationBase(object, metaclass=_SeriesLabelSubclassRegisteringMeta):
    """
    A base class for all mzPAF ion annotations.

    Provides shared functionality and common attributes for the data model.
    """

    __slots__ = (
        "series",
        "neutral_losses",
        "isotope",
        "adducts",
        "charge",
        "analyte_reference",
        "mass_error",
        "confidence",
        "rest",
        "is_auxiliary",
    )

    series_label = None
    _molecule_description_fields = {}

    series: str
    neutral_losses: List
    isotope: List[IsotopeVariant]
    adducts: List
    charge: int
    analyte_reference: Optional[int]
    mass_error: MassError
    confidence: float
    rest: Any
    is_auxiliary: bool

    def __init__(
        self,
        series,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        if isotope is None:
            isotope = 0
        if charge is None:
            charge = 1

        self.series = series
        self.neutral_losses = neutral_losses or []
        self.isotope = isotope
        self.adducts = adducts or []
        self.charge = charge
        self.analyte_reference = (
            int(analyte_reference) if analyte_reference is not None else None
        )
        self.mass_error = mass_error
        self.confidence = confidence
        self.rest = rest
        self.is_auxiliary = is_auxiliary

    def __hash__(self):
        return hash(self.serialize())

    def __eq__(self, other):
        return self.serialize() == str(other)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return self.serialize()

    def _format_ion(self) -> str:
        raise NotImplementedError()

    def serialize(self) -> str:
        """
        Serialize the data model back to an mzPAF formatted string.

        Returns
        -------
        str
        """
        parts = []
        if self.analyte_reference is not None:
            parts.append(f"{self.analyte_reference}@")
        parts.append(self._format_ion())
        if self.neutral_losses:
            parts.append(NeutralName.combine(self.neutral_losses, leading_sign=True))
        if self.isotope:
            parts.extend(map(str, self.isotope))
        if self.adducts:
            parts.append(
                "[{}]".format(combine_formula(self.adducts, leading_sign=False))
            )
        if self.charge != 0 and self.charge != 1:
            charge = abs(self.charge)
            parts.append(f"^{charge}")
        if self.mass_error is not None:
            parts.append("/")
            parts.append(self.mass_error.serialize())
        if self.confidence is not None:
            parts.append(f"*{self.confidence}")
        if self.rest is not None:
            parts.append("/")
            parts.append(self.rest)
        result = "".join(parts)
        if self.is_auxiliary:
            return f"&{result}"
        return result

    def __str__(self):
        return self.serialize()

    def _molecule_description(self) -> JSONDict:
        return {"series_label": self.series_label}

    def to_json(self, exclude_missing=False) -> JSONDict:
        """
        Convert the data model into a JSON-serializable :class:`dict`.

        Returns
        -------
        dict
        """
        # TODO: When neutral losses and adducts are formalized types, convert to string/JSON here
        d = {}
        skips = ("series", "rest", "is_auxiliary", "neutral_losses", "isotope")
        for key in IonAnnotationBase.__slots__:
            if key in skips:
                continue
            if key == "mass_error" and self.mass_error is not None:
                d[key] = self.mass_error.to_json()
            else:
                value = getattr(self, key)
                if (value is not None) or not exclude_missing:
                    d[key] = value
        d["neutral_losses"] = [str(s) for s in self.neutral_losses]
        d["isotope"] = [s.to_dict() for s in self.isotope if s]
        d["molecule_description"] = self._molecule_description()
        # if d['analyte_reference'] is None:
        #     d['analyte_reference'] =
        return d

    def _populate_from_dict(self, data) -> "IonAnnotationBase":
        # TODO: When neutral losses and adducts are formalized types, parse from string here
        for key, value in data.items():
            if key == "molecule_description":
                continue
            elif key == "mass_error" and value is not None:
                self.mass_error = MassError(value["value"], value["unit"])

            elif key == "isotope" and value is not None:
                if isinstance(value, str):
                    self.isotope.append(IsotopeVariant.parse(value))
                elif isinstance(value, int):
                    self.isotope.append(IsotopeVariant(value))
                elif isinstance(value, (list, tuple)):
                    self.isotope = []
                    for tok in value:
                        if isinstance(tok, str):
                            self.isotope.append(IsotopeVariant.parse(tok))
                        else:
                            self.isotope.append(IsotopeVariant.from_dict(tok))
                else:
                    warnings.warn(f"Failed to coerce {value} to isotopic variant")
            elif key == "neutral_losses" and value is not None:
                if isinstance(value, str):
                    self.neutral_losses = NeutralName.parse(value)
                elif isinstance(value, (list, tuple)):
                    self.neutral_losses = []
                    for tok in value:
                        self.neutral_losses.extend(NeutralName.parse(tok))
                else:
                    self.neutral_losses = []
                    warnings.warn(f"Failed to coerce {value} to neutral losses")
            else:
                setattr(self, key, value)
        self.rest = None
        self.is_auxiliary = False
        return self

    @classmethod
    def from_json(cls, data) -> "IonAnnotationBase":
        """
        Convert parsed JSON data back into formal data structures.

        Returns
        -------
        IonAnnotationBase
        """
        descr = data["molecule_description"]
        series_label = descr["series_label"]
        cls = cls._label_registry[series_label]
        self = cls.__new__(cls)
        self._populate_from_dict(data)
        return self


class PeptideFragmentIonAnnotation(IonAnnotationBase, _HasSequenceMixin):
    """
    Represent a canonical peptide backbone fragment, such as a b or y ion in the Roepstorff notation.

    Attributes
    ----------
    series : str
        The ion series
    position : int
        The ordinal bond position relative to the series terminal for this fragment
    sequence : str or :const:`None`
        The sequence this fragment is derived from if it is not a specified analyte from
        the spectrum metadata. Should be written in ProForma 2.0 notation.
    """

    __slots__ = ("position", "sequence")

    series_label = "peptide"

    _molecule_description_fields = {
        "series": "The peptide ion series this ion belongs to",
        "position": "The position from the appropriate terminal along the peptide this ion was fragmented at (starting with 1)",
        "sequence": "An arbitrary sequence that this ion came from, if not the primary sequence of the associated analyte",
    }

    position: int
    sequence: Optional[str]

    def __init__(
        self,
        series,
        position,
        sequence=None,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        super(PeptideFragmentIonAnnotation, self).__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )
        self.position = position
        self.sequence = sequence

    def _format_ion(self) -> str:
        if not self.sequence:
            return f"{self.series}{self.position}"
        return f"{self.series}{self.position}{{{self.sequence}}}"

    def _molecule_description(self) -> JSONDict:
        d = super()._molecule_description()
        d.update(
            {
                "series": self.series,
                "position": self.position,
            }
        )
        d["sequence"] = self.sequence
        return d

    def _populate_from_dict(self, data) -> IonAnnotationBase:
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.series = descr["series"]
        self.position = descr["position"]
        self.sequence = descr.get("sequence")
        return self


class InternalPeptideFragmentIonAnnotation(IonAnnotationBase, _HasSequenceMixin):
    """
    Represent an internal fragment of a peptide, containing neither of the original
    termini of the source peptide.

    Attributes
    ----------
    start_position : int
        The bond ordinal where the first break occurs relative to the N-terminal of the peptide.
    end_position : int
        The bond ordinal where the second break occurs relative to the N-terminal of the peptide.
    sequence : str or :const:`None`
        The sequence this fragment is derived from if it is not a specified analyte from
        the spectrum metadata. Should be written in ProForma 2.0 notation.
    """

    __slots__ = ("start_position", "end_position", "sequence")

    series_label = "internal"

    _molecule_description_fields = {
        "start_position": (
            "N-terminal amino acid residue of the fragment in the "
            "original peptide sequence (beginning with 1, counting from "
            "the N-terminus)"
        ),
        "end_position": (
            "C-terminal amino acid residue of the fragment in the original"
            " peptide sequence (beginning with 1, counting from the "
            "N-terminus)"
        ),
        "sequence": (
            "An arbitrary sequence that this ion came from, if not the primary"
            " sequence of the associated analyte"
        ),
    }

    start_position: int
    end_position: int
    sequence: str

    def __init__(
        self,
        series,
        start_position,
        end_position,
        sequence=None,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        super(InternalPeptideFragmentIonAnnotation, self).__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )
        self.start_position = start_position
        self.end_position = end_position
        self.sequence = sequence

    def _format_ion(self) -> str:
        if not self.sequence:
            return f"m{self.start_position}:{self.end_position}"
        return f"m{self.start_position}:{self.end_position}{{{self.sequence}}}"

    def _molecule_description(self) -> JSONDict:
        d = super()._molecule_description()
        d["start_position"] = self.start_position
        d["end_position"] = self.end_position
        d["sequence"] = self.sequence
        return d

    def _populate_from_dict(self, data) -> IonAnnotationBase:
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.start_position = descr["start_position"]
        self.end_position = descr["end_position"]
        self.sequence = descr.get("sequence")
        return self


class PrecursorIonAnnotation(IonAnnotationBase):
    """The intact, possibly charge-reduced, precursor ion."""

    __slots__ = ()

    series_label = "precursor"
    _molecule_description_fields = {}

    def __init__(
        self,
        series,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        super(PrecursorIonAnnotation, self).__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )

    def _format_ion(self):
        return "p"


class ImmoniumIonAnnotation(IonAnnotationBase):
    """
    Represent a single amino acid ionized with a neutral loss.

    Attributes
    ----------
    amino_acid : str
        The amino acid ionized
    modification : str or :const:`None`
        A modification still attached to the amino acid. Should be denoted using a controlled vocabulary term.
    """

    __slots__ = ("amino_acid", "modification")

    series_label = "immonium"
    _molecule_description_fields = {
        "amino_acid": "The amino acid represented by this immonium ion",
        "modification": "An optional modification that may be attached to this immonium ion",
    }

    amino_acid: str
    modification: str

    def __init__(
        self,
        series,
        amino_acid,
        modification=None,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        super(ImmoniumIonAnnotation, self).__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )
        self.amino_acid = amino_acid
        self.modification = modification

    def _format_ion(self):
        if self.modification is not None:
            modification = f"[{self.modification}]"
        else:
            modification = ""
        return f"I{self.amino_acid}{modification}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d["amino_acid"] = self.amino_acid
        if self.modification:
            d["modification"] = self.modification
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.amino_acid = descr["amino_acid"]
        self.modification = descr.get("modification")
        return self


class ReferenceIonAnnotation(IonAnnotationBase):
    """
    An ion which matches a referencea entity like a TMT reporter ion or a signature ion.

    Attribute
    ---------
    reference : str
        The reference identifier.
    """

    __slots__ = ("_reference", "reference_molecule")

    series_label = "reference"
    _molecule_description_fields = {"reference": "The molecule refernce identifier"}

    _reference: str
    reference: str
    reference_molecule: Optional[ReferenceMolecule]

    def __init__(
        self,
        series,
        reference,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        super(ReferenceIonAnnotation, self).__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )
        self._reference = None
        self.reference = reference

    @property
    def reference(self) -> str:
        return self._reference

    @reference.setter
    def reference(self, value):
        self._reference = value
        if value is not None:
            try:
                self.reference_molecule = ReferenceMolecule.get(value)
            except KeyError:
                warnings.warn(f"Could not find a reference entry for {value}")
                self.reference_molecule = None
        else:
            self.reference_molecule = None

    def _format_ion(self):
        return f"r[{self.reference}]"

    def _molecule_description(self):
        d = super()._molecule_description()
        d["reference"] = self.reference
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.reference = descr["reference"]
        return self


class NamedCompoundIonAnnotation(IonAnnotationBase):
    __slots__ = ("compound_name",)

    series_label = "named_compound"

    _molecule_description_fields = {
        "compound_name": "The name of the named compound ion being marked"
    }

    compound_name: str

    def __init__(
        self,
        series,
        label,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        super(NamedCompoundIonAnnotation, self).__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )
        self.compound_name = label

    def _format_ion(self):
        return f"_{{{self.compound_name}}}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d["compound_name"] = self.compound_name
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.compound_name = descr["compound_name"]
        return self


class FormulaAnnotation(IonAnnotationBase):
    __slots__ = ("formula",)

    series_label = "formula"
    _molecule_description_fields = {
        "formula": "The elemental formula of the ion being marked"
    }

    formula: str

    def __init__(
        self,
        series,
        formula,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        super(FormulaAnnotation, self).__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )
        self.formula = formula

    def _format_ion(self):
        return f"f{{{self.formula}}}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d["formula"] = self.formula
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.formula = descr["formula"]
        return self

    def to_composition(self) -> dict:
        return mod_comp(self.formula)


class SMILESAnnotation(IonAnnotationBase):
    __slots__ = ("smiles",)

    series_label = "smiles"
    _molecule_description_fields = {
        "smiles": "The SMILES definition of the ion being marked"
    }

    smiles: str

    def __init__(
        self,
        series,
        smiles,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        super().__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )
        self.smiles = smiles

    def _format_ion(self):
        return f"s{{{self.smiles}}}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d["smiles"] = self.smiles
        return d

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.smiles = descr["smiles"]
        return self


class Unannotated(IonAnnotationBase):
    series_label = "unannotated"
    _molecule_description_fields = {
        "unannotated_label": "A user-specified numeral label for an unannotated peak"
    }

    __slots__ = ("unannotated_label",)

    unannotated_label: str

    @classmethod
    def empty(cls):
        return cls(cls.series_label, None, mass_error=MassError())

    def __init__(
        self,
        series,
        unannotated_label,
        neutral_losses=None,
        isotope=None,
        adducts=None,
        charge=None,
        analyte_reference=None,
        mass_error=None,
        confidence=None,
        rest=None,
        is_auxiliary=None,
    ):
        self.unannotated_label = unannotated_label
        super().__init__(
            series,
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
            rest,
            is_auxiliary,
        )

    def _format_ion(self):
        if not self.unannotated_label:
            return "?"
        return f"?{self.unannotated_label}"

    def _molecule_description(self):
        d = super()._molecule_description()
        d["unannotated_label"] = self.unannotated_label
        return d

    def serialize(self) -> str:
        # mass_error is a required field in the data model, but it is meaningless for unannotated peaks,
        # so we mask it here
        mass_error = self.mass_error

        mask = mass_error.mass_error == 0
        if mask:
            self.mass_error = None

        val = super().serialize()

        if mask:
            self.mass_error = mass_error
        return val

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.unannotated_label = descr.get("unannotated_label")
        return self


class InvalidAnnotation(IonAnnotationBase):
    series_label = "!invalid!"

    content: str
    error: str

    def __init__(self, content: str, error: str):
        super().__init__(self.series_label)
        self.content = content
        self.error = error

    def serialize(self):
        return self.content

    def _molecule_description(self) -> JSONDict:
        return {
            "series_label": self.series_label,
            "content": self.content,
            "error": self.error,
        }

    def _populate_from_dict(self, data):
        super()._populate_from_dict(data)
        descr = data["molecule_description"]
        self.content = descr["content"]
        self.error = descr["error"]
        return self


def int_or_sign(string: str) -> int:
    if string == "+":
        return 1
    elif string == "-":
        return -1
    else:
        return int(string)


class AnnotationStringParser(object):
    """
    An annotation string parser using a specific parser pattern.

    This class organizes the parsing of common and type-specific components of a
    string into one or more peak annotations.

    Instances are :class:`Callable`, with :obj:`parse_annotation` being the reference
    parser.
    """

    pattern: Pattern

    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(
        self, annotation_string: str, *, wrap_errors: bool = False, **kwargs
    ) -> List[IonAnnotationBase]:
        """
        Parse a string into one or more :class:`IonAnnotationBase` instances.

        Parameters
        ----------
        annotation_string : str
            The string to be parsed
        wrap_errors : bool, optional
            Whether or not to capture parsing errors as :class:`InvalidAnnotation` or not. Defaults to :const:`False`.

        Returns
        -------
        list[:class:`IonAnnotationBase`] :
            The annotations parsed
        """
        try:
            return self.parse_annotation(annotation_string, **kwargs)
        except ValueError as err:
            if wrap_errors:
                return [InvalidAnnotation(annotation_string, str(err))]
            else:
                raise

    def _parse_string(
        self, annotation_string: str, **kwargs
    ) -> Tuple[re.Match, Dict[str, str]]:
        match = self.pattern.search(annotation_string)
        if match is None:
            raise ValueError(f"Invalid annotation string {annotation_string!r}")
        data = match.groupdict()
        return match, data

    def _coerce_isotope(self, data: Dict[str, str]) -> List[IsotopeVariant]:
        items = data.get("isotope") or ""
        it = isotope_pattern.finditer(items)
        isotope = []
        for grp in it:
            isotope.append(IsotopeVariant.from_parsed(grp.groupdict()))
        return isotope

    def _coerce_charge(self, data: Dict[str, str]) -> int:
        charge = data.get("charge", 1)
        if charge is None:
            charge = 1
        elif charge == 0:
            raise ValueError(
                f"The charge of an annotation cannot be zero (parsed {data['charge']})"
            )
        else:
            charge = int(charge)
        return charge

    def _coerce_adducts(self, data: Dict[str, str]) -> List[str]:
        adducts = tokenize_signed_symbol_list(data.get("adducts"))
        return adducts

    def _coerce_analyte_reference(self, data: Dict[str, str]) -> str:
        return data.get("analyte_reference", "1")

    def _coerce_neutral_losses(self, data: Dict[str, str]) -> List:
        tokens = NeutralName.parse(data.get("neutral_losses", ""))
        return tokens

    def _coerce_mass_error(self, data: Dict[str, str]) -> MassError:
        mass_error = data.get("mass_error")
        if mass_error is not None:
            mass_error = MassError(float(mass_error), data.get("mass_error_unit"))
        return mass_error

    def _coerce_confidence(self, data: Dict[str, str]) -> float:
        confidence = data.get("confidence")
        if confidence is not None:
            confidence = float(confidence)
            if confidence > 1.0:
                raise ValueError(
                    "A single peak interpretation's confidence cannot be greater than 1."
                )
        return confidence

    def parse_annotation(
        self, annotation_string: str, **kwargs
    ) -> List[IonAnnotationBase]:
        """
        Parse a string into one or more :class:`IonAnnotationBase` instances.

        Parameters
        ----------
        annotation_string : str
            The string to be parsed
        **kwargs
            Passed to the :meth:`_dispatch` which in turn creates :class:`IonAnnotationBase`
            instances

        Returns
        -------
        list[:class:`IonAnnotationBase`] :
            The annotations parsed
        """
        if not annotation_string:
            return []

        match, data = self._parse_string(annotation_string)

        is_auxiliary = bool(data.get("is_auxiliary"))
        adducts = self._coerce_adducts(data)
        charge = self._coerce_charge(data)
        isotope = self._coerce_isotope(data)

        # FIXME: ensure that neutral loss is not a plain mass, and tokenize separate blocks
        neutral_losses = self._coerce_neutral_losses(data)

        analyte_reference = self._coerce_analyte_reference(data)
        mass_error = self._coerce_mass_error(data)
        confidence = self._coerce_confidence(data)

        annotation = self._dispatch(
            annotation_string,
            data,
            adducts,
            charge,
            isotope,
            neutral_losses,
            analyte_reference,
            mass_error,
            confidence,
            **kwargs,
        )
        if is_auxiliary:
            annotation.is_auxiliary = True

        rest = annotation_string[match.end() :]
        if rest == "":
            return [annotation]
        else:
            if rest[0] != ",":
                raise ValueError(
                    f"Malformed trailing string {rest}, expected ',' for {annotation_string}"
                )
            else:
                rest = rest[1:]
            result = [annotation]
            result.extend(self.parse_annotation(rest, **kwargs))
            total_confidence = 0.0
            for annot in result:
                if annot.confidence is not None:
                    total_confidence += annot.confidence
            if total_confidence < 0 or total_confidence > (1 + 1e-3):
                raise ValueError(
                    f"The sum of all interpretations of a single peak's confidence cannot be greater than 1"
                    f" ({total_confidence}). {annotation_string}"
                )
            return result

    def _dispatch(
        self,
        annotation_string,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        if data.get("series"):
            return self._dispatch_peptide_fragment(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        elif data.get("series_internal"):
            return self._dispatch_internal_peptide_fragment(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        elif data.get("precursor"):
            return self._dispatch_precursor(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        elif data.get("immonium"):
            return self._dispatch_immonium(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        elif data.get("reference"):
            return self._dispatch_reference(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        elif data.get("named_compound"):
            return self._dispatch_named_compound(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        elif data.get("formula"):
            return self._dispatch_formula(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        elif data.get("smiles"):
            return self._dispatch_smiles(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        elif data.get("unannotated"):
            return self._dispatch_unannotated(
                data,
                neutral_losses=neutral_losses,
                isotope=isotope,
                adducts=adducts,
                charge=charge,
                analyte_reference=analyte_reference,
                mass_error=mass_error,
                confidence=confidence,
                **kwargs,
            )
        else:
            raise ValueError(
                f"Could not infer annotation type from {annotation_string}/{data}"
            )

    def _dispatch_peptide_fragment(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        return PeptideFragmentIonAnnotation(
            data["series"],
            int(data["ordinal"]),
            data["sequence_ordinal"],
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )

    def _dispatch_unannotated(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        if mass_error is None:
            mass_error = MassError(0)
        return Unannotated(
            None,
            data.get("unannotated_label"),
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )

    def _dispatch_internal_peptide_fragment(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        return InternalPeptideFragmentIonAnnotation(
            "internal",
            int(data["internal_start"]),
            int(data["internal_end"]),
            data["sequence_internal"],
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )

    def _dispatch_precursor(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        return PrecursorIonAnnotation(
            "precursor",
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )

    def _dispatch_immonium(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        return ImmoniumIonAnnotation(
            "immonium",
            data["immonium"],
            data["immonium_modification"],
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )

    def _dispatch_reference(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        return ReferenceIonAnnotation(
            "reference",
            (data["reference_label"]),
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )

    def _dispatch_named_compound(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        return NamedCompoundIonAnnotation(
            "named_compound",
            data["named_compound"],
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )

    def _dispatch_formula(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        return FormulaAnnotation(
            "formula",
            data["formula"],
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )

    def _dispatch_smiles(
        self,
        data,
        adducts,
        charge,
        isotope,
        neutral_losses,
        analyte_reference,
        mass_error,
        confidence,
        **kwargs,
    ):
        return SMILESAnnotation(
            "smiles",
            data["smiles"],
            neutral_losses,
            isotope,
            adducts,
            charge,
            analyte_reference,
            mass_error,
            confidence,
        )


#: The reference parser
parse_annotation = AnnotationStringParser(annotation_pattern)
