from dataclasses import dataclass
from enum import Flag, auto
from functools import cached_property
from typing import Counter, Mapping
import typing

from ...elements import ElementInfo, parse_composition

# type checking
if typing.TYPE_CHECKING:
    from ...fragment import IonType


class IonTypeProperty(Flag):
    """Flag enum for ion type properties."""

    NONE = 0
    FORWARD = auto()  # a, b, c
    BACKWARD = auto()  # x, y, z
    INTERNAL = auto()  # internal fragments
    INTACT = auto()  # precursor, neutral
    AA_SPECIFIC_FWD = auto()  # d, da, db
    AA_SPECIFIC_BWD = auto()  # v, w, wa, wb


@dataclass(frozen=True)
class FragmentIonInfo:
    """Information about a fragment ion"""

    id: str
    name: str
    formula: str | None
    monoisotopic_mass: float | None
    average_mass: float | None
    dict_composition: Mapping[str, int] | None
    properties: IonTypeProperty = IonTypeProperty.NONE

    @property
    def ion_type(self) -> "IonType":
        from .data import IonType

        if isinstance(self.id, IonType):
            return self.id
        return IonType[self.id]

    @property
    def is_forward(self) -> bool:
        """Check if ion is a forward ion type (a, b, c)"""
        return bool(self.properties & IonTypeProperty.FORWARD)

    @property
    def is_backward(self) -> bool:
        """Check if ion is a backward ion type (x, y, z)"""
        return bool(self.properties & IonTypeProperty.BACKWARD)

    @property
    def is_internal(self) -> bool:
        """Check if ion is an internal fragment"""
        return bool(self.properties & IonTypeProperty.INTERNAL)

    @property
    def is_intact(self) -> bool:
        """Check if ion is an intact ion (precursor, neutral)"""
        return bool(self.properties & IonTypeProperty.INTACT)

    @property
    def is_aa_specific_forward(self) -> bool:
        """Check if ion is an amino acid-specific forward ion (d, da, db)"""
        return bool(self.properties & IonTypeProperty.AA_SPECIFIC_FWD)

    @property
    def is_aa_specific_backward(self) -> bool:
        """Check if ion is an amino acid-specific backward ion (v, w, wa, wb)"""
        return bool(self.properties & IonTypeProperty.AA_SPECIFIC_BWD)

    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get the mass of the fragment ion"""
        if monoisotopic:
            if self.monoisotopic_mass is None:
                raise ValueError("Monoisotopic mass is not available for this ion type")
            return self.monoisotopic_mass
        else:
            if self.average_mass is None:
                raise ValueError("Average mass is not available for this ion type")
            return self.average_mass

    @cached_property
    def composition(self) -> Counter[ElementInfo]:
        """Get the composition as a Counter"""
        if self.dict_composition is None:
            raise ValueError("Composition is not available for this ion type")

        return Counter(parse_composition(dict(self.dict_composition)))

    def to_dict(self, float_precision: int = 6) -> dict[str, object]:
        """Convert the FragmentIonInfo to a dictionary"""
        return {
            "id": self.id,
            "name": self.name,
            "formula": self.formula,
            "monoisotopic_mass": round(self.monoisotopic_mass, float_precision)
            if self.monoisotopic_mass is not None
            else None,
            "average_mass": round(self.average_mass, float_precision)
            if self.average_mass is not None
            else None,
            "composition": self.dict_composition,
        }
