from dataclasses import dataclass
from enum import StrEnum
from functools import cached_property
from typing import Literal, Protocol, Self

from ..constants import IonType, IonTypeLiteral
from ..utils2 import get_label, get_number


class FragmentableAnnotation(Protocol):
    """Protocol defining the interface required for fragmentation."""

    @property
    def stripped_sequence(self) -> str: ...

    def contains_sequence_ambiguity(self) -> bool: ...

    def split(self) -> list[Self]: ...

    def slice(
        self,
        start: int | None,
        stop: int | None,
        inplace: bool = False,
    ) -> Self: ...

    def __len__(self) -> int: ...

    def serialize(
        self,
        include_plus: bool = False,
        precision: int | None = None,
    ) -> str: ...

    def set_charge(self, charge: int) -> Self: ...

    def mass(
        self,
        ion_type: IonTypeLiteral | IonType = IonType.PRECURSOR,
        monoisotopic: bool = True,
        isotope: int = 0,
        loss: float = 0.0,
        use_isotope_on_mods: bool = False,
        precision: int | None = None,
    ) -> float: ...


class FragmentReturnType(StrEnum):
    FRAGMENT = "fragment"
    MASS = "mass"
    MZ = "mz"
    LABEL = "label"
    MASS_LABEL = "mass-label"
    MZ_LABEL = "mz-label"


FragmentReturnLiteral = Literal[
    "fragment", "mass", "mz", "label", "mass-label", "mz-label"
]
FragmentLiteral = Literal["fragment"]
MassLiteral = Literal["mass"]
MzLiteral = Literal["mz"]
LabelLiteral = Literal["label"]
MassLabelLiteral = Literal["mass-label"]
MzLabelLiteral = Literal["mz-label"]


@dataclass(frozen=True)
class Fragment:
    """
    A dataclass for representing a peptide fragment ion.
    """

    charge: int
    ion_type: str
    start: int
    end: int
    monoisotopic: bool
    isotope: int
    loss: float
    parent_sequence: str
    mass: float
    neutral_mass: float
    mz: float
    sequence: str
    unmod_sequence: str
    internal: bool

    @cached_property
    def number(self) -> str:
        """
        Returns the number of the fragment, e.g., 2 for b2, 3 for y3, etc.

        :return: Number of the fragment.
        :rtype: str
        """
        return get_number(
            self.ion_type, len(self.parent_sequence), self.start, self.end
        )

    @cached_property
    def label(self) -> str:
        """
        Returns the label of the fragment, e.g., b2, y3i, etc.

        :return: Label of the fragment.
        :rtype: str
        """
        return get_label(
            self.ion_type, self.charge, self.number, self.loss, self.isotope
        )

    def __iter__(self):
        # Include regular attributes
        for key, value in self.__dict__.items():
            yield key, value

        # Explicitly include cached properties
        yield "label", self.label
        yield "number", self.number

    def to_dict(self):
        """
        Convert the Fragment object to a dictionary, including cached properties.
        """
        return dict(self)


FRAGMENT_RETURN_TYPING = Fragment | float | str | tuple[float, str]
