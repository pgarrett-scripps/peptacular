from enum import StrEnum
from typing import Literal, NamedTuple, Protocol, Self

from ..constants import PROTON_MASS, IonType, IonTypeLiteral
from ..funcs import get_label, get_number


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


class Fragment(NamedTuple):
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
    mz: float
    sequence: str
    parent_sequence: str
    parent_length: int
    internal: bool

    @property
    def priority(self) -> int:
        """
        Returns the priority of the fragment based on its position and type.
        Higher values indicate higher priority for display/analysis.
        """
        priority = 10  # Base priority

        # Penalize internal fragments (less common/useful)
        if self.internal:
            priority -= 50

        # Penalize neutral losses (less abundant)
        if self.loss != 0.0:
            priority -= 15

        # Penalize isotope peaks (less abundant than monoisotopic)
        if self.isotope != 0:
            priority -= 35

        return priority

    @property
    def number(self) -> str:
        """
        Returns the number of the fragment, e.g., 2 for b2, 3 for y3, etc.

        :return: Number of the fragment.
        :rtype: str
        """
        return get_number(self.ion_type, self.parent_length, self.start, self.end)

    @property
    def label(self) -> str:
        """
        Returns the label of the fragment, e.g., b2, y3i, etc.

        :return: Label of the fragment.
        :rtype: str
        """
        return get_label(
            self.ion_type, self.charge, self.number, self.loss, self.isotope
        )

    @property
    def mass(self) -> float:
        return self.mz * self.charge

    @property
    def neutral_mass(self) -> float:
        return self.mz * self.charge - (self.charge * PROTON_MASS)


FRAGMENT_RETURN_TYPING = Fragment | float | str | tuple[float, str]
