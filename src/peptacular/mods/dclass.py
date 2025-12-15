from dataclasses import dataclass
from functools import cached_property
from typing import Any, TypeVar, Self

from ..elements import ElementInfo, parse_composition


T = TypeVar("T", bound="OboEntity")


@dataclass(frozen=True, slots=True)
class OboEntity:
    """Base class for OBO file entities"""

    id: str
    name: str
    formula: str | None
    monoisotopic_mass: float | None
    average_mass: float | None
    dict_composition: dict[str, int] | None

    def __str__(self) -> str:
        return f"{self.name} ({self.formula})"

    @cached_property
    def composition(self) -> dict[ElementInfo, int] | None:
        """Get the composition as a dict of element symbols to counts"""
        return (
            parse_composition(self.dict_composition) if self.dict_composition else None
        )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(id={self.id}, name={self.name}, formula={self.formula}, "
            f"monoisotopic_mass={self.monoisotopic_mass}, average_mass={self.average_mass}, "
            f"composition={self.dict_composition})"
        )

    def update(self, **kwargs: Any) -> Self:
        """Return a new instance with updated fields"""
        return self.__class__(
            id=kwargs.get("id", self.id),
            name=kwargs.get("name", self.name),
            formula=kwargs.get("formula", self.formula),
            monoisotopic_mass=kwargs.get("monoisotopic_mass", self.monoisotopic_mass),
            average_mass=kwargs.get("average_mass", self.average_mass),
            dict_composition=kwargs.get("dict_composition", self.dict_composition),
        )

    def mass(self, monoisotopic: bool = True) -> float | None:
        """Get the mass of the entity"""
        return self.monoisotopic_mass if monoisotopic else self.average_mass

    def to_dict(self, float_format: str = "{:.6f}") -> dict[str, object]:
        """Convert the OboEntity to a dictionary"""
        return {
            "id": self.id,
            "name": self.name,
            "formula": self.formula,
            "monoisotopic_mass": float_format.format(self.monoisotopic_mass)
            if self.monoisotopic_mass is not None
            else None,
            "average_mass": float_format.format(self.average_mass)
            if self.average_mass is not None
            else None,
            "composition": self.dict_composition,
        }


class MonosaccharideInfo(OboEntity):
    """Class to store information about a monosaccharide"""

    pass


class UnimodInfo(OboEntity):
    """Class to store information about a Unimod modification"""

    pass


class PsimodInfo(OboEntity):
    """Class to store information about a PSI-MOD modification"""

    pass
