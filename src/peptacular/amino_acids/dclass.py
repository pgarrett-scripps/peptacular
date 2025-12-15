from collections import Counter
from dataclasses import dataclass
from functools import cached_property

from ..elements.lookup import parse_composition
from ..elements import ElementInfo


@dataclass(frozen=True)
class AminoAcidInfo:
    """Information about an amino acid"""

    id: str
    name: str
    three_letter_code: str
    formula: str | None
    monoisotopic_mass: float | None
    average_mass: float | None
    dict_composition: dict[str, int] | None
    is_mass_ambiguous: bool = False  # L / I are ambiguous but not mass ambiguous
    is_ambiguous: bool = False

    @cached_property
    def composition(self) -> Counter[ElementInfo] | None:
        """Get the composition as a Counter"""
        return (
            Counter(parse_composition(dict(self.dict_composition)))
            if self.dict_composition is not None
            else None
        )

    @property
    def one_letter_code(self) -> str:
        return self.id

    def get_mass(self, monoisotopic: bool = True) -> float | None:
        """Get the mass of the amino acid"""
        if monoisotopic:
            return self.monoisotopic_mass
        else:
            return self.average_mass

    def to_dict(self, float_format: str = "{:.6f}") -> dict[str, object]:
        """Convert the AminoAcidInfo to a dictionary"""
        return {
            "id": self.id,
            "name": self.name,
            "three_letter_code": self.three_letter_code,
            "formula": self.formula,
            "monoisotopic_mass": float_format.format(self.monoisotopic_mass)
            if self.monoisotopic_mass is not None
            else None,
            "average_mass": float_format.format(self.average_mass)
            if self.average_mass is not None
            else None,
            "composition": self.dict_composition,
        }
