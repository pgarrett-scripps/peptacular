"""Reference molecule information dataclass"""

from dataclasses import dataclass
from functools import cached_property
from collections import Counter
from ...elements import ElementInfo, parse_composition


@dataclass(frozen=True)
class RefMolInfo:
    """Information about a reference molecule"""

    name: str
    label_type: str
    molecule_type: str
    chemical_formula: str
    monoisotopic_mass: float  # calculated from formula
    average_mass: float  # calculated from formula
    dict_composition: dict[str, int]

    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get the mass of the molecule"""
        return self.monoisotopic_mass if monoisotopic else self.average_mass

    @cached_property
    def composition(self) -> Counter[ElementInfo]:
        """Get the composition as a Counter"""
        return Counter(parse_composition(dict(self.dict_composition)))

    def to_dict(self, float_format: str = "{:.6f}") -> dict[str, object]:
        """Convert the RefMolInfo to a dictionary"""
        return {
            "name": self.name,
            "label_type": self.label_type,
            "molecule_type": self.molecule_type,
            "chemical_formula": self.chemical_formula,
            "monoisotopic_mass": float_format.format(self.monoisotopic_mass),
            "average_mass": float_format.format(self.average_mass),
            "composition": self.dict_composition,
        }
