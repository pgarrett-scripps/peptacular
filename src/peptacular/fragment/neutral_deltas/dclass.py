from collections import Counter
from dataclasses import dataclass
from functools import cached_property

from ...elements.dclass import ElementInfo
from ...elements.lookup import parse_composition
from ...proforma_components.comps import ChargedFormula


@dataclass(frozen=True)
class NeutralDeltaInfo:
    formula: str
    name: str
    description: str
    amino_acids: frozenset[str]  # Set of single-letter amino acid codes
    monoisotopic_mass: float
    average_mass: float
    dict_composition: dict[str, int]

    def __hash__(self) -> int:
        return hash(self.name)

    def get_mass(self, monoisotopic: bool = True) -> float:
        """Get the mass of the fragment ion"""
        if monoisotopic:
            return self.monoisotopic_mass
        else:
            return self.average_mass

    @cached_property
    def composition(self) -> Counter[ElementInfo]:
        """Get the composition as a Counter"""
        return Counter(parse_composition(dict(self.dict_composition)))

    @cached_property
    def charged_formula(self) -> ChargedFormula:
        """Get the ChargedFormula representation of the neutral delta"""
        return ChargedFormula.from_composition(self.composition)

    def calculate_loss_sites(self, sequence: str) -> int:
        """Calculate the number of possible loss sites in a sequence"""
        return sum(1 for aa in sequence if aa in self.amino_acids)

    def to_dict(self, float_precision: int = 6) -> dict[str, object]:
        """Convert to dictionary representation"""
        return {
            "formula": self.formula,
            "name": self.name,
            "description": self.description,
            "amino_acids": list(self.amino_acids),
            "monoisotopic_mass": round(self.monoisotopic_mass, float_precision),
            "average_mass": round(self.average_mass, float_precision),
            "dict_composition": dict(self.dict_composition),
        }
