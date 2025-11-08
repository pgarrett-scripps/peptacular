import json
from dataclasses import asdict, dataclass, field
from importlib.resources import open_text
from typing import ClassVar, Dict, List, Optional

from ..chem.chem_util import parse_chem_formula


@dataclass(frozen=True)
class ReferenceMolecule:
    name: str
    molecule_type: str
    cv_term: Optional[str] = field(default=None)
    label_type: Optional[str] = field(default=None)
    monoisotopic_mass: float = field(default=None)
    average_mass: float = field(default=None)
    ion_mz: float = field(default=None)
    neutral_mass: Optional[float] = field(default=None)
    chemical_formula: Optional[str] = field(default=None)
    ion_chemical_formula: Optional[str] = field(default=None)
    references: List[str] = field(default_factory=list)

    _registry: ClassVar[Dict[str, "ReferenceMolecule"]] = None

    def mass(self, monoisotopic: bool = True) -> float:
        if monoisotopic:
            return self.monoisotopic_mass
        else:
            return self.average_mass

    def comp(self) -> dict:
        """Get chemical composition as a dictionary."""
        if self.chemical_formula is not None:
            composition = parse_chem_formula(self.chemical_formula)
            return composition.to_dict()
        else:
            return {}

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, state, **kwargs):
        return cls(**state)

    @classmethod
    def _load_registry(cls):
        with open_text("peptacular.data", "reference_molecules.json") as stream:
            data = json.load(stream)
        cls._registry = {}
        for k, v in data.items():
            v["name"] = k
            cls._registry[k] = cls.from_dict(v)

    @classmethod
    def get(cls, name: str) -> "ReferenceMolecule":
        if cls._registry is None:
            cls._load_registry()
        return cls._registry[name]

    @classmethod
    def is_reference(cls, name: str) -> bool:
        if cls._registry is None:
            cls._load_registry()
        return name in cls._registry


def load_json(stream) -> Dict[str, ReferenceMolecule]:
    if isinstance(stream, dict):
        payload = stream
    else:
        payload = json.load(stream)
    index = {}
    for name, state in payload.items():
        state["name"] = name
        index[name] = ReferenceMolecule.from_dict(state)
    return index
