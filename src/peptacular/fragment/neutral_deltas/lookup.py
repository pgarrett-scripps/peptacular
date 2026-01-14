from collections.abc import Iterable

from .data import NEUTRAL_DELTA_DICT, NeutralDelta
from .dclass import NeutralDeltaInfo


class NeutralDeltaLookup:
    def __init__(
        self, neutral_delta_data: dict[NeutralDelta, NeutralDeltaInfo]
    ) -> None:
        self._neutral_delta_data = neutral_delta_data

        self._delta_to_data: dict[NeutralDelta, NeutralDeltaInfo] = {}
        self._formula_to_data: dict[str, NeutralDeltaInfo] = {}
        self._name_to_data: dict[str, NeutralDeltaInfo] = {}

        for delta_key, delta_info in neutral_delta_data.items():
            delta_type = NeutralDelta(delta_key)
            self._delta_to_data[delta_type] = delta_info
            self._formula_to_data[delta_info.formula.lower()] = delta_info
            self._name_to_data[delta_info.name.lower()] = delta_info

    def query_delta(self, delta: NeutralDelta) -> NeutralDeltaInfo | None:
        """Query by NeutralDelta enum"""
        return self._delta_to_data.get(delta)

    def query_formula(self, formula: str) -> NeutralDeltaInfo | None:
        """Query by chemical formula (e.g., 'H2O', 'NH3')"""
        return self._formula_to_data.get(formula.lower())

    def query_name(self, name: str) -> NeutralDeltaInfo | None:
        """Query by neutral delta name (e.g., 'Water', 'Ammonia')"""
        return self._name_to_data.get(name.lower())

    def __getitem__(self, key: str | NeutralDelta) -> NeutralDeltaInfo:
        """Get neutral delta by formula, name, or NeutralDelta enum"""
        if isinstance(key, NeutralDelta):
            info = self.query_delta(key)
            if info is not None:
                return info
            raise KeyError(f"Neutral delta type '{key}' not found.")

        # Try formula first
        info = self.query_formula(key)
        if info is not None:
            return info

        # Then try name
        info = self.query_name(key)
        if info is not None:
            return info

        raise KeyError(f"Neutral delta '{key}' not found by formula or name.")

    def __contains__(self, key: str | NeutralDelta) -> bool:
        """Check if neutral delta exists"""
        try:
            self[key]
            return True
        except KeyError:
            return False

    def __iter__(self) -> Iterable[NeutralDeltaInfo]:
        """Iterate over all neutral delta infos"""
        return iter(self._neutral_delta_data.values())

    def __len__(self) -> int:
        """Get the number of neutral deltas"""
        return len(self._neutral_delta_data)

    def __repr__(self) -> str:
        return f"NeutralDeltaLookup({len(self)} neutral deltas)"


NEUTRAL_DELTA_LOOKUP = NeutralDeltaLookup(NEUTRAL_DELTA_DICT)
