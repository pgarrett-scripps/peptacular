from typing import Iterable
from .data import PROTEASES_DICT, Proteases
from .dclass import ProteaseInfo


class ProteaseLookup:
    def __init__(self, data: dict[Proteases, ProteaseInfo]) -> None:
        self.name_to_info: dict[str, ProteaseInfo] = {
            info.name: info for info in data.values()
        }

        # make keys lowercase for case-insensitive lookup
        self.id_to_info = {k.lower(): v for k, v in data.items()}
        self.name_to_info = {k.lower(): v for k, v in self.name_to_info.items()}

    def query_id(self, protease_id: str) -> ProteaseInfo | None:
        """Query by protease ID (e.g., 'trypsin', 'arg-c')"""
        return self.id_to_info.get(protease_id.lower())

    def query_name(self, name: str) -> ProteaseInfo | None:
        """Query by protease name (e.g., 'Trypsin', 'Arg-C')"""
        return self.name_to_info.get(name.lower())

    def __getitem__(self, key: str) -> ProteaseInfo:
        """Get protease by ID or name"""
        # Try name first (more specific)
        info = self.query_name(key)
        if info is not None:
            return info

        # Then try ID
        info = self.query_id(key)
        if info is not None:
            return info

        raise KeyError(f"Protease '{key}' not found by name or ID.")

    def __contains__(self, key: str) -> bool:
        """Check if protease exists"""
        try:
            self[key]
            return True
        except KeyError:
            return False

    def get(self, key: str) -> ProteaseInfo | None:
        """Get protease or None if not found"""
        try:
            return self[key]
        except KeyError:
            return None

    def __iter__(self) -> Iterable[ProteaseInfo]:
        """Iterator over all ProteaseInfo entries in the lookup."""
        return iter(self.id_to_info.values())


PROTEASE_LOOKUP = ProteaseLookup(PROTEASES_DICT)
