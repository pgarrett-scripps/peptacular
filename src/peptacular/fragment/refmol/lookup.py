from collections.abc import Iterable

from .data import REFMOL_DICT, RefMolID
from .dclass import RefMolInfo


class RefMolLookup:
    def __init__(self, refmol_data: dict[RefMolID, RefMolInfo]) -> None:
        self._refmol_data = refmol_data

        self._refmolid_to_data: dict[RefMolID, RefMolInfo] = {}
        self._name_to_data: dict[str, RefMolInfo] = {}
        self._label_type_to_data: dict[str, list[RefMolInfo]] = {}
        self._molecule_type_to_data: dict[str, list[RefMolInfo]] = {}

        for refmol_id, refmol_info in refmol_data.items():
            self._refmolid_to_data[refmol_id] = refmol_info
            self._name_to_data[refmol_info.name.lower()] = refmol_info

            # Group by label type
            if refmol_info.label_type:
                self._label_type_to_data.setdefault(
                    refmol_info.label_type.lower(), []
                ).append(refmol_info)

            # Group by molecule type
            if refmol_info.molecule_type:
                self._molecule_type_to_data.setdefault(
                    refmol_info.molecule_type.lower(), []
                ).append(refmol_info)

    def query_id(self, refmol_id: RefMolID) -> RefMolInfo | None:
        """Query by RefMolID enum"""
        return self._refmolid_to_data.get(refmol_id)

    def query_name(self, name: str) -> RefMolInfo | None:
        """Query by reference molecule name (e.g., 'TMT126', 'sidechain_A')"""
        return self._name_to_data.get(name.lower())

    def query_label_type(self, label_type: str) -> list[RefMolInfo]:
        """Query all molecules by label type (e.g., 'TMT', 'iTRAQ')"""
        return self._label_type_to_data.get(label_type.lower(), [])

    def query_molecule_type(self, molecule_type: str) -> list[RefMolInfo]:
        """Query all molecules by molecule type (e.g., 'reporter', 'sidechain', 'nucleobase')"""
        return self._molecule_type_to_data.get(molecule_type.lower(), [])

    def __getitem__(self, key: str | RefMolID) -> RefMolInfo:
        """Get reference molecule by ID or name"""
        if isinstance(key, RefMolID):
            info = self.query_id(key)
            if info is not None:
                return info
            raise KeyError(f"Reference molecule ID '{key}' not found.")

        # Try name
        info = self.query_name(key)
        if info is not None:
            return info

        raise KeyError(f"Reference molecule '{key}' not found by name.")

    def __contains__(self, key: str | RefMolID) -> bool:
        """Check if reference molecule exists"""
        try:
            self[key]
            return True
        except KeyError:
            return False

    def get(self, key: str | RefMolID) -> RefMolInfo | None:
        """Get reference molecule or None if not found"""
        try:
            return self[key]
        except KeyError:
            return None

    def __iter__(self) -> Iterable[RefMolInfo]:
        """Iterator over all RefMolInfo entries in the lookup."""
        return iter(self._refmol_data.values())


REFMOL_LOOKUP = RefMolLookup(REFMOL_DICT)
