from collections.abc import Iterable

from .data import ION_TYPE_DICT, IonType
from .dclass import FragmentIonInfo


class FragmentIonLookup:
    def __init__(self, fragment_ion_data: dict[IonType, FragmentIonInfo]) -> None:
        self._fragment_ion_data = fragment_ion_data

        self._iontype_to_data: dict[IonType, FragmentIonInfo] = {}
        self._id_to_data: dict[str, FragmentIonInfo] = {}
        self._name_to_data: dict[str, FragmentIonInfo] = {}

        for ion_id, ion_info in fragment_ion_data.items():
            ion_type = IonType(ion_id)
            self._iontype_to_data[ion_type] = ion_info
            self._id_to_data[ion_id.lower()] = ion_info
            self._name_to_data[ion_info.name.lower()] = ion_info

    def query_ion_type(self, ion_type: IonType) -> FragmentIonInfo | None:
        """Query by IonType enum"""
        return self._iontype_to_data.get(ion_type)

    def query_id(self, ion_id: str) -> FragmentIonInfo | None:
        """Query by fragment ion ID (e.g., 'a', 'b', 'y', 'z')"""
        return self._id_to_data.get(ion_id.lower())

    def query_name(self, name: str) -> FragmentIonInfo | None:
        """Query by fragment ion name (e.g., 'A-Ion', 'Y-Ion')"""
        return self._name_to_data.get(name.lower())

    def __getitem__(self, key: str | IonType) -> FragmentIonInfo:
        """Get fragment ion by ID, name, or IonType"""
        if isinstance(key, IonType):
            info = self.query_ion_type(key)
            if info is not None:
                return info
            raise KeyError(f"Fragment ion type '{key}' not found.")

        # Try ID first
        info = self.query_id(key)
        if info is not None:
            return info

        # Then try name
        info = self.query_name(key)
        if info is not None:
            return info

        raise KeyError(f"Fragment ion '{key}' not found by ID or name.")

    def __contains__(self, key: str | IonType) -> bool:
        """Check if fragment ion exists"""
        try:
            self[key]
            return True
        except KeyError:
            return False

    def get(self, key: str | IonType) -> FragmentIonInfo | None:
        """Get fragment ion or None if not found"""
        try:
            return self[key]
        except KeyError:
            return None

    def __iter__(self) -> Iterable[FragmentIonInfo]:
        """Iterator over all FragmentIonInfo entries in the lookup."""
        return iter(self._fragment_ion_data.values())


FRAGMENT_ION_LOOKUP = FragmentIonLookup(ION_TYPE_DICT)
