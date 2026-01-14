from collections.abc import Iterable
from functools import cached_property
from random import choice

from ..dclass import UnimodInfo, filter_infos
from .data import UNIMOD_MODIFICATIONS


class UnimodLookup:
    def __init__(self, data: dict[str, UnimodInfo]) -> None:
        self.id_to_info: dict[str, UnimodInfo] = data
        self.name_to_info: dict[str, UnimodInfo] = {
            info.name: info for info in data.values()
        }

        # make keys lowercase for case-insensitive lookup
        self.id_to_info = {k.lower(): v for k, v in self.id_to_info.items()}
        self.name_to_info = {k.lower(): v for k, v in self.name_to_info.items()}

    def query_id(self, mod_id: str | int) -> UnimodInfo | None:
        if isinstance(mod_id, int):
            mod_id = str(mod_id)

        # Strip UNIMOD: prefix if present (case-insensitive)
        if ":" in mod_id:
            prefix, value = mod_id.split(":", 1)
            if prefix.lower() == "unimod":
                mod_id = value

        return self.id_to_info.get(mod_id.lower())

    def query_name(self, name: str) -> UnimodInfo | None:
        # Strip UNIMOD: prefix if present (case-insensitive)
        if ":" in name:
            prefix, value = name.split(":", 1)
            if prefix.lower() == "u":
                name = value

        return self.name_to_info.get(name.lower())

    def query_mass(
        self, mass: float, tolerance: float = 0.01, monoisotopic: bool = True
    ) -> UnimodInfo | None:
        """Query Unimod modification by mass within a given tolerance."""
        matches: list[UnimodInfo] = []
        for info in self.id_to_info.values():
            mod_mass = info.monoisotopic_mass if monoisotopic else info.average_mass
            if mod_mass is not None and abs(mod_mass - mass) <= tolerance:
                matches.append(info)
        if len(matches) == 1:
            return matches[0]
        elif len(matches) > 1:
            # if all have the same composition, return the first one
            compositions = {
                tuple(
                    sorted(
                        m.dict_composition.items()
                        if m.dict_composition is not None
                        else []
                    )
                )
                for m in matches
            }
            if len(compositions) == 1:
                return matches[0]
            raise ValueError(
                f"Multiple UNIMOD modifications found for mass {mass} within tolerance {tolerance}: {[(m.id, m.monoisotopic_mass, m.formula) for m in matches]}"
            )
        return None

    def __getitem__(self, key: str) -> UnimodInfo:
        info = self.query_name(key)
        if info is not None:
            return info

        info = self.query_id(key)
        if info is not None:
            return info

        raise KeyError(f"Unimod modification '{key}' not found by name or ID.")

    def __contains__(self, key: str) -> bool:
        try:
            self[key]
            return True
        except KeyError:
            return False

    def get(self, key: str) -> UnimodInfo | None:
        try:
            return self[key]
        except KeyError:
            return None

    def __iter__(self) -> Iterable[UnimodInfo]:
        """Iterator over all UnimodInfo entries in the lookup."""
        return iter(self.name_to_info.values())

    def values(self) -> Iterable[UnimodInfo]:
        """Get all UnimodInfo entries in the lookup."""
        return self.name_to_info.values()

    def keys(self) -> Iterable[str]:
        """Get all keys (names) in the lookup."""
        return self.name_to_info.keys()

    @cached_property
    def _all_infos_tuple(self) -> tuple[UnimodInfo, ...]:
        """Cached tuple of all UnimodInfo entries."""
        return tuple(self.name_to_info.values())

    @cached_property
    def _infos_with_mass_tuple(self) -> tuple[UnimodInfo, ...]:
        """Cached tuple of UnimodInfo entries with monoisotopic mass."""
        return tuple(
            filter_infos(list(self.name_to_info.values()), has_monoisotopic_mass=True)
        )

    @cached_property
    def _infos_with_composition_tuple(self) -> tuple[UnimodInfo, ...]:
        """Cached tuple of UnimodInfo entries with composition."""
        return tuple(
            filter_infos(list(self.name_to_info.values()), has_composition=True)
        )

    @cached_property
    def _infos_with_mass_and_composition_tuple(self) -> tuple[UnimodInfo, ...]:
        """Cached tuple of UnimodInfo entries with both mass and composition."""
        return tuple(
            filter_infos(
                list(self.name_to_info.values()),
                has_monoisotopic_mass=True,
                has_composition=True,
            )
        )

    def choice(
        self, require_monoisotopic_mass: bool = True, require_composition: bool = True
    ) -> UnimodInfo:
        """Get a random UnimodInfo from the lookup."""
        if require_monoisotopic_mass and require_composition:
            valid_infos = self._infos_with_mass_and_composition_tuple
        elif require_monoisotopic_mass:
            valid_infos = self._infos_with_mass_tuple
        elif require_composition:
            valid_infos = self._infos_with_composition_tuple
        else:
            valid_infos = self._all_infos_tuple

        if not valid_infos:
            raise ValueError("No valid UnimodInfo entries found matching the criteria.")

        return choice(valid_infos)


UNIMOD_LOOKUP = UnimodLookup(UNIMOD_MODIFICATIONS)
