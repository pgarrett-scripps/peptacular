from collections import UserList
from typing import *

from ..proforma_dataclasses import Mod, convert_to_mod, ModValue


class ModList(UserList):
    """
    Represents a list of modifications in ProForma format.
    """

    def __init__(
        self, mods: Optional[Union[Iterable[ModValue], ModValue]] = None
    ) -> None:
        super().__init__()

        if mods is None:
            return

        if isinstance(mods, (str, int, float, Mod)):
            mods = [mods]

        self.extend(mods)

    def _validate_mod(self, mod: Mod) -> None:
        """Validate that the item is a Mod instance."""
        if not isinstance(mod, Mod):
            raise TypeError(f"Expected Mod instance, got {type(mod)}")

    def append(self, mod: ModValue) -> None:
        """Add a modification to the list."""

        if not isinstance(mod, Mod):
            mod = convert_to_mod(mod)

        self._validate_mod(mod)
        super().append(mod)

    def extend(self, mods: Iterable[ModValue]) -> None:
        """Extend the list with multiple modifications."""
        for mod in mods:
            if not isinstance(mod, Mod):
                mod = convert_to_mod(mod)
            self._validate_mod(mod)
            super().append(mod)

    def insert(self, index: int, mod: ModValue) -> None:
        """Insert a modification at the specified index."""
        if not isinstance(mod, Mod):
            mod = convert_to_mod(mod)
        self._validate_mod(mod)
        super().insert(index, mod)

    def __setitem__(
        self, index: Union[int, slice], mod: Union[ModValue, Iterable[ModValue]]
    ) -> None:
        """Set a modification at a specific index or slice."""
        if isinstance(index, slice):
            # Handle slice assignment
            converted_mods = []
            for m in mod:
                if not isinstance(m, Mod):
                    m = convert_to_mod(m)
                self._validate_mod(m)
                converted_mods.append(m)
            super().__setitem__(index, converted_mods)
        else:
            # Handle single item assignment
            if not isinstance(mod, Mod):
                mod = convert_to_mod(mod)
            self._validate_mod(mod)
            super().__setitem__(index, mod)

    def remove(self, mod: ModValue) -> None:
        """Remove the first occurrence of a modification."""
        if not isinstance(mod, Mod):
            mod = convert_to_mod(mod)
        super().remove(mod)

    def count(self, mod: ModValue) -> int:
        """Count occurrences of a modification."""
        if not isinstance(mod, Mod):
            mod = convert_to_mod(mod)
        return super().count(mod)

    def index(self, mod: ModValue, start: int = 0, stop: int = None) -> int:
        """Find the index of a modification."""
        if not isinstance(mod, Mod):
            mod = convert_to_mod(mod)
        if stop is None:
            return super().index(mod, start)
        return super().index(mod, start, stop)

    def __contains__(self, mod: ModValue) -> bool:
        """Check if a modification is in the list."""
        if not isinstance(mod, Mod):
            try:
                mod = convert_to_mod(mod)
            except (TypeError, ValueError):
                return False
        return super().__contains__(mod)

    def __add__(self, other: Union["ModList", list]) -> "ModList":
        """Concatenate with another ProformaModList or list."""
        if isinstance(other, list):
            other = ModList(other)
        elif not isinstance(other, ModList):
            raise TypeError(
                f"unsupported operand type(s) for +: 'ModList' and '{type(other).__name__}'"
            )

        result = ModList(self.data)
        result.extend(other)
        return result

    def __iadd__(self, other: Union["ModList", list]) -> "ModList":
        """In-place addition."""
        if isinstance(other, list):
            self.extend(other)
        else:
            self.extend(other.data)
        return self

    def __repr__(self) -> str:
        return f"ModList({self.data})"
