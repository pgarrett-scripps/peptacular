from __future__ import annotations
import copy
from collections import UserList
from collections.abc import Iterable
from typing import Any, Self, SupportsIndex
import warnings


from ...mod import Mod, MOD_VALUE_TYPES, setup_mod

# MOD_VALUE_TYPES = str | int | float
MODLIST_DATATYPE = MOD_VALUE_TYPES | Mod


class ModList(UserList[Mod]):
    """
    A list of modifications that automatically merges duplicates by aggregating multipliers.

    Accepts str, int, float, Mod instances, or iterables thereof.
    Internally stores everything as Mod instances with no duplicate values.
    """

    def __init__(self, data: Iterable[MODLIST_DATATYPE] | None = None) -> None:
        # Initialize empty UserList
        super().__init__()

        if data is not None:
            self.extend(data)

    def _normalize_input(self, item: MODLIST_DATATYPE) -> Mod:
        """Convert input to Mod instance"""
        return setup_mod(item)

    def _find_same_mod_index(self, mod: Mod) -> int | None:
        """Find index of existing mod with same value"""
        for i, m in enumerate(self.data):
            if m.val == mod.val:
                return i
        return None

    def _merge_or_append(self, item: Mod) -> None:
        """Merge with existing mod or append new one"""
        idx = self._find_same_mod_index(item)
        if idx is not None:
            self.data[idx].mult += item.mult
            # Remove if multiplier becomes zero
            if self.data[idx].mult == 0:
                del self.data[idx]
            elif self.data[idx].mult < 0:
                raise ValueError(
                    f"Modification multiplier cannot be negative: {self.data[idx].mult}"
                )
        else:
            if item.mult != 0:  # Only append if multiplier is non-zero
                self.data.append(item)

    # Override UserList methods to handle type conversion and merging

    def append(self, item: MODLIST_DATATYPE) -> None:
        """Add a modification, merging with existing mod if same value"""
        mod = self._normalize_input(item)
        self._merge_or_append(mod)

    def extend(self, other: Iterable[MODLIST_DATATYPE]) -> None:
        """Extend with modifications, merging duplicates"""
        if isinstance(other, Iterable):  # type: ignore
            for item in other:
                self.append(item)
        else:
            raise TypeError(f"Cannot extend ModList with {type(other)}")

    def insert(self, i: int, item: MODLIST_DATATYPE) -> None:
        """Insert behaves like append (merges duplicates, ignores index)"""
        warnings.warn(
            "ModList.insert() does not support index and behaves like append().",
            UserWarning,
            stacklevel=2,
        )
        self.append(item)

    def remove(self, item: MODLIST_DATATYPE) -> None:
        """Remove specified amount of modification (decrement multiplier)"""
        mod = self._normalize_input(item)
        idx = self._find_same_mod_index(mod)

        if idx is None:
            raise ValueError(f"{item} is not in list")

        self.data[idx].mult -= mod.mult
        if self.data[idx].mult == 0:
            del self.data[idx]
        elif self.data[idx].mult < 0:
            raise ValueError(
                f"Modification multiplier cannot be negative: {self.data[idx].mult}"
            )

    def discard(self, item: MODLIST_DATATYPE) -> None:
        """Remove modification without raising error if not found"""
        try:
            self.remove(item)
        except ValueError:
            pass

    def count(self, item: MODLIST_DATATYPE) -> int:
        """Count total multiplier for a modification"""
        mod = self._normalize_input(item)
        idx = self._find_same_mod_index(mod)
        return self.data[idx].mult if idx is not None else 0

    def index(
        self,
        item: MODLIST_DATATYPE,
        start: SupportsIndex = 0,
        stop: SupportsIndex | None = None,
    ) -> int:
        """Find index of modification by value"""
        mod = self._normalize_input(item)

        length = len(self.data)
        start_idx = start if isinstance(start, int) else 0
        stop_idx = stop if isinstance(stop, int) else length

        for i in range(start_idx, min(stop_idx, length)):
            if self.data[i].val == mod.val:
                return i
        raise ValueError(f"{item} is not in list")

    def __contains__(self, item: object) -> bool:
        """Check if modification value exists in list"""
        try:
            if isinstance(item, (str, int, float, Mod)):
                mod = self._normalize_input(item)
                return self._find_same_mod_index(mod) is not None
        except (TypeError, ValueError):
            pass
        return False

    def __add__(self, other: Iterable[MODLIST_DATATYPE] | Self) -> ModList:
        """Return new ModList combining both lists"""
        if isinstance(other, Iterable):  # type: ignore
            other_modlist: Self = self.__class__(other)
        elif isinstance(other, ModList):  # type: ignore
            other_modlist: Self = other
        else:
            raise TypeError(f"Cannot add ModList with {type(other)}")

        result = ModList()
        result.extend(self.data)
        result.extend(other_modlist.data)
        return result

    def __iadd__(self, other: Iterable[MODLIST_DATATYPE] | Self) -> ModList:
        """In-place addition"""
        if isinstance(other, Iterable):  # type: ignore
            other_modlist: Self = self.__class__(other)
        elif isinstance(other, ModList):  # type: ignore
            other_modlist: Self = other
        else:
            raise TypeError(f"Cannot add ModList with {type(other)}")

        self.extend(other_modlist.data)
        return self

    def __eq__(self, other: Any) -> bool:
        """Equality comparison as multisets (order-independent)"""
        if not isinstance(other, ModList):
            try:
                other = ModList(other)
            except (TypeError, ValueError):
                return False

        if len(self.data) != len(other.data):
            return False

        # Check each mod exists in other with same multiplier
        for mod in self.data:
            other_idx = other._find_same_mod_index(mod)
            if other_idx is None or other.data[other_idx].mult != mod.mult:
                return False
        return True

    def copy(self, deep: bool = True) -> ModList:
        """Create a copy of the ModList"""
        if deep:
            return copy.deepcopy(self)
        else:
            result = ModList()
            result.data = self.data.copy()
            return result

    def flatten(self, sort: bool = False) -> tuple[MOD_VALUE_TYPES, ...]:
        """
        Return tuple of all mod values, expanding multipliers

        Args:
            sort: Whether to sort the resulting values

        Returns:
            Tuple of mod values with multipliers expanded
        """
        values: list[MOD_VALUE_TYPES] = []
        for mod in self.data:
            values.extend([mod.val] * mod.mult)

        if sort:
            values.sort(key=str)

        return tuple(values)

    @property
    def has_mods(self) -> bool:
        """Check if list has any modifications"""
        return len(self.data) > 0

    @property
    def is_empty(self) -> bool:
        """Check if list is empty"""
        return len(self.data) == 0

    def __repr__(self) -> str:
        return f"ModList({self.data})"


ACCEPTED_MODLIST_INPUT_TYPES = (
    MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE] | ModList | None
)


def setup_mod_list(
    mods: ACCEPTED_MODLIST_INPUT_TYPES,
) -> ModList:
    """
    Helper function to set up a ModList from various input types.
    """

    if isinstance(mods, ModList):
        return mods

    mod_list: ModList = ModList()

    if mods is None:
        pass
    elif isinstance(mods, (str, int, float, Mod)):
        mod_list.append(setup_mod(mods))
    elif isinstance(mods, Iterable):  # type: ignore
        for mod in mods:
            mod_list.append(setup_mod(mod))
    else:
        raise TypeError(
            f"Invalid type for isotope_mods: {type(mods)}. Must be a string, int, float, Mod, a (value,int) tuple, or iterable of these."
        )
    return mod_list
