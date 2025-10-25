from __future__ import annotations

import warnings
from collections import UserList
from collections.abc import Iterable
from typing import Any, Self, SupportsIndex

from ...mod import MOD_VALUE_TYPES, Mod, setup_mod

MODLIST_DATATYPE = MOD_VALUE_TYPES | Mod


class ModList(UserList[Mod]):
    """
    A list of modifications that automatically merges duplicates by aggregating multipliers.
    """

    def __init__(
        self,
        data: Iterable[MODLIST_DATATYPE] | None = None,
        allow_dups: bool = True,
        stackable: bool = True,
    ) -> None:
        super().__init__()
        self.allow_dups = allow_dups
        self.stackable = stackable

        if allow_dups is False and stackable is True:
            warnings.warn(
                "ModList is non-duplicable but stackable. Think long and hard about this one...",
                UserWarning,
                stacklevel=2,
            )

        if data is not None:
            self.extend(data)

    def _normalize_input(self, item: MODLIST_DATATYPE) -> Mod:
        """Convert input to Mod instance - optimized to check Mod first"""
        if isinstance(item, Mod):
            return item
        return setup_mod(item)

    def _find_same_mod_index(self, mod: Mod) -> int | None:
        """Find index of existing mod with same value"""
        mod_val = mod.val
        for i, m in enumerate(self.data):
            if m.val == mod_val:
                return i
        return None

    def _merge_or_append(self, item: Mod) -> None:
        """Merge with existing mod or append new one"""
        idx = self._find_same_mod_index(item)
        if idx is not None:
            if not self.allow_dups and item.mult > 0:
                raise ValueError(
                    f"Cannot add modification {item} to non-stackable ModList"
                )

            self.data[idx] = Mod(self.data[idx].val, self.data[idx].mult + item.mult)
            if self.data[idx].mult == 0:
                del self.data[idx]
            elif self.data[idx].mult < 0:
                raise ValueError(
                    f"Modification multiplier cannot be negative: {self.data[idx].mult}"
                )
        else:
            if item.mult != 0:
                self.data.append(item)

    def append(self, item: MODLIST_DATATYPE) -> None:
        """Add a modification, merging with existing mod if same value"""
        # Fast path: if already a Mod, skip normalization
        if isinstance(item, Mod):
            self._merge_or_append(item)
        else:
            mod = setup_mod(item)
            self._merge_or_append(mod)

    def extend(self, other: Iterable[MODLIST_DATATYPE]) -> None:
        """Extend with modifications, merging duplicates"""
        # Fast path for ModList
        if isinstance(other, ModList):
            for item in other.data:
                self._merge_or_append(item)
        elif isinstance(other, Iterable):  # type: ignore
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

        # self.data[idx].mult -= mod.mult (now frozen)
        self.data[idx] = Mod(self.data[idx].val, self.data[idx].mult - mod.mult)

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
        result = ModList()
        result.extend(self.data)

        if isinstance(other, ModList):
            result.extend(other.data)
        elif isinstance(other, Iterable):  # type: ignore
            result.extend(other)
        else:
            raise TypeError(f"Cannot add ModList with {type(other)}")

        return result

    def __iadd__(self, other: Iterable[MODLIST_DATATYPE] | Self) -> ModList:
        """In-place addition"""
        if isinstance(other, ModList):
            self.extend(other.data)
        elif isinstance(other, Iterable):  # type: ignore
            self.extend(other)
        else:
            raise TypeError(f"Cannot add ModList with {type(other)}")

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

        for mod in self.data:
            other_idx = other._find_same_mod_index(mod)
            if other_idx is None or other.data[other_idx].mult != mod.mult:
                return False
        return True

    # ModList - optimized shallow copy
    def copy(self, deep: bool = True) -> ModList:
        """Create a copy of the ModList - optimized to bypass __init__"""
        result = object.__new__(ModList)
        result.allow_dups = self.allow_dups
        result.stackable = self.stackable

        # Micro-optimization: avoid .copy() overhead for empty lists
        if self.data:
            result.data = self.data.copy()
        else:
            result.data = []

        return result

    def flatten(self, sort: bool = False) -> tuple[MOD_VALUE_TYPES, ...]:
        """Return tuple of all mod values, expanding multipliers"""
        values: list[MOD_VALUE_TYPES] = []
        for mod in self.data:
            values.extend([mod.val] * mod.mult)

        if sort:
            values.sort(key=str)

        return tuple(values)

    def serialize(
        self,
        brackets: str,
        sort: bool = False,
        include_plus: bool = False,
        precision: int | None = None,
    ) -> str:
        """Serialize the ModList into a string representation."""
        elems: list[str] = []
        for mod in self.data:
            if self.stackable:
                elems.append(
                    mod.serialize(
                        brackets, include_plus=include_plus, precision=precision
                    )
                )
            else:
                for _ in range(mod.mult):
                    elems.append(
                        Mod(mod.val, 1).serialize(
                            brackets, include_plus=include_plus, precision=precision
                        )
                    )

        if sort:
            elems.sort()

        return "".join(elems)

    @property
    def has_mods(self) -> bool:
        return len(self.data) > 0

    @property
    def is_empty(self) -> bool:
        return len(self.data) == 0

    def __repr__(self) -> str:
        return f"ModList({self.data})"


ACCEPTED_MODLIST_INPUT_TYPES = (
    MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE] | ModList | None
)


def setup_mod_list(
    mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE] | ModList | None,
    allow_dups: bool = False,
    stackable: bool = False,
) -> ModList:
    """Helper function to set up a ModList from various input types."""
    if isinstance(mods, ModList):
        return mods

    mod_list: ModList = ModList(allow_dups=allow_dups, stackable=stackable)

    if mods is None:
        pass
    elif isinstance(mods, (str, int, float, Mod)):
        mod_list.append(setup_mod(mods))
    elif isinstance(mods, Iterable):  # type: ignore
        for mod in mods:
            mod_list.append(setup_mod(mod))
    else:
        raise TypeError(f"Invalid type for isotope_mods: {type(mods)}")
    return mod_list


def populate_mod_list(
    modlist: ModList,
    mods: MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE] | ModList | None,
) -> None:
    """Helper function to populate an existing ModList from various input types."""

    if mods is None:
        return

    if isinstance(mods, ModList):
        modlist.extend(mods.data)
    elif isinstance(mods, (str, int, float, Mod)):
        modlist.append(setup_mod(mods))
    elif isinstance(mods, Iterable):  # type: ignore
        for mod in mods:
            modlist.append(setup_mod(mod))
    else:
        raise TypeError(f"Invalid type for isotope_mods: {type(mods)}")
