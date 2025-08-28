from __future__ import annotations
import copy
from collections import UserDict
from collections.abc import Iterable, Mapping
from typing import Any

from ...mod import setup_mod
from .modlist import ModList, setup_mod_list, MODLIST_DATATYPE

MODDICT_VALUE_TYPES = ModList | MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE]


class ModDict(UserDict[int, ModList]):
    """
    A dictionary mapping non-negative integer positions to ModList instances.

    - Keys must be non-negative integers (>= 0).
    - Values are automatically converted to ModList instances.
    - Empty ModLists are automatically removed.
    - Supports flexible input types for values.

    Examples:
        md = ModDict()
        md[1] = "Oxidation"                    # stores ModList(["Oxidation"])
        md[2] = [15.994915, "Phospho"]         # stores ModList([15.994915, "Phospho"])
        md[3] = ModList([Mod("X", 2)])         # stores as-is

        # Bulk initialization
        md = ModDict({1: "phospho", 2: ["acetyl", 42.0]})
    """

    def __init__(
        self, initial: Mapping[int, MODDICT_VALUE_TYPES] | None = None
    ) -> None:
        super().__init__()
        if initial is not None:
            self.update(initial)

    def _validate_key(self, key: Any) -> None:
        """Validate that key is a non-negative integer"""
        if not isinstance(key, int):
            raise TypeError(f"ModDict keys must be integers, got {type(key).__name__}")
        if key < 0:
            raise ValueError("ModDict keys must be non-negative integers (>= 0)")

    def _normalize_value(self, value: MODDICT_VALUE_TYPES) -> ModList:
        """Convert value to ModList instance"""
        if isinstance(value, ModList):

            # ensure that modlsit has same attributes
            if not value.allow_dups or value.stackable:
                return setup_mod_list(
                    list(value), allow_dups=True, stackable=False
                )

            return value
        return setup_mod_list(value, allow_dups=True, stackable=False)

    def _cleanup_empty(self) -> None:
        """Remove any keys with empty ModLists"""
        empty_keys = [k for k, v in self.data.items() if not v]
        for k in empty_keys:
            del self.data[k]

    def __getitem__(self, key: int) -> ModList:
        """Get ModList for key, with validation"""
        self._validate_key(key)
        try:
            return self.data[key]
        except KeyError:
            raise KeyError(f"Key {key} not found in ModDict")

    def __setitem__(self, key: int, value: MODDICT_VALUE_TYPES | None) -> None:
        """Set value for key, converting to ModList and handling empty values"""
        self._validate_key(key)

        if value is None:
            # If value is None, remove the key if it exists
            if key in self.data:
                del self.data[key]
            return

        # Convert value to ModList
        modlist = self._normalize_value(value)

        # Only store non-empty ModLists
        if modlist:
            self.data[key] = modlist
        elif key in self.data:
            # Remove key if it exists and we're setting it to empty
            del self.data[key]

    def __delitem__(self, key: int) -> None:
        """Delete key with validation"""
        self._validate_key(key)
        del self.data[key]

    def __contains__(self, key: object) -> bool:
        """Check if key exists with validation"""
        if not isinstance(key, int) or key < 0:
            return False
        return key in self.data

    def get(self, key: object, default: Any = None) -> Any:
        """Get ModList for key with optional default"""
        # Use a broad signature compatible with UserDict.get while preserving ModDict behavior.
        if not isinstance(key, int) or key < 0:  # type: ignore[arg-type]
            if default is None:
                return None
            return self._normalize_value(default)

        if key in self.data:
            return self.data[key]

        if default is None:
            return None
        return self._normalize_value(default)

    def pop(self, key: int, default: Any = None) -> Any:
        """Pop and return ModList for key (or normalized default).

        Uses a flexible signature compatible with MutableMapping.pop: when a default
        is provided it will be normalized to a ModList before returning.
        """
        # Validate key as a non-negative integer
        self._validate_key(key)

        try:
            return self.data.pop(key)
        except KeyError:
            if default is None:
                raise
            return self._normalize_value(default)

    def setdefault(
        self, key: int, default: MODDICT_VALUE_TYPES | None = None
    ) -> ModList:
        """Get ModList for key, setting default if key doesn't exist"""
        self._validate_key(key)

        if key in self.data:
            return self.data[key]

        # Set and return default value
        default_modlist = (
            self._normalize_value(default) if default is not None else ModList()
        )

        if default_modlist:  # Only set if non-empty
            self.data[key] = default_modlist
            return self.data[key]
        else:
            return default_modlist

    def update(
        self, other: object | None = None, **kwargs: MODDICT_VALUE_TYPES
    ) -> None:
        """
        Update ModDict with another mapping or iterable of (key, value) pairs.
        Values are normalized to ModList via __setitem__.

        The parameter 'other' is intentionally typed as 'object | None' to avoid
        narrowing the signature from the base MutableMapping.update, while the
        runtime behavior remains the same.
        """
        if other is not None:
            if hasattr(other, "keys"):
                # Mapping-like object
                for key in other:  # type: ignore[call-arg]
                    self[key] = other[key]  # type: ignore[index]
            else:
                # Iterable of (key, value) pairs
                for key, value in other:  # type: ignore[misc]
                    self[key] = value

        # Handle keyword arguments
        for key, value in kwargs.items():
            try:
                int_key = int(key)
                self[int_key] = value
            except ValueError:
                raise TypeError(
                    f"Keyword argument key '{key}' cannot be converted to int"
                )

    def __eq__(self, other: Any) -> bool:
        """Equality comparison for ModDict objects"""
        if not isinstance(other, ModDict):
            try:
                other = setup_mod_dict(other)
            except (TypeError, ValueError):
                return False

        if set(self.keys()) != set(other.keys()):
            return False

        for key in self.keys():
            if self[key] != other[key]:
                return False
        return True

    def __repr__(self) -> str:
        return f"ModDict({dict(self.data)})"

    def copy(self, deep: bool = True) -> ModDict:
        """Create a copy of the ModDict"""
        if deep:
            return copy.deepcopy(self)
        else:
            result = ModDict()
            for key, value in self.data.items():
                result.data[key] = value.copy(deep=False)
            return result

    def merge(self, other: Mapping[int, MODDICT_VALUE_TYPES]) -> None:
        """
        Merge another ModDict or compatible mapping into this ModDict.
        Values at the same position will be combined into the same ModList.

        Args:
            other: Another ModDict, mapping, or iterable of (key, value) pairs

        Examples:
            md1 = ModDict({1: ["Oxidation"], 2: [15.994915]})
            md2 = ModDict({1: ["Phospho"], 3: ["Acetyl"]})
            md1.merge(md2)
            # md1 now contains: {1: ["Oxidation", "Phospho"], 2: [15.994915], 3: ["Acetyl"]}
        """

        # Convert to standardized format
        if isinstance(other, ModDict):
            items = other.items()
        elif isinstance(other, Mapping):  # type: ignore
            # Mapping-like object
            items = [
                (k, setup_mod_list(v, allow_dups=True, stackable=False))
                for k, v in other.items()
            ]
        else:
            raise TypeError(f"Cannot merge ModDict with {type(other)}")

        for key, value in items:
            self._validate_key(key)

            # Convert value to ModList
            new_modlist = self._normalize_value(value)

            if key in self.data:
                # Extend existing ModList with new values
                self.data[key] += new_modlist
            else:
                # Set new key with the ModList
                self[key] = new_modlist

        # Clean up any empty ModLists that might have been created
        self._cleanup_empty()

    def extend_at_key(self, key: int, value: MODDICT_VALUE_TYPES) -> None:
        """Extend ModList at key with additional values"""
        self._validate_key(key)

        if key in self.data:
            self.data[key] += self._normalize_value(value)
        else:
            self[key] = value

    def append_at_key(self, key: int, value: MODLIST_DATATYPE) -> None:
        """Append single value to ModList at key"""
        self._validate_key(key)

        if key in self.data:
            self.data[key].append(value)
        else:
            self[key] = [value]

    def remove_from_key(self, key: int, value: MODLIST_DATATYPE) -> None:
        """Remove value from ModList at key"""
        self._validate_key(key)

        if key not in self.data:
            raise KeyError(f"Key {key} not found in ModDict")

        self.data[key].remove(value)

        # Clean up if ModList becomes empty
        if not self.data[key]:
            del self.data[key]

    def discard_from_key(self, key: int, value: MODLIST_DATATYPE) -> None:
        """Remove value from ModList at key without raising error if not found"""
        try:
            self.remove_from_key(key, value)
        except (KeyError, ValueError):
            pass

    def keys_with_value(self, value: MODLIST_DATATYPE) -> list[int]:
        """Return list of keys that contain the specified value"""
        mod = setup_mod(value)
        return [key for key, modlist in self.items() if mod in modlist]

    def positions_with_mods(self) -> list[int]:
        """Return sorted list of positions that have modifications"""
        return sorted(self.keys())

    def is_empty(self) -> bool:
        """Check if ModDict has no modifications"""
        return len(self.data) == 0

    @property
    def has_mods(self) -> bool:
        """Check if ModDict has any modifications"""
        return len(self.data) > 0


ACCEPTED_MODDICT_INPUT_TYPES = Mapping[int, MODDICT_VALUE_TYPES] | ModDict | None


def setup_mod_dict(initial: ACCEPTED_MODDICT_INPUT_TYPES = None) -> ModDict:
    """
    Initialize a ModDict from various input types.

    Args:
        initial: Initial data - can be a ModDict, mapping, iterable of (key, value) pairs, or None

    Returns:
        A ModDict instance

    Examples:
        setup_mod_dict({1: "phospho", 2: ["acetyl", 42.0]})
        setup_mod_dict([(1, "phospho"), (2, ["acetyl", 42.0])])
        setup_mod_dict(existing_moddict)
    """

    if isinstance(initial, ModDict):
        return initial

    return ModDict(initial)
