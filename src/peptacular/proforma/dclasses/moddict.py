from __future__ import annotations
from collections import UserDict
from collections.abc import Iterable, Mapping
from typing import Any

from ...mod import setup_mod
from .modlist import ModList, setup_mod_list, MODLIST_DATATYPE

MODDICT_VALUE_TYPES = ModList | MODLIST_DATATYPE | Iterable[MODLIST_DATATYPE]


class ModDict(UserDict[int, ModList]):
    """Dictionary mapping positions to ModList instances."""

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
        """Convert value to ModList instance - optimized"""
        # Fast path: already a ModList with correct settings
        if isinstance(value, ModList):
            if value.allow_dups and not value.stackable:
                return value
            # Need to convert
            return setup_mod_list(value.data, allow_dups=True, stackable=False)

        # Convert other types
        return setup_mod_list(value, allow_dups=True, stackable=False)

    def _cleanup_empty(self) -> None:
        """Remove any keys with empty ModLists"""
        empty_keys = [k for k, v in self.data.items() if not v]
        for k in empty_keys:
            del self.data[k]

    def __getitem__(self, key: int) -> ModList:
        """Get ModList for key, with validation"""
        self._validate_key(key)
        return self.data[key]

    def __setitem__(self, key: int, value: MODDICT_VALUE_TYPES | None) -> None:
        """Set value for key, converting to ModList and handling empty values"""
        self._validate_key(key)

        if value is None:
            if key in self.data:
                del self.data[key]
            return

        modlist = self._normalize_value(value)

        if modlist:
            self.data[key] = modlist
        elif key in self.data:
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
        if not isinstance(key, int) or key < 0:
            return None if default is None else self._normalize_value(default)

        if key in self.data:
            return self.data[key]

        return None if default is None else self._normalize_value(default)

    def pop(self, key: int, default: Any = None) -> Any:
        """Pop and return ModList for key"""
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

        default_modlist = (
            self._normalize_value(default) if default is not None else ModList()
        )

        if default_modlist:
            self.data[key] = default_modlist
            return self.data[key]
        else:
            return default_modlist

    def update(
        self, other: object | None = None, **kwargs: MODDICT_VALUE_TYPES
    ) -> None:
        """Update ModDict with another mapping or iterable"""
        if other is not None:
            if hasattr(other, "keys"):
                for key in other:  # type: ignore
                    self[key] = other[key]  # type: ignore
            else:
                for key, value in other:  # type: ignore
                    self[key] = value

        for key, value in kwargs.items():
            try:
                self[int(key)] = value
            except ValueError:
                raise TypeError(
                    f"Keyword argument key '{key}' cannot be converted to int"
                )

    def __eq__(self, other: Any) -> bool:
        """Equality comparison"""
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
        """Create a copy of the ModDict - optimized to bypass __init__"""
        result = object.__new__(ModDict)

        # Fast path: empty dict
        if not self.data:
            result.data = {}
            return result

        # Copy with pre-allocated dict
        result.data = {key: modlist.copy() for key, modlist in self.data.items()}

        return result

    def merge(self, other: Mapping[int, MODDICT_VALUE_TYPES]) -> None:
        """Merge another ModDict into this one"""
        if isinstance(other, ModDict):
            items = other.items()
        elif isinstance(other, Mapping):  # type: ignore
            items = [
                (k, setup_mod_list(v, allow_dups=True, stackable=False))
                for k, v in other.items()
            ]
        else:
            raise TypeError(f"Cannot merge ModDict with {type(other)}")

        for key, value in items:
            self._validate_key(key)
            new_modlist = self._normalize_value(value)

            if key in self.data:
                self.data[key] += new_modlist
            else:
                self[key] = new_modlist

        self._cleanup_empty()

    def extend_at_key(self, key: int, value: MODDICT_VALUE_TYPES) -> None:
        """Extend ModList at key"""
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

        if not self.data[key]:
            del self.data[key]

    def discard_from_key(self, key: int, value: MODLIST_DATATYPE) -> None:
        """Remove value from ModList at key without raising error"""
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

    def clean_empty_lists(self):
        for key in list(self.data.keys()):
            if not self.data[key]:
                del self.data[key]


ACCEPTED_MODDICT_INPUT_TYPES = Mapping[int, MODDICT_VALUE_TYPES] | ModDict | None


def setup_mod_dict(initial: ACCEPTED_MODDICT_INPUT_TYPES = None) -> ModDict:
    """Initialize a ModDict from various input types."""
    if isinstance(initial, ModDict):
        return initial
    return ModDict(initial)
