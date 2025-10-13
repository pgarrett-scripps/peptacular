from __future__ import annotations
from functools import lru_cache
from typing import Any

MOD_VALUE_TYPES = str | int | float


@lru_cache(maxsize=512)
def _serialize_mod_cached(
    val: MOD_VALUE_TYPES,
    val_type: str,  # Add this!
    mult: int,
    brackets: str,
    include_plus: bool,
    precision: int | None,
) -> str:
    """Global cache for mod serialization."""
    # Serialize val (val_type is just for cache key differentiation)
    if isinstance(val, float):
        if precision is not None:
            val_str = f"{val:.{precision}f}"
        else:
            val_str = str(val)
            if "." not in val_str and "e" not in val_str.lower():
                val_str += ".0"
    else:
        val_str = str(val)

    # Add plus if needed
    if include_plus and isinstance(val, (int, float)) and val > 0:
        val_str = f"+{val_str}"

    # Format with multiplier
    if mult > 1:
        return f"{brackets[0]}{val_str}{brackets[1]}^{mult}"
    else:
        return f"{brackets[0]}{val_str}{brackets[1]}"


class Mod:
    """
    A modification with optional multiplier
    """

    def __init__(self, val: MOD_VALUE_TYPES, mult: int = 1):
        self.__val = val
        self.__mult = mult

    @property
    def val(self) -> MOD_VALUE_TYPES:
        return self.__val

    @property
    def mult(self) -> int:
        return self.__mult

    def serialize(
        self,
        brackets: str,
        include_plus: bool = False,
        precision: int | None = None,
    ) -> str:
        """Serialize the mod into a string (globally cached)."""
        return _serialize_mod_cached(
            self.val,
            type(self.val).__name__,  # Add type name to cache key!
            self.mult,
            brackets,
            include_plus,
            precision,
        )

    def __hash__(self) -> int:
        return hash((self.val, self.mult))

    def __repr__(self) -> str:
        # keep str quotes for strings but not for numbers
        if isinstance(self.val, str):
            return f"Mod('{self.val}', {self.mult})"
        return f"Mod({self.val}, {self.mult})"

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, (str, float, int)):
            other = Mod(other, 1)
        if self.val != other.val:
            return False
        if self.mult != other.mult:
            return False
        return True

    def __lt__(self, other: Any) -> bool:
        """
        Compare two mods (Doesn't Really matter, just for sorting)
        """
        if isinstance(other, (str, float, int)):
            other = Mod(other, 1)
        if str(self.val) < str(other.val):
            return True
        if self.val == other.val:
            return self.mult < other.mult
        return False

    def copy(self, deep: bool = True) -> Mod:
        """
        Create a copy of the mod
        """
        return self

    def to_dict(self) -> dict[str, str | float | int]:
        """
        Convert the mod to a dictionary
        """
        return {"val": self.val, "mult": self.mult}


def setup_mod(mod: MOD_VALUE_TYPES | Mod) -> Mod:
    if isinstance(mod, Mod):
        return mod
    if isinstance(mod, (str, int, float)):  # type: ignore
        return Mod(mod, 1)
    raise TypeError(f"Invalid mod input: {mod}")
