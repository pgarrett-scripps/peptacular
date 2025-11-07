from __future__ import annotations

import sys
from functools import lru_cache
from typing import Any, NamedTuple

MOD_VALUE_TYPES = str | int | float


@lru_cache(maxsize=512)
def _serialize_mod_cached(
    val: MOD_VALUE_TYPES,
    val_type: str,
    mult: int,
    brackets: str | None,
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

    if brackets is None:
        brackets_tuple: tuple[str, str] = ("", "")
        if mult > 1:
            raise ValueError("Brackets must be provided if multiplier is greater than 1.")
    else:
        if len(brackets) != 2:
            raise ValueError("Brackets string must be of length 2.")
        brackets_tuple = (brackets[0], brackets[1])

    # Format with multiplier
    if mult > 1:
        return f"{brackets_tuple[0]}{val_str}{brackets_tuple[1]}^{mult}"
    else:
        return f"{brackets_tuple[0]}{val_str}{brackets_tuple[1]}"


class Mod(NamedTuple):
    """
    A modification with optional multiplier
    """
    val: MOD_VALUE_TYPES
    mult: int = 1

    def serialize(
        self,
        brackets: str | None,
        include_plus: bool = False,
        precision: int | None = None,
    ) -> str:
        """Serialize the mod into a string (globally cached)."""
        return _serialize_mod_cached(
            self.val,
            type(self.val).__name__,
            self.mult,
            brackets,
            include_plus,
            precision,
        )
    
    @staticmethod
    def parse(mod_str: str, brackets: str | None) -> "Mod":
        """
        Parse a mod string into a Mod object (returns cached instance)
        """
        # Handle multiplier
        if brackets is not None and len(brackets) != 2:
            raise ValueError("Brackets string must be of length 2.")
        
        mult = 1
        val_str = mod_str
        
        # Remove brackets if present
        if brackets is not None:
            open_br, close_br = brackets[0], brackets[1]
            if mod_str.startswith(open_br) and close_br in mod_str:
                # Find the closing bracket and check for multiplier
                close_idx = mod_str.rfind(close_br)
                val_str = mod_str[len(open_br):close_idx]
                
                # Check for ^multiplier after closing bracket
                remainder = mod_str[close_idx + len(close_br):]
                if remainder.startswith('^'):
                    mult = int(remainder[1:])
            else:
                val_str = mod_str
        
        # Remove leading plus sign if present
        if val_str.startswith('+'):
            val_str = val_str[1:]
        
        # Try to parse as float or int
        try:
            if '.' in val_str or 'e' in val_str.lower():
                val = float(val_str)
            else:
                val = int(val_str)
        except ValueError:
            # If parsing fails, treat as string
            val = val_str
        
        return get_mod(val, mult)

    def __hash__(self) -> int:
        return hash((self.val, self.mult))

    def __repr__(self) -> str:
        # keep str quotes for strings but not for numbers
        if isinstance(self.val, str):
            return f"Mod('{self.val}', {self.mult})"
        return f"Mod({self.val}, {self.mult})"

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, (str, float, int)):
            other = get_mod(other, 1)
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
            other = get_mod(other, 1)
        if str(self.val) < str(other.val):
            return True
        if self.val == other.val:
            return self.mult < other.mult
        return False

    def copy(self, deep: bool = True) -> "Mod":
        """
        Create a copy of the mod (returns same cached instance)
        """
        return self

    def to_dict(self) -> dict[str, str | float | int]:
        """
        Convert the mod to a dictionary
        """
        return {"val": self.val, "mult": self.mult}


@lru_cache(maxsize=2048)
def _get_mod_impl(val: MOD_VALUE_TYPES, mult: int, val_type: type) -> Mod:
    """
    Internal implementation with type-aware caching.
    
    The val_type parameter ensures int(1) and float(1.0) don't collide in cache.
    """
    # Intern strings for better memory sharing
    if isinstance(val, str):
        val = sys.intern(val)
    return Mod(val, mult)


def get_mod(val: MOD_VALUE_TYPES, mult: int = 1) -> Mod:
    """
    Get a cached Mod instance. Returns the same object for identical values.
    
    This ensures memory efficiency by reusing Mod instances across the application.
    String values are automatically interned for better memory sharing.
    
    Note: Cache is type-aware - get_mod(1, 1) and get_mod(1.0, 1) return different Mods.
    """
    return _get_mod_impl(val, mult, type(val))


def setup_mod(mod: MOD_VALUE_TYPES | Mod) -> Mod:
    """
    Convert input to a Mod, using cached instances for memory efficiency.
    """
    if isinstance(mod, Mod):
        # Canonicalize through cache to ensure same object is used
        return get_mod(mod.val, mod.mult)
    if isinstance(mod, (str, int, float)):
        return get_mod(mod, 1)
    raise TypeError(f"Invalid mod input: {mod}")


def clear_mod_cache() -> None:
    """
    Clear the mod cache. Useful for testing or memory management.
    """
    _get_mod_impl.cache_clear()
    _serialize_mod_cached.cache_clear()


def get_mod_cache_info() -> dict[str, Any]:
    """
    Get cache statistics for debugging and monitoring.
    
    Returns:
        Dictionary with cache hits, misses, size, and maxsize for both caches
    """
    mod_info = _get_mod_impl.cache_info()
    serialize_info = _serialize_mod_cached.cache_info()
    
    return {
        "mod_cache": {
            "hits": mod_info.hits,
            "misses": mod_info.misses,
            "size": mod_info.currsize,
            "maxsize": mod_info.maxsize,
            "hit_rate": mod_info.hits / (mod_info.hits + mod_info.misses) if (mod_info.hits + mod_info.misses) > 0 else 0.0
        },
        "serialize_cache": {
            "hits": serialize_info.hits,
            "misses": serialize_info.misses,
            "size": serialize_info.currsize,
            "maxsize": serialize_info.maxsize,
            "hit_rate": serialize_info.hits / (serialize_info.hits + serialize_info.misses) if (serialize_info.hits + serialize_info.misses) > 0 else 0.0
        }
    }