from __future__ import annotations
import copy
from typing import Any

MOD_VALUE_TYPES = str | int | float


class Mod:
    """
    A modification with optional multiplier
    """

    def __init__(self, val: MOD_VALUE_TYPES, mult: int = 1):
        self.val = val
        self.mult = mult

    def _serialize_val(self, precision: int | None = None) -> str:
        s = ""
        if precision is not None:
            if isinstance(self.val, float):
                s = f"{self.val:.{precision}f}"
            elif isinstance(self.val, int):
                s = str(self.val)
            else:
                s = str(self.val)
        else:
            s = str(self.val)

        return s

    def serialize(
        self,
        brackets: str,
        include_plus: bool = False,
        precision: int | None = None,
    ) -> str:
        """
        Serialize the mod into a string
        """
        # Determine if the value is positive and prefix '+' for positive numbers
        if include_plus is True:
            val_str = (
                f"+{self._serialize_val(precision)}"
                if isinstance(self.val, (int, float)) and self.val > 0
                else self._serialize_val(precision)
            )
        else:
            val_str = self._serialize_val(precision)

        # Return the formatted string based on the multiplier value
        return (
            f"{brackets[0]}{val_str}{brackets[1]}^{self.mult}"
            if self.mult > 1
            else f"{brackets[0]}{val_str}{brackets[1]}"
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

        if deep:
            return copy.deepcopy(self)
        return copy.copy(self)

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
