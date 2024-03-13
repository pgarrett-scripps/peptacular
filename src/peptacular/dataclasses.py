from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import List

from peptacular.util import convert_type


@dataclass
class Mod:
    """
    A modification with optional multiplier
    """
    val: int | float | str
    mult: int

    def __post_init__(self):
        self.val = convert_type(self.val)

    def flatten(self) -> List[int | float | str]:
        return [self.val] * self.mult

    def serialize(self, brackets: str) -> str:
        return f"{brackets[0]}{self.val}{brackets[1]}^{self.mult}" if self.mult > 1 else \
            f"{brackets[0]}{self.val}{brackets[1]}"

    def __hash__(self):
        return hash((self.val, self.mult))

    def __repr__(self):
        # keep str quotes for strings but not for numbers
        if isinstance(self.val, str):
            return f"Mod('{self.val}', {self.mult})"
        return f"Mod({self.val}, {self.mult})"

    def dict(self):
        return {
            "val": self.val,
            "mult": self.mult
        }


@dataclass
class Interval:
    """
    A sequence interval with optional modifications
    """
    start: int
    end: int
    ambiguous: bool
    mods: List[Mod] | None = None

    def dict(self):
        result = {
            "start": self.start,
            "end": self.end,
            "ambiguous": self.ambiguous,
            "mods": copy.deepcopy(self.mods)
        }
        return result

    def __repr__(self):
        return f"Interval({self.start}, {self.end}, {self.ambiguous}, {self.mods})"