from __future__ import annotations

import copy
from collections import Counter
from dataclasses import dataclass
from typing import List, Any

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

    def __eq__(self, other: Mod) -> bool:
        if self.val != other.val:
            return False
        if self.mult != other.mult:
            return False
        return True

    def __lt__(self, other: Mod) -> bool:
        """
        Compare two mods (Doesn't Really matter, just for sorting)
        """
        if str(self.val) < str(other.val):
            return True
        if self.val == other.val:
            return self.mult < other.mult
        return False


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

    def __eq__(self, other: Interval) -> bool:
        if self.start != other.start:
            return False
        if self.end != other.end:
            return False
        if self.ambiguous != other.ambiguous:
            return False
        if not are_mods_equal(self.mods, other.mods):
            return False
        return True

    def __hash__(self):
        return hash((self.start, self.end, self.ambiguous, tuple(sorted(self.mods)) if self.mods else None))

    def has_mods(self) -> bool:
        return self.mods is not None


def are_mods_equal(mods1: List[Mod] | None, mods2: List[Mod] | None) -> bool:
    """
    Check if two lists of mods are equal

    :param mods1: List of mods
    :type mods1: List[Mod]

    :param mods2: List of mods
    :type mods2: List[Mod]

    :return: True if the lists are equal, False otherwise
    :rtype: bool

    .. code-block:: python

        >>> are_mods_equal([Mod('phospho', 1)], [Mod('phospho', 1)])
        True

        >>> are_mods_equal([Mod('phospho', 1)], [Mod('acetyl', 1)])
        False

        >>> are_mods_equal([Mod('phospho', 1)], [Mod('phospho', 2)])
        False

        >>> are_mods_equal([Mod('phospho', 1), Mod('acetyl', 2)], [Mod('acetyl', 2), Mod('phospho', 1)])
        True

        >>> are_mods_equal(None, None)
        True

        >>> are_mods_equal(None, [Mod('phospho', 1)])
        False

    """

    if mods1 is None and mods2 is None:
        return True

    if mods1 is None and mods2 is not None:
        return False

    if mods1 is not None and mods2 is None:
        return False

    return Counter(mods1) == Counter(mods2)


def are_intervals_equal(intervals1: List[Interval] | None, intervals2: List[Interval] | None) -> bool:
    """
    Check if two lists of intervals are equal

    :param intervals1: List of intervals
    :type intervals1: List[Interval]

    :param intervals2: List of intervals
    :type intervals2: List[Interval]

    :return: True if the lists are equal, False otherwise
    :rtype: bool

    .. code-block:: python

        >>> are_intervals_equal([Interval(1, 2, False, [Mod(1, 1)])], [Interval(1, 2, False, [Mod(1, 1)])])
        True

        >>> intervals1 = [Interval(1, 2, False, [Mod(1, 1)]), Interval(3, 4, False, [Mod(1, 1)])]
        >>> intervals2 = [Interval(1, 2, False, [Mod(1, 1)]), Interval(3, 4, False, [Mod(1, 1)])]
        >>> are_intervals_equal(intervals1, intervals2)
        True

        >>> intervals1 = [Interval(1, 20, False, [Mod(1, 1)]), Interval(3, 4, False, [Mod(1, 1)])]
        >>> intervals2 = [Interval(1, 2, False, [Mod(1, 1)]), Interval(3, 4, False, [Mod(1, 1)])]
        >>> are_intervals_equal(intervals1, intervals2)
        False

        >>> are_intervals_equal(None, None)
        True

        >>> are_intervals_equal(None, [Interval(1, 2, False, [Mod(1, 1)])])
        False

        >>> are_intervals_equal([Interval(1, 2, False, None)], [Interval(1, 2, False, [Mod(1, 1)])])
        False

    """

    if intervals1 is None and intervals2 is None:
        return True

    if intervals1 is None and intervals2 is not None:
        return False

    if intervals1 is not None and intervals2 is None:
        return False

    if len(intervals1) != len(intervals2):
        return False

    return Counter(intervals1) == Counter(intervals2)
