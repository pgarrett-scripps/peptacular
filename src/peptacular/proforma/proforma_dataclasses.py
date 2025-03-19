"""
proforma_dataclasses.py
"""

import copy
from collections import Counter
from dataclasses import dataclass
from typing import List, Union, Any, Optional, Dict

from peptacular.util import convert_type


@dataclass
class Mod:
    """
    A modification with optional multiplier
    """
    val: Union[str, float, int]
    mult: int

    def __post_init__(self) -> None:
        self.val = convert_type(self.val)

    def flatten(self) -> List[Union[str, float, int]]:
        """
        Flatten the mod into a list of repeated values based on the multiplier
        """
        return [self.val] * self.mult

    def serialize(self, brackets: str, include_plus: bool = False) -> str:
        """
        Serialize the mod into a string
        """
        # Determine if the value is positive and prefix '+' for positive numbers
        if include_plus is True:
            val_str = f"+{self.val}" if isinstance(self.val, (int, float)) and self.val > 0 else str(self.val)
        else:
            val_str = str(self.val)

        # Return the formatted string based on the multiplier value
        return f"{brackets[0]}{val_str}{brackets[1]}^{self.mult}" if self.mult > 1 else \
            f"{brackets[0]}{val_str}{brackets[1]}"

    def __hash__(self) -> int:
        return hash((self.val, self.mult))

    def __repr__(self) -> str:
        # keep str quotes for strings but not for numbers
        if isinstance(self.val, str):
            return f"Mod('{self.val}', {self.mult})"
        return f"Mod({self.val}, {self.mult})"

    def dict(self) -> Dict[str, Union[str, float, int]]:
        """
        Convert the mod to a dictionary
        """
        return {
            "val": self.val,
            "mult": self.mult
        }

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


@dataclass
class Interval:
    """
    A sequence interval with optional modifications
    """
    start: int
    end: int
    ambiguous: bool
    mods: Optional[List[Mod]] = None

    def dict(self) -> Dict[str, Any]:
        """
        Convert the interval to a dictionary
        """
        result = {
            "start": self.start,
            "end": self.end,
            "ambiguous": self.ambiguous,
            "mods": copy.deepcopy(self.mods)
        }
        return result

    def __repr__(self) -> str:
        return f"Interval({self.start}, {self.end}, {self.ambiguous}, {self.mods})"

    def __eq__(self, other: Any) -> bool:
        if self.start != other.start:
            return False
        if self.end != other.end:
            return False
        if self.ambiguous != other.ambiguous:
            return False
        if not are_mods_equal(self.mods, other.mods):
            return False
        return True

    def __hash__(self) -> int:
        return hash((self.start, self.end, self.ambiguous, tuple(sorted(self.mods)) if self.mods else None))

    def has_mods(self) -> bool:
        """
        Check if the interval has modifications
        """
        return self.mods is not None


def are_mods_equal(mods1: Optional[List[Mod]], mods2: Optional[List[Mod]]) -> bool:
    """
    Check if two lists of mods are equal

    :param mods1: List of mods. Can be None.
    :type mods1: Optional[List[Mod]]

    :param mods2: List of mods. Can be None.
    :type mods2: Optional[List[Mod]]

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


def are_intervals_equal(intervals1: Optional[List[Interval]], intervals2: Optional[List[Interval]]) -> bool:
    """
    Check if two lists of intervals are equal

    :param intervals1: List of intervals. Can be None.
    :type intervals1: Optional[List[Interval]]

    :param intervals2: List of intervals. Can be None.
    :type intervals2: Optional[List[Interval]]

    :return: True if the lists are equal, False otherwise
    :rtype: bool

    .. code-block:: python

        >>> are_intervals_equal([Interval(1, 2, False, [Mod(1, 1)])], [Interval(1, 2, False, [Mod(1, 1)])])
        True

        >>> ints1 = [Interval(1, 2, False, [Mod(1, 1)]), Interval(3, 4, False, [Mod(1, 1)])]
        >>> ints2 = [Interval(1, 2, False, [Mod(1, 1)]), Interval(3, 4, False, [Mod(1, 1)])]
        >>> are_intervals_equal(ints1, ints2)
        True

        >>> ints1 = [Interval(1, 20, False, [Mod(1, 1)]), Interval(3, 4, False, [Mod(1, 1)])]
        >>> ints2 = [Interval(1, 2, False, [Mod(1, 1)]), Interval(3, 4, False, [Mod(1, 1)])]
        >>> are_intervals_equal(ints1, ints2)
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
