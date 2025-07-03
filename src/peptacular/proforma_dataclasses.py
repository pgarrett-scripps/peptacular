"""
proforma_dataclasses.py
"""

import copy
from collections import Counter
from dataclasses import dataclass
from typing import List, Union, Any, Optional, Dict
from typing import List, Dict, Any, Union, Tuple

from .utils2 import convert_type



@dataclass
class Mod:
    """
    A modification with optional multiplier
    """

    val: Union[str, float, int]
    mult: int = 1
    tag: Optional[str] = None
    loc_score: Optional[float] = None
    dislpay_val: bool = True

    def __post_init__(self) -> None:

        if isinstance(self.val, str):
            self.val = convert_type(self.val)

    def flatten(self) -> List[Union[str, float, int]]:
        """
        Flatten the mod into a list of repeated values based on the multiplier
        """
        return [self.val] * self.mult

    def _serialize_val(self, precision: Optional[float] = None) -> str:
        s = ""
        if self.dislpay_val:
            if precision is not None:
                if isinstance(self.val, float):
                    s = f"{self.val:.{precision}f}"
                elif isinstance(self.val, int):
                    s = str(self.val)
                else:
                    s = str(self.val)
            else:
                s = str(self.val)

        if self.tag is not None:
            s = f"{s}#{self.tag}"
            if self.loc_score is not None:
                s = f"{s}({self.loc_score})"

        return s

    def serialize(
        self,
        brackets: str,
        include_plus: bool = False,
        precision: Optional[float] = None,
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

    def dict(self) -> Dict[str, Union[str, float, int]]:
        """
        Convert the mod to a dictionary
        """
        return {"val": self.val, "mult": self.mult}

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

    def copy(self, deep: bool = True) -> "Mod":
        """
        Create a copy of the mod
        """

        if deep:
            return copy.deepcopy(self)

        return copy.copy(self)


class Interval:
    """
    A sequence interval with optional modifications
    """

    def __init__(
        self, start: int, end: int, ambiguous: bool, mods: Optional[List[Mod]] = None
    ) -> None:
        """
        Initialize an Interval object

        :param start: Start position of the interval (inclusive)
        :type start: int

        :param end: End position of the interval (exclusive)
        :type end: int

        :param ambiguous: Whether the interval is ambiguous
        :type ambiguous: bool

        :param mods: Optional list of modifications for the interval
        :type mods: Optional[List[Mod]]
        """

        self.start = start
        self.end = end
        self.ambiguous = ambiguous
        if mods is not None:
            self._mods = fix_list_of_mods(mods) # type: ignore
        else:
            self._mods = []

    @property
    def mods(self) -> List[Mod]:
        """
        Get the modifications of the interval
        """
        return self._mods

    @mods.setter
    def mods(self, value: Optional[List[Mod]]):
        """
        Set the modifications of the interval
        """
        if value is None:
            self._mods = []
        else:
            self._mods = fix_list_of_mods(value) # type: ignore

    def has_mods(self) -> bool:
        """
        Check if the interval has modifications
        """
        return len(self.mods) > 0

    def dict(self) -> Dict[str, Any]:
        """
        Convert the interval to a dictionary
        """
        result: Dict[str, Any] = {
            "start": self.start,
            "end": self.end,
            "ambiguous": self.ambiguous,
            "mods": copy.deepcopy(self.mods),
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
        return hash(
            (
                self.start,
                self.end,
                self.ambiguous,
                tuple(sorted(self.mods)) if self.mods else None,
            )
        )

    def copy(self, deep: bool = True) -> "Interval":
        """
        Create a copy of the interval
        """
        if deep:
            return copy.deepcopy(self)
        return copy.copy(self)


Span = Tuple[int, int, int]

# Chem Composition Type
ChemComposition = Dict[str, Union[int, float]]

ModIndex = Union[int, str]
ModValue = Union[str, int, float, Mod]
ModDictValue = Union[List[Mod], List[Interval], int]
ModDict = Dict[ModIndex, ModDictValue]

ACCEPTED_MOD_INPUT = Union[List[ModValue], ModValue]
INTERVAL_VALUE = Union[Tuple[int, int, bool, Union[ACCEPTED_MOD_INPUT, None]], Interval]
ACCEPTED_INTERVAL_INPUT = Union[List[INTERVAL_VALUE], INTERVAL_VALUE]


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

    if mods1 is None or mods2 is None:
        return mods1 is mods2

    return Counter(mods1) == Counter(mods2)


def are_intervals_equal(
    intervals1: Optional[List[Interval]], intervals2: Optional[List[Interval]]
) -> bool:
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

    if len(intervals1) != len(intervals2): # type: ignore
        return False

    return Counter(intervals1) == Counter(intervals2)


def convert_to_mod(mod: ModValue) -> Mod:
    """
    Convert the input mod to a Mod instance.

    :param mod:
    :return:

    .. code-block:: python

        >>> convert_to_mod('phospho')
        Mod('phospho', 1)

        >>> convert_to_mod(3.0)
        Mod(3.0, 1)

        >>> convert_to_mod(Mod('phospho', 1))
        Mod('phospho', 1)

    """
    if isinstance(mod, Mod):
        return mod

    if isinstance(mod, (str, int, float)): # type: ignore
        return Mod(mod, 1)

    raise ValueError(f"Invalid mod input: {mod}")


def fix_list_of_list_of_mods(
    mods: Union[List[List[ModValue]], List[ModValue], ModValue],
) -> List[List[Mod]]:
    """
    Convert the input mods to a list of lists of Mod instances.

    :param mods: Either a list of lists of mods, a list of mods, or a single mod.
    :type mods: Union[List[List[ModValue]], List[ModValue], ModValue]

    :return: List of lists of Mod instances
    :return: List[List[Mod]]

    .. code-block:: python

        >>> fix_list_of_list_of_mods('phospho')
        [[Mod('phospho', 1)]]

        >>> fix_list_of_list_of_mods([3.0])
        [[Mod(3.0, 1)]]

        >>> fix_list_of_list_of_mods(['phospho'])
        [[Mod('phospho', 1)]]

        >>> fix_list_of_list_of_mods(['phospho', 1])
        [[Mod('phospho', 1), Mod(1, 1)]]

        >>> fix_list_of_list_of_mods([['phospho']])
        [[Mod('phospho', 1)]]

        >>> fix_list_of_list_of_mods([['phospho', 1]])
        [[Mod('phospho', 1), Mod(1, 1)]]

        >>> fix_list_of_list_of_mods([['phospho', 1], ['acetyl', 2]])
        [[Mod('phospho', 1), Mod(1, 1)], [Mod('acetyl', 1), Mod(2, 1)]]
    """

    # Case 1: mods is a ModValue (str or Mod)
    if isinstance(mods, (str, int, float, Mod)):
        return [[convert_to_mod(mods)]]

    # Case 2: mods is a SingleModList
    if isinstance(mods, list) and all( # type: ignore
        isinstance(mod, (str, int, float, Mod)) for mod in mods
    ):
        return [fix_list_of_mods(mods)] # type: ignore

    # Case 3: mods is a MultiModList
    if isinstance(mods, list):  # type: ignore
        return [fix_list_of_mods(sublist) for sublist in mods]

    raise ValueError(f"Invalid mod input: {mods}")


def remove_empty_list_of_list_of_mods(
    mods: List[List[Mod]],
) -> Union[List[List[Mod]], None]:
    """
    Remove empty lists from a list of lists of Mod instances.

    :param mods: List of lists of Mod instances
    :type mods: List[List[Mod]]

    :return: List of lists of Mod instances
    :rtype: List[List[Mod]]

    .. code-block:: python

        >>> remove_empty_list_of_list_of_mods([[Mod('phospho', 1)], [Mod('acetyl', 1)], []])
        [[Mod('phospho', 1)], [Mod('acetyl', 1)]]

        >>> remove_empty_list_of_list_of_mods([[Mod('phospho', 1)], [], [Mod('acetyl', 1)]])
        [[Mod('phospho', 1)], [Mod('acetyl', 1)]]

        >>> remove_empty_list_of_list_of_mods([[], [], []])

    """

    new_mods2: List[List[Mod]] = []
    for mod_list in mods:
        if mod_list:
            new_mods2.append(mod_list)

    if new_mods2:
        return new_mods2

    return None


def fix_list_of_mods(mods: Union[List[ModValue], ModValue]) -> List[Mod]:
    """
    Convert the input mods to a list of lists of Mod instances.

    :param mods: Either a list of mods or a single mod.
    :type mods: Union[List[ModValue], ModValue]

    :return: List of Mod instances
    :rtype: List[Mod]

    .. code-block:: python

        >>> fix_list_of_mods('phospho')
        [Mod('phospho', 1)]

        >>> fix_list_of_mods([3.0])
        [Mod(3.0, 1)]

        >>> fix_list_of_mods([])
        []

        >>> fix_list_of_mods(['phospho'])
        [Mod('phospho', 1)]

        >>> fix_list_of_mods(['phospho', Mod(1, 1)])
        [Mod('phospho', 1), Mod(1, 1)]

    """

    # Case 1: mods is a ModValue (str or Mod)
    if isinstance(mods, (str, int, float, Mod)):
        return [convert_to_mod(mods)]

    # Case 2: mods is a list of ModValues
    if isinstance(mods, list): # type: ignore

        # No change needed if all elements are already Mod instances
        if all(isinstance(mod, Mod) for mod in mods):
            return mods # type: ignore

        for i, _ in enumerate(mods):
            mods[i] = convert_to_mod(mods[i])

        return mods # type: ignore

    # Case 3: mods is an invalid input
    raise ValueError(f"Invalid mod input: {mods}")


def remove_empty_list_of_mods(mods: List[Mod]) -> Union[List[Mod], None]:
    """
    Remove empty lists from a list of Mod instances.

    :param mods: List of Mod instances
    :type mods: List[Mod]

    :return: List of Mod instances
    :rtype: List[Mod]

    .. code-block:: python

        >>> remove_empty_list_of_mods([Mod('phospho', 1), Mod('acetyl', 1), Mod('acetyl', 1)])
        [Mod('phospho', 1), Mod('acetyl', 1), Mod('acetyl', 1)]

        >>> remove_empty_list_of_mods([])

    """

    return mods if mods else None


def fix_dict_of_mods(
    mods: Dict[Any, Union[List[ModValue], ModValue]],
) -> Dict[Any, List[Mod]]:
    """
    Convert the input mods to a dictionary of lists of Mod instances. Mainly used to convert internal mod input to
    the correct format. This will not work when used on the whole mod dict.

    :param mods: Dictionary of mods
    :type mods: Dict[Any, Union[List[ModValue], ModValue]]

    :return: Dictionary of lists of Mod instances
    :rtype: Dict[Any, List[Mod]]

    .. code-block:: python

        >>> fix_dict_of_mods({1: 'phospho'})
        {1: [Mod('phospho', 1)]}

        >>> fix_dict_of_mods({2: [3.0]})
        {2: [Mod(3.0, 1)]}

        >>> fix_dict_of_mods({2: []})
        {2: []}

    """

    return {k: fix_list_of_mods(v) for k, v in mods.items()}


def fix_interval_input(interval: INTERVAL_VALUE) -> Interval:
    """
    Convert the input intervals to a list of Interval instances.

    :param interval:
    :return:

    .. code-block:: python

        >>> fix_interval_input((1, 2, False, 'phospho'))
        Interval(1, 2, False, [Mod('phospho', 1)])

        >>> fix_interval_input((1, 2, False, [3.0]))
        Interval(1, 2, False, [Mod(3.0, 1)])

        >>> fix_interval_input((1, 2, False, ['phospho']))
        Interval(1, 2, False, [Mod('phospho', 1)])

        >>> fix_interval_input(Interval(1, 2, False, [Mod('phospho', 1)]))
        Interval(1, 2, False, [Mod('phospho', 1)])

    """

    # parse mod
    if isinstance(interval, Interval):
        return interval

    mods = fix_list_of_mods(interval[3]) if interval[3] else None
    return Interval(interval[0], interval[1], interval[2], mods)


def fix_intervals_input(
    intervals: Union[List[INTERVAL_VALUE], INTERVAL_VALUE],
) -> List[Interval]:
    """
    Convert the input intervals to a list of Interval instances.

    :param intervals: List of intervals
    :param intervals: Union[List[IntervalValue], IntervalValue]

    :return: List of Interval instances
    :rtype: List[Interval]

    .. code-block:: python

        >>> fix_intervals_input([(1, 2, False, 'phospho')])
        [Interval(1, 2, False, [Mod('phospho', 1)])]

        >>> fix_intervals_input((1, 2, False, [3.0]))
        [Interval(1, 2, False, [Mod(3.0, 1)])]

        >>> fix_intervals_input([(1, 2, False, ['phospho'])])
        [Interval(1, 2, False, [Mod('phospho', 1)])]

        >>> fix_intervals_input([Interval(1, 2, False, [Mod('phospho', 1)])])
        [Interval(1, 2, False, [Mod('phospho', 1)])]

    """

    if isinstance(intervals, (Tuple, Interval)):
        return [fix_interval_input(intervals)]

    if isinstance(intervals, list): # type: ignore
        if all(isinstance(interval, Interval) for interval in intervals):
            return intervals  # type: ignore

        return [fix_interval_input(interval) for interval in intervals]

    raise ValueError(f"Invalid interval input: {intervals}")
