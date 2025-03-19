"""
This module provides functions to convert various input formats for modifications and intervals into standardized
data structures.
"""

from typing import List, Dict, Any, Union, Tuple

from peptacular.proforma.proforma_dataclasses import Mod, Interval

ModIndex = Union[int, str]
ModValue = Union[str, int, float, Mod]
ModDictValue = Union[List[Mod], List[Interval], int]
ModDict = Dict[ModIndex, ModDictValue]

ACCEPTED_MOD_INPUT = Union[List[ModValue], ModValue]
INTERVAL_VALUE = Union[Tuple[int, int, bool, Union[ACCEPTED_MOD_INPUT, None]], Interval]
ACCEPTED_INTERVAL_INPUT = Union[List[INTERVAL_VALUE], INTERVAL_VALUE]


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
    # Check if mod is already a Mod instance to avoid unnecessary conversion.
    if isinstance(mod, Mod):
        return mod

    if isinstance(mod, (str, int, float)):
        return Mod(mod, 1)

    raise ValueError(f"Invalid mod input: {mod}")


def fix_list_of_list_of_mods(mods: Union[List[List[ModValue]], List[ModValue], ModValue]) -> List[List[Mod]]:
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
    if isinstance(mods, list) and all(isinstance(mod, (str, int, float, Mod)) for mod in mods):
        return [fix_list_of_mods(mods)]

    # Case 3: mods is a MultiModList
    if isinstance(mods, list):
        return [fix_list_of_mods(sublist) for sublist in mods]

    raise ValueError(f"Invalid mod input: {mods}")


def remove_empty_list_of_list_of_mods(mods: List[List[Mod]]) -> Union[List[List[Mod]], None]:
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

    new_mods2 = []
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
    if isinstance(mods, list):

        # No change needed if all elements are already Mod instances
        if all(isinstance(mod, Mod) for mod in mods):
            return mods

        for i, _ in enumerate(mods):
            mods[i] = convert_to_mod(mods[i])

        return mods

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


def fix_dict_of_mods(mods: Dict[Any, Union[List[ModValue], ModValue]]) -> Dict[Any, List[Mod]]:
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


def fix_intervals_input(intervals: Union[List[INTERVAL_VALUE], INTERVAL_VALUE]) -> List[Interval]:
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

    if isinstance(intervals, list):
        if all(isinstance(interval, Interval) for interval in intervals):
            return intervals

        return [fix_interval_input(interval) for interval in intervals]

    raise ValueError(f"Invalid interval input: {intervals}")
