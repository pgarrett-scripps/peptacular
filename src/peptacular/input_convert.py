from typing import List, Dict, Any, Union

from peptacular.types import ModValue, IntervalValue
from peptacular.proforma.proforma_dataclasses import Mod, Interval


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
    return Mod(mod, 1)


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
    elif isinstance(mods, list) and (not mods or isinstance(mods[0], (str, int, float, Mod))):
        return [list(map(convert_to_mod, mods))]

    # Case 3: mods is a MultiModList
    elif isinstance(mods, list):
        return [list(map(convert_to_mod, sublist)) for sublist in mods]

    return []


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

        >>> fix_list_of_mods(['phospho'])
        [Mod('phospho', 1)]

        >>> fix_list_of_mods(['phospho', 1])
        [Mod('phospho', 1), Mod(1, 1)]

    """

    # Case 1: mods is a ModValue (str or Mod)
    if isinstance(mods, (str, int, float, Mod)):
        return [convert_to_mod(mods)]

    # Case 2: mods is a SingleModList
    elif isinstance(mods, list) and (not mods or isinstance(mods[0], (str, int, float, Mod))):
        return list(map(convert_to_mod, mods))

    return []


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

    """

    return {k: fix_list_of_mods(v) for k, v in mods.items()}


def fix_interval_input(interval: IntervalValue) -> Interval:
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


def fix_intervals_input(intervals: Union[List[IntervalValue], IntervalValue]) -> List[Interval]:
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

    if not isinstance(intervals, List):
        intervals = [intervals]

    new_intervals = []
    for i, interval in enumerate(intervals):
        new_intervals.append(fix_interval_input(interval))

    return new_intervals
