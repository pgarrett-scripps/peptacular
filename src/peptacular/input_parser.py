from typing import List, Tuple, Dict

from peptacular.types import ModValue, Mod
from peptacular.dataclasses import Interval
from peptacular.types import ACCEPTED_MOD_INPUT


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


def fix_list_of_list_of_mods(mods: List[List[ModValue]] | List[ModValue] | ModValue) -> List[List[Mod]]:
    """
    Convert the input mods to a list of lists of Mod instances.

    :param mods:
    :return:

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


def fix_list_of_mods(mods: List[ModValue] | ModValue) -> List[Mod]:
    """
    Convert the input mods to a list of lists of Mod instances.

    :param mods:
    :return:

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


def fix_dict_of_mods(mods: Dict[int | str, List[ModValue] | ModValue]) -> Dict[int | str, List[Mod]]:
    """
    Convert the input mods to a dictionary of lists of Mod instances.

    :param mods:
    :return:

    .. code-block:: python

        >>> fix_dict_of_mods({'phospho': 'phospho'})
        {'phospho': [Mod('phospho', 1)]}

        >>> fix_dict_of_mods({'phospho': [3.0]})
        {'phospho': [Mod(3.0, 1)]}

        >>> fix_dict_of_mods({'phospho': ['phospho']})
        {'phospho': [Mod('phospho', 1)]}

        >>> fix_dict_of_mods({'phospho': ['phospho', 1]})
        {'phospho': [Mod('phospho', 1), Mod(1, 1)]}

    """

    return {k: fix_list_of_mods(v) for k, v in mods.items()}


def fix_interval_input(interval: Tuple[int, int, bool, ACCEPTED_MOD_INPUT | None] | Interval) -> Interval:
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

    """

    # parse mod

    mods = fix_list_of_mods(interval[3]) if interval[3] else None
    return Interval(interval[0], interval[1], interval[2], mods)
