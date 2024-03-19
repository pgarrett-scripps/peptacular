from __future__ import annotations

from typing import Dict

from peptacular.chem_util import write_chem_formula
from peptacular.mods.mod_setup import _glycan_comp, _parse_glycan_formula  # To avoid circular import
from peptacular.types import Chem_Composition


def write_glycan_formula(glycan_dict: Chem_Composition, sep: str = None) -> str:
    """
    Writes a glycan dictionary to a string.

    :param glycan_dict: A dictionary containing the glycan components and their counts.
    :type glycan_dict: dict
    :param sep: The separator to use between the glycan component and its count.
    :type sep: str

    :return: A glycan formula string.
    :rtype: str

    .. code-block:: python

            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu': 1})
            'HexNAc2Hex3Neu1'

            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            'HexNAc2Hex3Neu5Gc1'

            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': -1})
            'HexNAc2Hex3Neu5Gc-1'

            >>> write_glycan_formula({})
            ''

            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1}, sep=' ')
            'HexNAc 2 Hex 3 Neu5Gc 1'

    """

    if sep is not None:
        return sep.join([f'{component}{sep}{count}' for component, count in glycan_dict.items()])

    # Convert the dictionary to a list of strings
    glycan_str = [f'{component}{count}' for component, count in glycan_dict.items()]

    # Join the list of strings into a single string
    return ''.join(glycan_str)


def glycan_comp(glycan: Dict[str, int] | str) -> Chem_Composition:
    """
    Converts a glycan dictionary to a chemical formula.

    :param glycan: A dictionary containing the glycan components and their counts, or a glycan formula string.
    :type glycan: dict | str

    :return: A chemical formula string.
    :rtype: str

    .. code-block:: python

            >>> _glycan_comp({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            {'C': 45, 'H': 73, 'N': 3, 'O': 34}

            >>> _glycan_comp({'HexNAc': 2.2, 'Hex': 3, 'Neu5Gc': 1})
            {'C': 46.6, 'H': 75.6, 'N': 3.2, 'O': 35.0}

    """
    return _glycan_comp(glycan)


def parse_glycan_formula(formula: str) -> Chem_Composition:
    """
    Parses a glycan sequence into its constituent parts.

    :param formula: A glycan formula string.
    :type formula: str

    :raises UnknownGlycanError: If the glycan formula contains an unknown glycan.

    :return: A dictionary containing the glycan components and their counts.
    :rtype: dict

    .. code-block:: python

            >>> _parse_glycan_formula('HexNAc2Hex3Neu1')
            {'HexNAc': 2, 'Hex': 3, 'Neu': 1}

            >>> _parse_glycan_formula('HexNAc2Hex3Neu5Gc1')
            {'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1}

            >>> _parse_glycan_formula('HexNAc2Hex3Neu5Gc-1')
            {'HexNAc': 2, 'Hex': 3, 'Neu5Gc': -1}

            >>> _parse_glycan_formula('HexNAc2.2Hex3.9Neu1')
            {'HexNAc': 2.2, 'Hex': 3.9, 'Neu': 1}

            >>> _parse_glycan_formula('')
            {}

            >>> _parse_glycan_formula('Hex')
            {'Hex': 1}

            # This will raise an UnknownGlycanError
            >>> _parse_glycan_formula('HeSNAc2Hex3Neu5Gc1X')
            Traceback (most recent call last):
            peptacular.errors.UnknownGlycanError: Unknown glycan: HeSNAc2Hex3Neu5Gc1X

    """
    return _parse_glycan_formula(formula)


def convert_glycan_formula_to_chem_formula(glycan: Dict[str, int] | str) -> str:
    """
    Converts a glycan dictionary to a chemical formula.

    :param glycan: A dictionary containing the glycan components and their counts, or a glycan formula string.
    :type glycan: dict | str

    :return: A chemical formula string.
    :rtype: str

    .. code-block:: python

            >>> convert_glycan_formula_to_chem_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            'C45H73N3O34'

    """

    return write_chem_formula(glycan_comp(glycan))
