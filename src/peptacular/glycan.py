"""
Glycan.py
"""

from typing import Union

from peptacular.types import ChemComposition
from peptacular.chem.chem_util import write_chem_formula
from peptacular.errors import InvalidGlycanFormulaError
from peptacular.mods.mod_db_setup import _glycan_comp, _parse_glycan_formula  # To avoid circular import


def write_glycan_formula(glycan_dict: ChemComposition, sep: str = '') -> str:
    """
    Writes a glycan dictionary to a string.

    :param glycan_dict: A dictionary containing the glycan components and their counts.
    :type glycan_dict: ChemComposition
    :param sep: The separator to use between the glycan component and its count. Default is an empty string.
    :type sep: str

    :return: A glycan formula string.
    :rtype: str

    .. code-block:: python

            # Int Counts
            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu': 1})
            'HexNAc2Hex3Neu1'

            # Float counts
            >>> write_glycan_formula({'HexNAc': 2.1, 'Hex': 3, 'Neu5Gc': 1})
            'HexNAc2.1Hex3Neu5Gc1'

            # Negative counts
            >>> write_glycan_formula({'HexNAc': -2.1, 'Hex': 3, 'Neu5Gc': -1})
            'HexNAc-2.1Hex3Neu5Gc-1'

            # Empty input
            >>> write_glycan_formula({})
            ''

            # Using a custom separator
            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1}, sep=' ')
            'HexNAc 2 Hex 3 Neu5Gc 1'

    """
    return sep.join([f'{component}{sep}{count}' for component, count in glycan_dict.items()])


def glycan_comp(glycan: Union[ChemComposition, str]) -> ChemComposition:
    """
    Converts a glycan dictionary to a chemical formula.

    :param glycan: A dictionary containing the glycan components and their counts, or a glycan formula string.
    :type glycan: Union[ChemComposition, str]

    :raises InvalidGlycanFormulaError: If the glycan formula is invalid.

    :return: A dictionary containing the glycan components and their counts.
    :rtype: ChemComposition

    .. code-block:: python

            # Int Counts
            >>> glycan_comp({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            {'C': 45, 'H': 73, 'N': 3, 'O': 34}

            # Float counts
            >>> glycan_comp({'HexNAc': 2.2, 'Hex': 3, 'Neu5Gc': 1})
            {'C': 46.6, 'H': 75.6, 'N': 3.2, 'O': 35.0}

            # Negative counts
            >>> glycan_comp({'HexNAc': -2, 'Hex': -3, 'Neu5Gc': -1})
            {'C': -45, 'H': -73, 'N': -3, 'O': -34}

            # Unknown Glycan
            >>> glycan_comp('HexXX')
            Traceback (most recent call last):
            peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "HexXX". Unknown glycan: "XX"!

            # Invalid Count
            >>> glycan_comp('Hex2.1.')
            Traceback (most recent call last):
            peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "Hex2.1.". Invalid count: "2.1."!

    """
    return _glycan_comp(glycan)


def parse_glycan_formula(formula: str, sep: str = '') -> ChemComposition:
    """
    Parses a glycan sequence into its constituent parts.

    :param formula: A glycan formula string.
    :type formula: str
    :param sep: The separator to use between the glycan component and its count. Default is an empty string.
    :type sep: str

    :raises InvalidGlycanFormulaError: If the glycan formula is invalid.

    :return: A dictionary containing the glycan components and their counts.
    :rtype: ChemComposition

    .. code-block:: python

            # Int Counts
            >>> parse_glycan_formula('HexNAc2Hex3Neu1')
            {'HexNAc': 2, 'Hex': 3, 'Neu': 1}

            # Float counts
            >>> parse_glycan_formula('HexNAc2.2Hex3.9Neu1')
            {'HexNAc': 2.2, 'Hex': 3.9, 'Neu': 1}

            # Negative counts
            >>> parse_glycan_formula('HexNAc-2Hex3Neu-1')
            {'HexNAc': -2, 'Hex': 3, 'Neu': -1}

            # Empty string
            >>> parse_glycan_formula('')
            {}

            # Unknown Glycan
            >>> parse_glycan_formula('HexXX')
            Traceback (most recent call last):
            peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "HexXX". Unknown glycan: "XX"!

            # Invalid Count
            >>> parse_glycan_formula('Hex2.1.')
            Traceback (most recent call last):
            peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "Hex2.1.". Invalid count: "2.1."!

            # Using a custom separator
            >>> parse_glycan_formula('HexNAc 2 Hex 3 Neu5Gc 1', sep=' ')
            {'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1}

    """

    try:
        return _parse_glycan_formula(formula, sep)
    except InvalidGlycanFormulaError as err:
        raise InvalidGlycanFormulaError(formula, err.msg) from err


def convert_glycan_formula_to_chem_formula(glycan: Union[ChemComposition, str]) -> str:
    """
    Converts a glycan string or dictionary to a chemical formula.

    :param glycan: A dictionary containing the glycan components and their counts, or a glycan formula string.
    :type glycan: ChemComposition

    :raises InvalidGlycanFormulaError: If the glycan formula is invalid.

    :return: A chemical formula string.
    :rtype: str

    .. code-block:: python

            # Int Counts
            >>> convert_glycan_formula_to_chem_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            'C45H73N3O34'

            # Float counts
            >>> convert_glycan_formula_to_chem_formula({'HexNAc': 2.2, 'Hex': 3, 'Neu5Gc': 1})
            'C46.6H75.6N3.2O35.0'

            # Negative counts
            >>> convert_glycan_formula_to_chem_formula({'HexNAc': -2.1, 'Hex': 3, 'Neu5Gc': -1})
            'C-9.8H-14.3N-3.1O-4.5'

            # Using a string
            >>> convert_glycan_formula_to_chem_formula('HexNAc2Hex3Neu5Gc1')
            'C45H73N3O34'

            # Unknown glycan
            >>> convert_glycan_formula_to_chem_formula('HexXX')
            Traceback (most recent call last):
            peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "HexXX". Unknown glycan: "XX"!

            # Invalid Count
            >>> convert_glycan_formula_to_chem_formula('Hex2.1.')
            Traceback (most recent call last):
            peptacular.errors.InvalidGlycanFormulaError: Error parsing glycan formula: "Hex2.1.". Invalid count: "2.1."!

    """

    try:
        return write_chem_formula(glycan_comp(glycan))
    except InvalidGlycanFormulaError as err:
        raise InvalidGlycanFormulaError(glycan, err.msg) from err
