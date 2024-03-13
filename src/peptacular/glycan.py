from peptacular.util import convert_type
from peptacular.types import Chem_Composition
from peptacular.errors import UnknownGlycanError
from peptacular.constants import MONOSACCHARIDES_NAMES_SORTED


def parse_glycan_formula(formula: str) -> Chem_Composition:
    """
    Parses a glycan sequence into its constituent parts.

    :param formula: A glycan formula string.
    :type formula: str

    :raises UnknownGlycanError: If the glycan formula contains an unknown glycan.

    :return: A dictionary containing the glycan components and their counts.
    :rtype: dict

    .. code-block:: python

            >>> parse_glycan_formula('HexNAc2Hex3Neu1')
            {'HexNAc': 2, 'Hex': 3, 'Neu': 1}

            >>> parse_glycan_formula('HexNAc2Hex3Neu5Gc1')
            {'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1}

            >>> parse_glycan_formula('HexNAc2Hex3Neu5Gc-1')
            {'HexNAc': 2, 'Hex': 3, 'Neu5Gc': -1}

            >>> parse_glycan_formula('HexNAc2.2Hex3.9Neu1')
            {'HexNAc': 2.2, 'Hex': 3.9, 'Neu': 1}

            >>> parse_glycan_formula('')
            {}

            # This will raise an UnknownGlycanError
            >>> parse_glycan_formula('HeSNAc2Hex3Neu5Gc1X')
            Traceback (most recent call last):
            peptacular.errors.UnknownGlycanError: Unknown glycan: HeSNAc2Hex3Neu5Gc1X

    """

    d = {}

    if formula == '':
        return d

    original_formula = formula

    while formula != '':
        for glycan_name in MONOSACCHARIDES_NAMES_SORTED:
            if formula.startswith(glycan_name):
                formula = formula[len(glycan_name):]
                count = ''

                # get the count (can have +- or digits)
                for c in formula:
                    if c.isdigit() or c in '+-.':
                        count += c
                    else:
                        break

                # remove the count from the formula
                formula = formula[len(count):]

                # add the count to the dictionary
                d[glycan_name] = convert_type(count)

                break
        else:
            raise UnknownGlycanError(original_formula)

    return d


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

