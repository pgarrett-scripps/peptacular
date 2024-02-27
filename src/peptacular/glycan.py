from collections import Counter
from typing import Dict, Union

from peptacular.chem import parse_chem_formula, write_chem_formula
from peptacular.constants import MONOSACCHARIDES_NAMES_SORTED, MONOSACCHARIDE_ID_TO_COMPOSITIONS, \
    MONOSACCHARIDE_NAME_TO_ID


def parse_glycan_formula(formula: str) -> Dict[str, int]:
    """
    Parses a glycan sequence into its constituent parts.

    :param formula: A glycan formula string.
    :type formula: str

    :return: A dictionary containing the glycan components and their counts.
    :rtype: dict

    .. code-block:: python

            >>> parse_glycan_formula('HexNAc2Hex3Neu1')
            {'HexNAc': 2, 'Hex': 3, 'Neu': 1}

            >>> parse_glycan_formula('HexNAc2Hex3Neu5Gc1')
            {'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1}

            >>> parse_glycan_formula('HexNAc2Hex3Neu5Gc-1')
            {'HexNAc': 2, 'Hex': 3, 'Neu5Gc': -1}

    """

    d = {}

    while formula != '':
        for glycan_name in MONOSACCHARIDES_NAMES_SORTED:
            if formula.startswith(glycan_name):
                formula = formula[len(glycan_name):]
                count = ''

                # get the count (can have +- or digits)
                for c in formula:
                    if c.isdigit() or c in '+-':
                        count += c
                    else:
                        break

                # remove the count from the formula
                formula = formula[len(count):]

                # add the count to the dictionary
                d[glycan_name] = int(count)

                break
        else:
            raise ValueError(f'Could not parse glycan: {formula}')

    return d


def write_glycan_formula(glycan_dict: Dict[str, int]) -> str:
    """
    Writes a glycan dictionary to a string.

    :param glycan_dict: A dictionary containing the glycan components and their counts.
    :type glycan_dict: dict

    :return: A glycan formula string.
    :rtype: str

    .. code-block:: python

            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu': 1})
            'HexNAc2Hex3Neu1'

            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            'HexNAc2Hex3Neu5Gc1'

            >>> write_glycan_formula({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': -1})
            'HexNAc2Hex3Neu5Gc-1'

    """

    # Convert the dictionary to a list of strings
    glycan_str = [f'{component}{count}' for component, count in glycan_dict.items()]

    # Join the list of strings into a single string
    return ''.join(glycan_str)


def convert_glycan_formula_to_chem_formula(glycan: Union[Dict[str, int], str]) -> str:
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

    if isinstance(glycan, str):
        glycan = parse_glycan_formula(glycan)

    counts = Counter()
    for component, count in glycan.items():
        monosaccharide_id = MONOSACCHARIDE_NAME_TO_ID[component]
        chem_formula = parse_chem_formula(MONOSACCHARIDE_ID_TO_COMPOSITIONS[monosaccharide_id])
        for element, element_count in chem_formula.items():
            counts[element] += element_count * count

    return write_chem_formula(counts)
