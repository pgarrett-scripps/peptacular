from __future__ import annotations

import re
from collections import Counter
from typing import Dict

from peptacular import constants
from peptacular.constants import HILL_ORDER
from peptacular.errors import InvalidFormulaError, UnknownElementError
from peptacular.types import Chem_Composition
from peptacular.util import convert_type


def parse_chem_formula(formula: str) -> Dict[str, int]:
    """
    Parses a chemical formula and returns a dictionary with the element and their counts.

    :param formula: The chemical formula.
    :type formula: str

    :raises InvalidFormulaError: If the formula is invalid.

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        >>> parse_chem_formula('C6H12O6666')
        {'C': 6, 'H': 12, 'O': 6666}

        >>> parse_chem_formula('[D6]H12O6666')
        {'D': 6, 'H': 12, 'O': 6666}

        >>> parse_chem_formula('C6H12O-6666')
        {'C': 6, 'H': 12, 'O': -6666}

        >>> parse_chem_formula('[13C6666]H12O-6')
        {'13C': 6666, 'H': 12, 'O': -6}

        >>> parse_chem_formula('13C 6 H 12 O 6')
        {'13C': 6, 'H': 12, 'O': 6}

        >>> parse_chem_formula('C 6 Ce H 12.2 O 6 D')
        {'C': 6, 'Ce': 1, 'H': 12.2, 'O': 6, 'D': 1}

        >>> parse_chem_formula('C 1 1H -1.2 2H 3 D')
        {'C': 1, '1H': -1.2, '2H': 3, 'D': 1}

        # floats
        >>> parse_chem_formula('C4.45N8.22H59.99[13C34]S0.04O16.33')
        {'C': 4.45, 'N': 8.22, 'H': 59.99, '13C': 34, 'S': 0.04, 'O': 16.33}

        # floats
        >>> parse_chem_formula('C6H12O6ssss')
        Traceback (most recent call last):
        peptacular.errors.InvalidFormulaError: Cannot parse formula: C6H12O6ssss

    """

    if ' ' in formula:
        return _parse_split_chem_formula(formula)

    counters = []
    for component in _split_chem_formula(formula):
        if component.startswith('['):  # Isotope notation
            counters.append(_parse_isotope_component(component[1:-1]))
        else:
            counters.append(_parse_condensed_chem_formula(component))

    result_counter = Counter()
    for counter in counters:
        for key, value in counter.items():
            result_counter[key] += value

    return dict(result_counter)


def write_chem_formula(formula: Dict[str, int], sep: str = None, hill_order: bool = False) -> str:
    """
    Writes a chemical formula from a dictionary with the element and their counts.

    :param formula: A dictionary with the element and their counts.
    :type formula: dict
    :param sep: The separator to use between element and counts. Default is ''.
    :type sep: str
    :param hill_order: Whether to use Hill notation. Default is False.
    :type hill_order: bool

    :return: The chemical formula.
    :rtype: str

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> write_chem_formula({'C': 6, 'H': 12, 'O': 6})
        'C6H12O6'

        >>> write_chem_formula({'C': 6, 'H': 12, 'O': 6, 'S': 0})
        'C6H12O6'

        >>> write_chem_formula({'13C': 6, 'H': 12, 'O': 6})
        '[13C6]H12O6'

        >>> write_chem_formula({'13C': 6, 'H': 12, 'O': 6}, sep=' ')
        '13C 6 H 12 O 6'

        >>> write_chem_formula({'13C': 6, 'H': 12, 'O': 6}, sep='|')
        '13C|6|H|12|O|6'

        >>> write_chem_formula({'O': 6, 'H': 12, '13C': 6, 'C': 10}, hill_order=True)
        'C10[13C6]H12O6'

    """

    if hill_order:
        formula = dict(sorted(formula.items(), key=lambda item: HILL_ORDER[item[0]]))

    if sep is not None:
        return sep.join([f'{k}{sep}{v}' for k, v in formula.items() if v != 0])

    s = ''
    for k, v in formula.items():
        if v == 0 or k == '':
            continue
        if k[0].isdigit() or k == 'D' or k == 'T':
            s += f'[{k}{v}]'
        else:
            s += f'{k}{v}'

    return s


def _split_chem_formula(formula: str) -> list[str]:
    """
    Splits a chemical formula into its components. The formula is assumed to be proForma2.0 compliant, wherase the
    isotope notation for an element and its count is enclosed in square brackets and the formula contains no
    whitespace.

    :param formula: The chemical formula.
    :type formula: str

    :return: The components of the chemical formula.
    :rtype: list[str]

    .. code-block:: python

        # Split a chemical formula into its components.
        >>> _split_chem_formula('C6H12O6')
        ['C6H12O6']

        >>> _split_chem_formula('C6H12O-6')
        ['C6H12O-6']

        >>> _split_chem_formula('[13C6]H12O-6')
        ['[13C6]', 'H12O-6']

        >>> _split_chem_formula('[13C6]C6H12O6[13C6]')
        ['[13C6]', 'C6H12O6', '[13C6]']

    """
    components = []
    i = 0
    while i < len(formula):
        if formula[i] == '[':
            component_start = i
            component_end = formula.index(']', component_start)
            components.append(formula[component_start:component_end + 1])
            i = component_end + 1
        else:
            component_start = i
            while i < len(formula) and formula[i] not in '[]':
                i += 1
            components.append(formula[component_start:i])
    return components


def _parse_isotope_component(formula: str) -> Dict[str, int]:
    """
    Parses an isotope notation and returns the element and the mass number.

    :param formula: The isotope notation.
    :type formula: str

    :raises InvalidFormulaError: If the formula is invalid.

    :return: A tuple containing the element and the mass number.
    :rtype: tuple

    .. code-block:: python

        # Parse an isotope notation.
        >>> _parse_isotope_component('13C6')
        {'13C': 6}

        >>> _parse_isotope_component('D6')
        {'D': 6}

        >>> _parse_isotope_component('13C-6')
        {'13C': -6}

        >>> _parse_isotope_component('13C')
        {'13C': 1}

        >>> _parse_isotope_component('13Ce333')
        {'13Ce': 333}

        >>> _parse_isotope_component('13Ce1.3')
        {'13Ce': 1.3}

        >>> _parse_isotope_component('')
        {}

        >>> _parse_isotope_component('12323')
        Traceback (most recent call last):
        peptacular.errors.InvalidFormulaError: Cannot parse formula: 12323

        # example of invalid isotope notation
        >>> _parse_isotope_component('13C-')
        Traceback (most recent call last):
        peptacular.errors.InvalidFormulaError: Cannot parse formula: 13C-

    """
    counts = {}

    if formula == '':
        return counts

    if formula[0] == 'D' or formula[0] == 'T':
        return _parse_condensed_chem_formula(formula)

    match = re.match(r'([0-9]*)([A-Za-z]+)(-?\d*\.?\d*)', formula)
    if match:
        mass_number_and_element = match.group(1) + match.group(2)  # Combine mass number (if any) and element symbol
        count = convert_type(match.group(3)) if match.group(3) else 1  # Get the count, defaulting to 1 if not provided
        if isinstance(count, str):
            raise InvalidFormulaError(formula)
        counts.setdefault(mass_number_and_element, 0)
        counts[mass_number_and_element] += count
    else:
        raise InvalidFormulaError(formula)

    return counts


def _parse_condensed_chem_formula(formula: str) -> Dict[str, int]:
    """
    Parses a chemical formula and returns a dictionary with the element and their counts.

    :param formula: The chemical formula.
    :type formula: str

    :raises InvalidFormulaError: If the formula is invalid.

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> _parse_condensed_chem_formula('C6H12O6')
        {'C': 6, 'H': 12, 'O': 6}

        >>> _parse_condensed_chem_formula('C6H12O6C6')
        {'C': 12, 'H': 12, 'O': 6}

        >>> _parse_condensed_chem_formula('C6H12O-6')
        {'C': 6, 'H': 12, 'O': -6}

        >>> _parse_condensed_chem_formula('C6H12O6Ce333')
        {'C': 6, 'H': 12, 'O': 6, 'Ce': 333}

        >>> _parse_condensed_chem_formula('C6H12.1O6C6Ce333.1')
        {'C': 12, 'H': 12.1, 'O': 6, 'Ce': 333.1}

        >>> _parse_condensed_chem_formula('')
        {}

        # example of invalid formula
        >>> _parse_condensed_chem_formula('123')
        Traceback (most recent call last):
        peptacular.errors.InvalidFormulaError: Cannot parse formula: 123

        # example of invalid formula
        >>> _parse_condensed_chem_formula('C6H12O-')
        Traceback (most recent call last):
        peptacular.errors.InvalidFormulaError: Cannot parse formula: C6H12O-

        # example of invalid formula
        >>> _parse_condensed_chem_formula('C6H12O6ssss')
        Traceback (most recent call last):
        peptacular.errors.InvalidFormulaError: Cannot parse formula: C6H12O6ssss

    """

    element_counts = {}

    if formula == '':
        return element_counts

    # Find element and their counts in non-isotope components
    matches = list(re.finditer(r'([A-Z][a-z]*)(-?\d*\.?\d*)', formula))

    if len(matches) == 0:
        raise InvalidFormulaError(formula)

    matched_chars = 0
    for match in matches:
        element = match.group(1)
        count_str = match.group(2)
        matched_chars += len(element) + len(count_str)
        count = convert_type(count_str) if count_str else 1
        if isinstance(count, str):
            raise InvalidFormulaError(formula)
        element_counts.setdefault(element, 0)
        element_counts[element] += count

    if matched_chars != len(formula):
        raise InvalidFormulaError(formula)

    return element_counts


def _parse_split_chem_formula(formula: str) -> Dict[str, int]:
    """
    Parses a chemical formula and returns a dictionary with the element and their counts.

    :param formula: The chemical formula.
    :type formula: str

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        >>> _parse_split_chem_formula('C 6 H 12 O 6')
        {'C': 6, 'H': 12, 'O': 6}

        >>> _parse_split_chem_formula('C 1 1H -1 2H 3 D')
        {'C': 1, '1H': -1, '2H': 3, 'D': 1}

        >>> _parse_split_chem_formula('C 6 Ce H 12 O 6 T')
        {'C': 6, 'Ce': 1, 'H': 12, 'O': 6, 'T': 1}

        >>> _parse_split_chem_formula('C 6 Ce H 12.2 O 6 T')
        {'C': 6, 'Ce': 1, 'H': 12.2, 'O': 6, 'T': 1}

    """

    element_counts = {}
    components = formula.split(' ')

    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    i = 0
    while i < len(components):
        element = components[i]
        # Check if the next component is a number (positive or negative)
        if i + 1 < len(components) and is_number(components[i + 1]):
            count = convert_type(components[i + 1])
            i += 2  # Skip the next item since it's the count
        else:
            count = 1  # Default count is 1
            i += 1

        # Add or update the element count in the dictionary
        if element in element_counts:
            element_counts[element] += count
        else:
            element_counts[element] = count

    return element_counts


def _parse_chem_comp_from_proforma_str(chem_str: str) -> Dict[str, int]:
    """
    Parse a chemical formula string and return its composition.

    :param chem_str: The chemical formula string to parse.
    :type chem_str: str

    :raises InvalidFormulaError: If the formula is invalid.

    :return: The composition of the chemical formula.
    :rtype: dict

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> _parse_chem_comp_from_proforma_str('C6H12O6')
        {'C': 6, 'H': 12, 'O': 6}

        >>> _parse_chem_comp_from_proforma_str('C6H12O-6')
        {'C': 6, 'H': 12, 'O': -6}

        >>> _parse_chem_comp_from_proforma_str('[13C6]H12O-6')
        {'13C': 6, 'H': 12, 'O': -6}

        >>> _parse_chem_comp_from_proforma_str('13C 6 H 12 O 6')
        {'13C': 6, 'H': 12, 'O': 6}

    """

    if chem_str.lower().startswith('formula:'):
        chem_str = ''.join(chem_str.split(':')[1:])

    return parse_chem_formula(chem_str)


def chem_mass(formula: Chem_Composition | str, monoisotopic: bool = True, precision: int | None = None) -> float:
    """
    Calculate the mass of a chemical formula.

    :param formula: The chemical formula.
    :type formula: dict
    :param monoisotopic: Whether to use monoisotopic masses.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to.
    :type precision: int

    :raises UnknownElementError: If the chemical formula contains an unknown element.

    :return: The mass of the chemical formula.
    :rtype: float

    .. code-block:: python

        # Calculate the mass of a chemical formula.
        >>> chem_mass({'C': 6, 'H': 12, 'O': 6}, precision=3)
        180.063

        >>> chem_mass({'13C': 6, 'H': 12, 'O': 6}, precision=3)
        186.084

        >>> chem_mass({'C': 6, 'H': 12, 'O': 6}, monoisotopic=False, precision=3)
        180.156

        # Use average masses for all elements except for iosotopes.
        >>> chem_mass({'13C': 6, 'H': 12, 'O': 6}, monoisotopic=False, precision=3)
        186.112

        # Use average masses for all elements except for iosotopes.
        >>> chem_mass({'13C': 6, 'D': 12, 'O': 6}, monoisotopic=False, precision=3)
        198.186

        # Example Error
        >>> chem_mass({'C': 6, 'H': 12, 'O': 6, 'X': 1}, precision=3)
        Traceback (most recent call last):
        peptacular.errors.UnknownElementError: Unknown element: X

    """

    if isinstance(formula, str):
        formula = parse_chem_formula(formula)

    m = 0.0
    for element, count in formula.items():

        if element not in constants.ISOTOPIC_ATOMIC_MASSES:
            raise UnknownElementError(element)

        if monoisotopic is True:
            m += constants.ISOTOPIC_ATOMIC_MASSES[element] * count
        else:
            if element[0].isdigit() or element == 'D' or element == 'T':  # element is a isotope
                m += constants.ISOTOPIC_ATOMIC_MASSES[element] * count
            else:
                m += constants.AVERAGE_ATOMIC_MASSES[element] * count

    if precision is not None:
        m = round(m, precision)

    return m
