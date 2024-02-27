import re
from collections import Counter
from typing import Dict

from peptacular.constants import AA_COMPOSITIONS, ION_TYPE_COMPOSITION_ADJUSTMENTS


def _split_formula(formula: str) -> list[str]:
    """
    Splits a chemical formula into its components. The formula is assumed to be proForma2.0 compliant, wherase the
    isotpoe notation for an element and its count is enclosed in square brackets and the formula contains no
    whitespace.

    :param formula: The chemical formula.
    :type formula: str

    :return: The components of the chemical formula.
    :rtype: list[str]

    .. code-block:: python

        # Split a chemical formula into its components.
        >>> _split_formula('C6H12O6')
        ['C6H12O6']

        >>> _split_formula('C6H12O-6')
        ['C6H12O-6']

        >>> _split_formula('[13C6]H12O-6')
        ['[13C6]', 'H12O-6']

        >>> _split_formula('[13C6]C6H12O6[13C6]')
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


def _parse_isotope(formula: str) -> Counter[str, int]:
    """
    Parses an isotope notation and returns the element and the mass number.

    :param formula: The isotope notation.
    :type formula: str

    :return: A tuple containing the element and the mass number.
    :rtype: tuple

    .. code-block:: python

        # Parse an isotope notation.
        >>> _parse_isotope('13C6')
        Counter({'13C': 6})

        >>> _parse_isotope('13C-6')
        Counter({'13C': -6})

        >>> _parse_isotope('13C')
        Counter({'13C': 1})

        >>> _parse_isotope('13Ce333')
        Counter({'13Ce': 333})

    """
    counts = Counter()
    match = re.match(r'([0-9]*)([A-Za-z]+)(-?\d*)', formula)
    if match:
        mass_number_and_element = match.group(1) + match.group(2)  # Combine mass number (if any) and element symbol
        count = int(match.group(3)) if match.group(3) else 1  # Get the count, defaulting to 1 if not provided
        counts[mass_number_and_element] = count
    else:
        raise ValueError(f'Invalid isotope notation: {formula}')

    return counts


def _parse_formula(formula: str) -> Counter[str, int]:
    """
    Parses a chemical formula and returns a dictionary with the element and their counts.

    :param formula: The chemical formula.
    :type formula: str

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> _parse_formula('C6H12O6')
        Counter({'H': 12, 'C': 6, 'O': 6})

        >>> _parse_formula('C6H12O6C6')
        Counter({'C': 12, 'H': 12, 'O': 6})

        >>> _parse_formula('C6H12O-6')
        Counter({'H': 12, 'C': 6, 'O': -6})

        >>> _parse_formula('C6H12O6Ce333')
        Counter({'Ce': 333, 'H': 12, 'C': 6, 'O': 6})


    """
    element_counts = Counter()

    # Find element and their counts in non-isotope components
    matches = re.finditer(r'([A-Z][a-z]*)(-?\d*)', formula)
    for match in matches:
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        element_counts[element] += count

    return element_counts


def parse_chem_formula(formula: str) -> Dict[str, int]:
    """
    Parses a chemical formula and returns a dictionary with the element and their counts.

    :param formula: The chemical formula.
    :type formula: str

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> parse_chem_formula('C6H12O6')
        {'C': 6, 'H': 12, 'O': 6}

        >>> parse_chem_formula('C6H12O6C6')
        {'C': 12, 'H': 12, 'O': 6}

        >>> parse_chem_formula('C6H12O-6')
        {'C': 6, 'H': 12, 'O': -6}

        >>> parse_chem_formula('C6H12O6Ce333')
        {'C': 6, 'H': 12, 'O': 6, 'Ce': 333}

        >>> parse_chem_formula('[13C6]H12O-6')
        {'13C': 6, 'H': 12, 'O': -6}

        >>> parse_chem_formula('([13C6]H12O6')
        {'13C': 6, 'H': 12, 'O': 6}

        >>> parse_chem_formula('13C 6 H 12 O 6')
        {'13C': 6, 'H': 12, 'O': 6}

        >>> parse_chem_formula('[13Ce333]H12O6')
        {'13Ce': 333, 'H': 12, 'O': 6}

        >>> parse_chem_formula('C 1 1H -1 2H 3')
        {'C': 1, '1H': -1, '2H': 3}

    """

    element_counts = {}
    if ' ' in formula:
        components = formula.split(' ')
        for i in range(0, len(components), 2):
            element_counts[components[i]] = int(components[i + 1])
        return element_counts

    counters = []
    for component in _split_formula(formula):
        if component.startswith('['):  # Isotope notation
            counters.append(_parse_isotope(component[1:-1]))
        else:
            counters.append(_parse_formula(component))

    result_counter = Counter()
    for counter in counters:
        for key, value in counter.items():
            result_counter[key] += value

    return dict(result_counter)


def write_chem_formula(formula: Dict[str, int], sep: str = '') -> str:
    """
    Writes a chemical formula from a dictionary with the element and their counts.

    :param formula: A dictionary with the element and their counts.
    :type formula: dict
    :param sep: The separator to use between element and counts. Default is ''.
    :type sep: str

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

    """

    if sep == ' ':
        return ' '.join([f'{k} {v}' for k, v in formula.items() if v != 0])

    s = ''
    for k, v in formula.items():
        if v == 0:
            continue
        if k[0].isdigit():
            s += f'[{k}{v}]'
        else:
            s += f'{k}{v}'

    return s


def calculate_sequence_composition(sequence: str, ion_type: str) -> Dict[str, int]:
    """
    Calculate the composition of a sequence.

    :param sequence: The sequence.
    :type sequence: str
    :param ion_type: The ion type.
    :type ion_type: str

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> calculate_sequence_composition('PEPTIDE', 'y')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> calculate_sequence_composition('PEPTIDE', 'b')
        {'C': 34, 'H': 51, 'N': 7, 'O': 14}

        >>> calculate_sequence_composition('I', 'by')
        {'C': 6, 'H': 12, 'N': 1, 'O': 1}

    """

    composition = {}
    for aa in sequence:
        aa_comp = AA_COMPOSITIONS[aa]
        for k, v in aa_comp.items():
            composition[k] = composition.get(k, 0) + v

    ion_type_adjustment = ION_TYPE_COMPOSITION_ADJUSTMENTS[ion_type]

    for k, v in ion_type_adjustment.items():
        composition[k] = composition.get(k, 0) + v

    return composition
