"""
chem_util.py
"""

from typing import Union, Optional, List

from peptacular.constants import HILL_ORDER, ISOTOPIC_ATOMIC_MASSES, ELECTRON_MASS, PROTON_MASS, \
    NEUTRON_MASS, AVERAGE_ATOMIC_MASSES, ISOTOPE_COMPONENT_PATTERN, CONDENSED_CHEM_FORMULA_PATTERN
from peptacular.errors import InvalidChemFormulaError
from peptacular.util import convert_type

from peptacular.types import ChemComposition


def parse_chem_formula(formula: str, sep: str = '') -> ChemComposition:
    """
    Parses a chemical formula and returns a dict mapping elements to their counts.

    :param formula: The chemical formula.
    :type formula: str
    :param sep: The separator to use between element and counts. Default is ''.
    :type sep: str

    :raises InvalidFormulaError: If the formula is invalid.

    :return: A dictionary with the element and their counts.
    :rtype: ChemComposition

    .. code-block:: python

        >>> parse_chem_formula('C6H12O6666')
        {'C': 6, 'H': 12, 'O': 6666}

        >>> parse_chem_formula('[D6]H12O6666')
        {'D': 6, 'H': 12, 'O': 6666}

        >>> parse_chem_formula('C6H12O-6666')
        {'C': 6, 'H': 12, 'O': -6666}

        >>> parse_chem_formula('[13C6666]H12O-6')
        {'13C': 6666, 'H': 12, 'O': -6}

        >>> parse_chem_formula('13C 6 H 12 O 6', sep=' ')
        {'13C': 6, 'H': 12, 'O': 6}

        >>> parse_chem_formula('C 6 Ce H 12.2 O 6 D', sep=' ')
        {'C': 6, 'Ce': 1, 'H': 12.2, 'O': 6, 'D': 1}

        >>> parse_chem_formula('C 1 1H -1.2 2H 3 D', sep=' ')
        {'C': 1, '1H': -1.2, '2H': 3, 'D': 1}

        # Example using floats
        >>> parse_chem_formula('C4.45N8.22H59.99[13C34]S0.04O16.33')
        {'C': 4.45, 'N': 8.22, 'H': 59.99, '13C': 34, 'S': 0.04, 'O': 16.33}

        # Example using floats and electron
        >>> parse_chem_formula('C4.45N8.22H59.99[13C34]S0.04e-16.33')
        {'C': 4.45, 'N': 8.22, 'H': 59.99, '13C': 34, 'S': 0.04, 'e': -16.33}

        # Example Error with invalid formula
        >>> parse_chem_formula('C6H12O6ssss')
        Traceback (most recent call last):
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "C6H12O6ssss". Cannot parse: "ssss"!

    """

    if sep != '':
        return _parse_split_chem_formula(formula, sep)

    comps = []
    for component in _split_chem_formula(formula):
        if component.startswith('['):  # Isotope notation
            try:
                comps.append(_parse_isotope_component(component[1:-1]))
            except InvalidChemFormulaError as err:
                raise InvalidChemFormulaError(formula, err.msg) from err
        else:
            try:
                comps.append(_parse_condensed_chem_formula(component))
            except InvalidChemFormulaError as err:
                raise InvalidChemFormulaError(formula, err.msg) from err

    combined_comp = {}
    for comp in comps:
        for key, value in comp.items():
            combined_comp[key] = combined_comp.get(key, 0) + value

    return combined_comp


def write_chem_formula(composition: ChemComposition, sep: str = '', hill_order: bool = False) -> str:
    """
    Writes a chemical formula from a dict mapping elements to their counts.

    :param composition: A dictionary with the element and their counts.
    :type: composition: ChemComposition
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

        # formula containg floats and electron
        >>> write_chem_formula({'C': 4.45, 'N': 8.22, 'H': 59.99, '13C': 34, 'S': 0.04, 'e': -16.33}, hill_order=True)
        'C4.45[13C34]H59.99N8.22S0.04e-16.33'

    """

    if hill_order:
        composition = dict(sorted(composition.items(), key=lambda item: HILL_ORDER.get(item[0], 10_000)))

    if sep != '':
        return sep.join([f'{k}{sep}{v}' for k, v in composition.items() if v != 0])

    s = ''
    for k, v in composition.items():
        if v == 0 or k == '':
            continue
        if k[0].isdigit() or k == 'D' or k == 'T':
            s += f'[{k}{v}]'
        else:
            s += f'{k}{v}'

    return s


# Must keep here to avoid circular dep with mod_db (since mod_db uses mass_calc, and mass_calc uses chem_util)
def chem_mass(formula: Union[ChemComposition, str],
              monoisotopic: bool = True,
              precision: Optional[int] = None,
              sep: str = '') -> float:
    """
    Calculate the mass of a chemical formula or composition.

    :param formula: The chemical formula or composition.
    :type formula: Union[ChemComposition, str]
    :param monoisotopic: Whether to use monoisotopic masses. Default is True.
    :type monoisotopic: bool
    :param precision: The number of decimal places to round the mass to. Default is None.
    :type precision: Optional[int]
    :param sep: The separator to use between element and counts. Default is ''.
    :type sep: str

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

        # Complex example using floats and particles
        >>> chem_mass('C4.45N8.22H59.99[13C34]S0.04e-16.33', monoisotopic=False, precision=3)
        672.437

        # Example Error
        >>> chem_mass("C6X2", precision=3)
        Traceback (most recent call last):
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "{'C': 6, 'X': 2}". Unknown element: "X"!

    """

    if isinstance(formula, str):
        formula = parse_chem_formula(formula, sep)

    m = 0.0
    for element, count in formula.items():

        if element not in ISOTOPIC_ATOMIC_MASSES:
            if element == 'e':
                m += ELECTRON_MASS * count
            elif element == 'p':
                m += PROTON_MASS * count
            elif element == 'n':
                m += NEUTRON_MASS * count
            else:
                raise InvalidChemFormulaError(formula, f'Unknown element: "{element}"!')

            continue

        if monoisotopic is True:
            m += ISOTOPIC_ATOMIC_MASSES[element] * count
        else:
            if element[0].isdigit() or element == 'D' or element == 'T':  # element is a isotope
                m += ISOTOPIC_ATOMIC_MASSES[element] * count
            else:
                m += AVERAGE_ATOMIC_MASSES[element] * count

    if precision is not None:
        m = round(m, precision)

    return m



def _split_chem_formula(formula: str) -> List[str]:
    """
    Splits a chemical formula into its components. The formula is assumed to be proForma2.0 compliant, wherase the
    isotope notation for an element and its count is enclosed in square brackets and the formula contains no
    whitespace.

    :param formula: The chemical formula.
    :type formula: str

    :return: The components of the chemical formula, split by instances of the isotope notation.
    :rtype: List[str]

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


def _parse_isotope_component(formula: str) -> ChemComposition:
    """
    Parses an isotope notation and returns the element and the mass number.

    :param formula: The isotope notation.
    :type formula: str

    :raises InvalidFormulaError: If the formula is invalid.

    :return: A dictionary with the element and their counts.
    :rtype: ChemComposition

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

        >>> _parse_isotope_component('12323') # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "12323". ...

        # example of invalid isotope notation
        >>> _parse_isotope_component('13C-') # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "13C-". Invalid count: "-"!

    """
    counts = {}

    if formula == '':
        return counts

    if formula[0] == 'D' or formula[0] == 'T':
        return _parse_condensed_chem_formula(formula)

    match = ISOTOPE_COMPONENT_PATTERN.match(formula)
    if match:
        mass_number_and_element = match.group(1) + match.group(2)  # Combine mass number (if any) and element symbol
        count = convert_type(match.group(3)) if match.group(3) else 1  # Get the count, defaulting to 1 if not provided
        if isinstance(count, str):
            raise InvalidChemFormulaError(formula, f'Invalid count: "{count}"!')
        counts.setdefault(mass_number_and_element, 0)
        counts[mass_number_and_element] += count
    else:
        raise InvalidChemFormulaError(formula, f'Invalid isotope notation: "{formula}"!')

    return counts


def _parse_condensed_chem_formula(formula: str) -> ChemComposition:
    """
    Parses a chemical formula and returns a dict mapping elements to their counts.

    :param formula: The chemical formula.
    :type formula: str

    :raises InvalidFormulaError: If the formula is invalid.

    :return: A dictionary with the element and their counts.
    :rtype: ChemComposition

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

        # also works for electron, proton, and neutron
        >>> _parse_condensed_chem_formula('C6e-2n-2p2')
        {'C': 6, 'e': -2, 'n': -2, 'p': 2}

        # example of invalid formula
        >>> _parse_condensed_chem_formula('123')
        Traceback (most recent call last):
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "123". Cannot parse: "123"!

        # example of invalid formula
        >>> _parse_condensed_chem_formula('C6H12O-')
        Traceback (most recent call last):
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "C6H12O-". Invalid count: "-"!

        # example of invalid formula
        >>> _parse_condensed_chem_formula('C6H12O6ssss')
        Traceback (most recent call last):
        peptacular.errors.InvalidChemFormulaError: Error parsing chem formula: "C6H12O6ssss". Cannot parse: "ssss"!

    """

    element_counts = {}

    if formula == '':
        return element_counts

    # Find element and their counts in non-isotope components
    matches = list(CONDENSED_CHEM_FORMULA_PATTERN.finditer(formula))

    if len(matches) == 0:
        raise InvalidChemFormulaError(formula, f'Cannot parse: "{formula}"!')

    matched_chars = 0
    for match in matches:
        element = match.group(1)
        count_str = match.group(2)
        matched_chars += len(element) + len(count_str)
        count = convert_type(count_str) if count_str else 1
        if isinstance(count, str):
            raise InvalidChemFormulaError(formula, f'Invalid count: "{count}"!')
        element_counts.setdefault(element, 0)
        element_counts[element] += count

    if matched_chars != len(formula):
        raise InvalidChemFormulaError(formula, f'Cannot parse: "{formula[matched_chars:]}"!')

    return element_counts


def _parse_split_chem_formula(formula: str, sep: str) -> ChemComposition:
    """
    Parses a chemical formula and returns a dict mapping elements to their counts.

    :param formula: The chemical formula.
    :type formula: str
    :param sep: The separator to use between element and counts.
    :type sep: str

    :return: A dictionary with the element and their counts.
    :rtype: ChemComposition

    .. code-block:: python

        >>> _parse_split_chem_formula('C 6 H 12 O 6', ' ')
        {'C': 6, 'H': 12, 'O': 6}

        >>> _parse_split_chem_formula('C 1 1H -1 2H 3 D', ' ')
        {'C': 1, '1H': -1, '2H': 3, 'D': 1}

        >>> _parse_split_chem_formula('C 6 Ce H 12 O 6 T', ' ')
        {'C': 6, 'Ce': 1, 'H': 12, 'O': 6, 'T': 1}

        >>> _parse_split_chem_formula('C 6 Ce H 12.2 O 6 T', ' ')
        {'C': 6, 'Ce': 1, 'H': 12.2, 'O': 6, 'T': 1}

        >>> _parse_split_chem_formula('C 6 Ce H 12.2 O 6 D', ' ')
        {'C': 6, 'Ce': 1, 'H': 12.2, 'O': 6, 'D': 1}

        >>> _parse_split_chem_formula('e 13 n 12 p 6', ' ')
        {'e': 13, 'n': 12, 'p': 6}

    """

    element_counts = {}
    components = formula.split(sep)

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
