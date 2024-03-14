"""
chem.py contains functions for parsing and writing chemical formulas, and for calculating the composition of a sequence
with/and modifications.
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Dict, Union, List, Optional

from peptacular.sequence.sequence import parse_single_sequence
from peptacular import constants
from peptacular.constants import AA_COMPOSITIONS, ION_TYPE_COMPOSITION_ADJUSTMENTS, \
    MONOSACCHARIDE_ID_TO_COMPOSITIONS, MONOSACCHARIDE_NAME_TO_ID
from peptacular.errors import InvalidFormulaError, InvalidCompositionError, AmbiguousAminoAcidError, \
    UnknownAminoAcidError, DeltaMassCompositionError
from peptacular.glycan import parse_glycan_formula
from peptacular.mod_db import parse_unimod_comp, parse_psi_comp, is_psi_mod_str, is_unimod_str, parse_xlmod_comp, \
    parse_resid_comp, is_resid_str, is_xlmod_str, is_gno_str, parse_gno_comp
from peptacular.sequence.proforma import parse_static_mods, ProFormaAnnotation, Mod, parse_isotope_mods
from peptacular.types import Chem_Composition
from peptacular.util import convert_type


def mod_comp(mod: str | Mod) -> Dict[str, int]:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str

    :raises InvalidCompositionError: If the modification string is invalid.

    :return: The parsed composition.
    :rtype: dict

    .. code-block:: python

        >>> mod_comp('Acetyl|INFO:newly discovered')
        {'H': 2, 'C': 2, 'O': 1}

        >>> mod_comp('Acetyl|Obs:+42.010565')
        {'H': 2, 'C': 2, 'O': 1}

    """

    if isinstance(mod, Mod):
        return mod_comp(mod.val)

    if isinstance(mod, (float, int)):
        raise InvalidCompositionError(mod)

    mods = mod.split('|')
    for m in mods:
        m = _parse_modification_composition(m)
        if m is not None:
            return m

    raise InvalidCompositionError(mod)


def estimate_comp(neutral_mass: float,
                  isotopic_mods: Optional[List[str] | List[Mod]] = None) -> Dict[str, float]:
    """
    Estimate the number of each element in a molecule based on its molecular mass using the averagine model.

    :param neutral_mass: The total neutral mass of the molecule.
    :type neutral_mass: float
    :param isotopic_mods: The isotopic modifications.
    :type isotopic_mods: list

    :return: The estimated number of each element in the molecule.
    :rtype: Dict[str, float]

    .. python::

        # Example usage
        >>> estimate_comp(1000)['C']
        44.468334554796016

        >>> estimate_comp(1000, ['13C'])['13C']
        44.468334554796016

    """

    composition = {atom: ratio * neutral_mass / constants.ISOTOPIC_AVERAGINE_MASS for atom, ratio in
                   constants.AVERAGINE_RATIOS.items()}

    if isotopic_mods:
        composition = _apply_isotope_mods_to_composition(composition, isotopic_mods)

    return composition


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

        >>> parse_chem_formula('[D]6H12O6666')
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
        formula = dict(sorted(formula.items(), key=lambda item: constants.HILL_ORDER[item[0]]))

    if sep is not None:
        return sep.join([f'{k}{sep}{v}' for k, v in formula.items() if v != 0])

    s = ''
    for k, v in formula.items():
        if v == 0:
            continue
        if k[0].isdigit() or k == 'D' or k == 'T':
            s += f'[{k}{v}]'
        else:
            s += f'{k}{v}'

    return s


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

def glycan_comp(glycan: Dict[str, int] | str) -> Chem_Composition:
    """
    Converts a glycan dictionary to a chemical formula.

    :param glycan: A dictionary containing the glycan components and their counts, or a glycan formula string.
    :type glycan: dict | str

    :return: A chemical formula string.
    :rtype: str

    .. code-block:: python

            >>> glycan_comp({'HexNAc': 2, 'Hex': 3, 'Neu5Gc': 1})
            {'C': 45, 'H': 73, 'N': 3, 'O': 34}

            >>> glycan_comp({'HexNAc': 2.2, 'Hex': 3, 'Neu5Gc': 1})
            {'C': 46.6, 'H': 75.6, 'N': 3.2, 'O': 35.0}

    """

    if isinstance(glycan, str):
        glycan = parse_glycan_formula(glycan)  # raises UnknownGlycanError

    counts = Counter()
    for component, count in glycan.items():

        monosaccharide_id = MONOSACCHARIDE_NAME_TO_ID[component]
        chem_formula = parse_chem_formula(MONOSACCHARIDE_ID_TO_COMPOSITIONS[monosaccharide_id])
        for element, element_count in chem_formula.items():
            counts[element] += element_count * count

    return dict(counts)


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


def _parse_glycan_comp_from_proforma_str(glycan_str: str) -> Dict[str, int]:
    """
    Parse a glycan string and return its mass.

    :param glycan_str: The glycan string to parse.
    :type glycan_str: str

    :raises UnknownGlycanError: If the glycan string contains an unknown monosaccharide.
    :raises InvalidFormulaError: If the formula is invalid.

    :return: The mass of the glycan.
    :rtype: float

    .. code-block:: python

        #  Get Composition
        >>> _parse_glycan_comp_from_proforma_str('HexNAc2Hex3Neu1')
        {'C': 43, 'H': 71, 'N': 3, 'O': 32}

        # Using a glycan name
        >>> _parse_glycan_comp_from_proforma_str('HexNAc')
        {'C': 8, 'H': 13, 'N': 1, 'O': 5}

        # Using a glycan ID
        >>> _parse_glycan_comp_from_proforma_str('6BAAE1B1')
        {'C': 3, 'H': 4, 'O': 2}

        # Using a glycan ID
        >>> _parse_glycan_comp_from_proforma_str('Glycan:6BAAE1B1')
        {'C': 3, 'H': 4, 'O': 2}

        # Invalid glycan ID
        >>> _parse_glycan_comp_from_proforma_str('6BAAXE121B1')
        Traceback (most recent call last):
        peptacular.errors.UnknownGlycanError: Unknown glycan: 6BAAXE121B1

        # Invalid glycan str
        >>> _parse_glycan_comp_from_proforma_str('HeSNAc2Hex3Neu1')
        Traceback (most recent call last):
        peptacular.errors.UnknownGlycanError: Unknown glycan: HeSNAc2Hex3Neu1

    """

    if glycan_str.lower().startswith('glycan:'):
        glycan_str = ''.join(glycan_str.split(':')[1:])

    glycan_id = constants.MONOSACCHARIDE_NAME_TO_ID.get(glycan_str, glycan_str)  # get id if possible
    glycan_comp = constants.MONOSACCHARIDE_ID_TO_COMPOSITIONS.get(glycan_id, None)  # raises UnknownGlycanError

    if glycan_comp is not None:
        return parse_chem_formula(glycan_comp)
    else:
        return parse_chem_formula(convert_glycan_formula_to_chem_formula(glycan_str))


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


def _parse_modification_composition(mod: str) -> Union[None, Dict[str, int]]:
    """
    Parses a modification composition and returns a dictionary with the element and their counts.

    :param mod: The modification composition.
    :type mod: str

    :raises InvalidCompositionError: If the modification composition string is invalid.
    :raises UnknownGlycanError: If the glycan formula contains an unknown glycan.
    :raises InvalidFormulaError: If the formula is invalid.
    :raises UnknownModificationError: If the modification is unknown.
    :raises NotImplementedError: If the modification is not implemented.

    :return: A dictionary with the element and their counts.
    :rtype: dict

    .. code-block:: python

        # Calculate the mass of a peptide sequence.
        >>> _parse_modification_composition('U:2')
        {'H': 1, 'N': 1, 'O': -1}

        >>> _parse_modification_composition('Formula:[13C4]H12')
        {'13C': 4, 'H': 12}

        >>> _parse_modification_composition('Glycan:HexNAc2Hex3Neu5Gc1')
        {'C': 45, 'H': 73, 'N': 3, 'O': 34}

        >>> _parse_modification_composition('1')

    """

    if isinstance(mod, (int, float)):  # cannot get composition from a delta mass
        return None

    converted_mod = convert_type(mod)

    if isinstance(converted_mod, (int, float)):  # cannot get composition from a delta mass
        return None

    # localization fix
    if isinstance(mod, str) and '#' in mod:
        if mod.startswith('#'):  # for localized positions, return empty comp
            return {}
        else:  # for only the original declaration consider the mass modification
            mod = mod.split('#')[0]

    mod_lower = mod.lower()
    if mod_lower.startswith('glycan:'):
        return _parse_glycan_comp_from_proforma_str(mod)

    elif is_gno_str(mod):
        # not implemented
        return parse_chem_formula(parse_gno_comp(mod))

    elif is_xlmod_str(mod):
        # not implemented
        return parse_chem_formula(parse_xlmod_comp(mod))

    elif is_resid_str(mod):
        # not implemented
        return parse_chem_formula(parse_resid_comp(mod))

    elif mod_lower.startswith('info:'):  # cannot get comp for info
        pass

    elif mod_lower.startswith('obs:'):  # cannot get comp for observed mass
        pass

    elif is_psi_mod_str(mod):  # psi-mod
        return parse_chem_formula(parse_psi_comp(mod))

    elif is_unimod_str(mod):  # unimod
        return parse_chem_formula(parse_unimod_comp(mod))

    # chemical formula
    elif mod_lower.startswith('formula:'):
        return _parse_chem_comp_from_proforma_str(mod)


def _parse_delta_mass(delta_mass: str) -> float:
    """
    Parse the delta mass.

    :param delta_mass: The delta mass.
    :type delta_mass: str
    :return: The parsed delta mass.
    :rtype: float

    .. code-block:: python

        >>> _parse_delta_mass('42.0')
        42.0

        >>> _parse_delta_mass('Acetyl')

        >>> _parse_delta_mass('UniMod:1')

        >>> _parse_delta_mass('UniMod:+1')
        1.0

        >>> _parse_delta_mass('Obs:42.0')
        42.0
    """

    if isinstance(delta_mass, (int, float)):
        return delta_mass

    mass = None
    for mod_str in delta_mass.split('|'):

        # check if the mod contains a localization score
        if '#' in mod_str:
            mod_str = mod_str.split('#')[0]

        try:
            mass = float(mod_str)
            break
        except ValueError:
            pass

        mod_lower = mod_str.lower()
        if mod_lower.startswith('unimod:') or mod_lower.startswith('u:') or \
                mod_lower.startswith('psi-mod:') or mod_lower.startswith('mod:') or mod_lower.startswith('m:') or \
                mod_lower.startswith('resid:') or mod_lower.startswith('r:') or \
                mod_lower.startswith('xlmod:') or mod_lower.startswith('x:') or \
                mod_lower.startswith('gno:') or mod_lower.startswith('g:'):

            mod_str = mod_str.split(':')[1]

            if mod_str.startswith('+') or mod_str.startswith('-'):
                mass = float(mod_str)
                break

        if mod_str.lower().startswith('obs:'):
            mod_str = mod_str.split(':')[1]
            try:
                mass = float(mod_str)
                break
            except ValueError:
                pass

    return mass


def _get_sequence_comp(sequence: str | ProFormaAnnotation, ion_type: str) -> Dict[str, int]:
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
        >>> _get_sequence_comp('PEPTIDE', 'y')
        {'C': 34, 'H': 53, 'N': 7, 'O': 15}

        >>> _get_sequence_comp('PEPTIDE', 'b')
        {'C': 34, 'H': 51, 'N': 7, 'O': 14}

        >>> _get_sequence_comp('<H>PEPTIDE', 'b')
        {'C': 34, 'H': 51, 'N': 7, 'O': 14}

        >>> _get_sequence_comp('{Unimod:2}PEPTIDE', 'p')
        {'C': 34, 'H': 54, 'N': 8, 'O': 14}

        >>> _get_sequence_comp('<13C>PEPTIDE', 'b')
        {'H': 51, 'N': 7, 'O': 14, '13C': 34}

        >>> _get_sequence_comp('PEPTIDE[Unimod:2]', 'y')
        {'C': 34, 'H': 54, 'N': 8, 'O': 14}

        >>> _get_sequence_comp('<[Unimod:2]@T>PEPTIDE', 'y')
        {'C': 34, 'H': 54, 'N': 8, 'O': 14}

        >>> _get_sequence_comp('<13C>PEPTIDE[Unimod:213413]', 'b')
        Traceback (most recent call last):
        peptacular.errors.UnknownModificationError: Unknown modification: Unimod:213413

        >>> _get_sequence_comp('I', 'by')
        {'C': 6, 'H': 12, 'N': 1, 'O': 1}

        # Ambiguous amino acid
        >>> _get_sequence_comp('B', 'by')
        Traceback (most recent call last):
        peptacular.errors.AmbiguousAminoAcidError: Ambiguous amino acid: B

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    if ion_type == 'p':
        ion_type = 'y'
    else:
        _ = annotation.pop_labile_mods()

    if 'B' in annotation.sequence:
        raise AmbiguousAminoAcidError('B')

    if 'J' in annotation.sequence:
        raise AmbiguousAminoAcidError('J')

    # Get the composition of the base sequence
    composition = {}
    for aa in annotation.sequence:
        try:
            aa_comp = AA_COMPOSITIONS[aa]
        except KeyError:
            raise UnknownAminoAcidError(aa)
        for k, v in aa_comp.items():
            composition[k] = composition.get(k, 0) + v

    # Apply the adjustments for the ion type
    for k, v in ION_TYPE_COMPOSITION_ADJUSTMENTS[ion_type].items():
        composition[k] = composition.get(k, 0) + v

    # Pop isotopic mods
    isotopic_mods = annotation.pop_isotope_mods()
    static_mods = annotation.pop_static_mods()

    mod_compositions = [composition]

    # Loop over unknown mods and apply them
    if annotation.has_unknown_mods():
        for unknown_mod in annotation.unknown_mods:
            mod_compositions.append(mod_comp(unknown_mod))

    # loop over intervals and apply them
    if annotation.has_intervals():
        for interval in annotation.intervals:
            if interval.has_mods():
                for interval_mod in interval.mods:
                    mod_compositions.append(mod_comp(interval_mod))

    # apply labile mods
    if annotation.has_labile_mods():
        for labile_mod in annotation.labile_mods:
            mod_compositions.append(mod_comp(labile_mod))

    # apply N-term mods
    if annotation.has_nterm_mods():
        for nterm_mod in annotation.nterm_mods:
            mod_compositions.append(mod_comp(nterm_mod))

    # apply C-term mods
    if annotation.has_cterm_mods():
        for cterm_mod in annotation.cterm_mods:
            mod_compositions.append(mod_comp(cterm_mod))

    # apply internal mods
    if annotation.has_internal_mods():
        for k, internal_mods in annotation.internal_mods.items():
            for internal_mod in internal_mods:
                mod_compositions.append(mod_comp(internal_mod))

    # Apply static mods
    if static_mods:
        static_map = parse_static_mods(static_mods)

        for aa, mod in static_map.items():
            aa_count = annotation.sequence.count(aa)
            for m in mod:
                mod_composition = mod_comp(m.val)

                mod_composition = {k: v * aa_count for k, v in mod_composition.items()}
                mod_compositions.append(mod_composition)

    # Merge the compositions
    composition = {}
    for comp in mod_compositions:
        for k, v in comp.items():
            composition[k] = composition.get(k, 0) + v

    # Apply isotopic mods
    if isotopic_mods:
        composition = _apply_isotope_mods_to_composition(composition, isotopic_mods)

    return composition


def _parse_mod_delta_mass_only(mod: str | Mod) -> Union[float, None]:
    """
    Parse a modification string.

    :param mod: The modification string.
    :type mod: str

    :raises ValueError: If the modification string is invalid.

    :return: The parsed composition.
    :rtype: dict

    .. code-block:: python

        >>> _parse_mod_delta_mass_only('Acetyl|INFO:newly discovered')

        >>> _parse_mod_delta_mass_only('Acetyl|Obs:+42.010565')

        >>> _parse_mod_delta_mass_only('Obs:+42.010565')
        42.010565

        >>> _parse_mod_delta_mass_only('R:+1|INFO:Fantastic')
        1.0

    """

    if isinstance(mod, float) or isinstance(mod, int):
        return mod

    if isinstance(mod, Mod):
        return _parse_mod_delta_mass_only(mod.val)

    mods = mod.split('|')
    for m in mods:
        try:
            m = _parse_modification_composition(m)
        except DeltaMassCompositionError:
            continue
        if m is not None:
            return None

    for m in mods:
        m = _parse_delta_mass(m)
        if m is not None:
            return m

    raise ValueError(f'Invalid modification: {mod}')


def _apply_isotope_mods_to_composition(composition: Dict[str, int], isotopic_mods: List[Mod]) -> Dict[str, int]:
    """
    Apply isotopic modifications to a composition.

    :param composition: The composition.
    :type composition: dict
    :param isotopic_mods: The isotopic modifications.
    :type isotopic_mods: list

    :return: The modified composition.
    :rtype: dict

    .. code-block:: python

        # Example usage
        >>> _apply_isotope_mods_to_composition({'C': 6, 'H': 12, 'O': 6}, ['13C'])
        {'H': 12, 'O': 6, '13C': 6}

        >>> _apply_isotope_mods_to_composition({'C': 6, 'H': 12, 'O': 6}, ['13C', '15N'])
        {'H': 12, 'O': 6, '13C': 6}

    """

    isotope_map = parse_isotope_mods(isotopic_mods)

    for element, isotope_label in isotope_map.items():
        if element in composition:

            if element == isotope_label:
                continue

            # Check if the isotope label already exists in the composition
            if isotope_label in composition:
                # Add the count of the original element to the isotope label count
                composition[isotope_label] += composition[element]
            else:
                # If the isotope label doesn't exist, create it with the count of the original element
                composition[isotope_label] = composition[element]

            del composition[element]

    return composition
