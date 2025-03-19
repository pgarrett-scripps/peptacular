"""
ProForma Parser Module
"""

import copy
import itertools
import random
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from functools import wraps
from typing import List, Dict, Optional, Generator, Any, Callable, Tuple, Union
from typing import Counter as CounterType
import regex as re

from peptacular.proforma.input_convert import fix_list_of_mods, fix_dict_of_mods, fix_intervals_input, \
    ACCEPTED_MOD_INPUT, ACCEPTED_INTERVAL_INPUT, ModValue, INTERVAL_VALUE
from peptacular.proforma.proforma_dataclasses import Mod, Interval, are_mods_equal, are_intervals_equal
from peptacular.constants import AMINO_ACIDS, AMBIGUOUS_AMINO_ACIDS, MASS_AMBIGUOUS_AMINO_ACIDS, ADDUCT_PATTERN, \
    ISOTOPE_NUM_PATTERN
from peptacular.errors import ProFormaFormatError
from peptacular.types import ChemComposition


def _parse_modifications(proforma_sequence: str, opening_bracket: str = '[', closing_bracket: str = ']') -> List[Mod]:
    """
    Parse modifications from a proforma sequence with support for nested brackets.

    :param proforma_sequence: The proforma sequence to parse.
    :type proforma_sequence: str
    :param opening_bracket: The opening bracket character.
    :type opening_bracket: str
    :param closing_bracket: The closing bracket character.
    :type closing_bracket: str

    :return: The parsed modifications.
    :rtype: List[Mod]

    .. code-block:: python

        >>> _parse_modifications('PEPTIDE[1]@C')
        [Mod(1, 1)]

        >>> _parse_modifications('PEPTIDE[Formula:[13C]H23]@S,D')
        [Mod('Formula:[13C]H23', 1)]

    """
    mods = []
    position = 0
    length = len(proforma_sequence)
    while position < length:
        if proforma_sequence[position] == opening_bracket:
            mod, new_position = _parse_modification(proforma_sequence[position:], opening_bracket, closing_bracket)
            mods.append(mod)
            position += new_position  # Update position based on the length of the parsed modification.
        else:
            position += 1  # Move to the next character if it's not an opening bracket.
    return mods


def _parse_modification(proforma_sequence: str,
                        opening_bracket: str = '[',
                        closing_bracket: str = ']') -> Tuple[Mod, int]:
    """
    Parse a single modification from a proforma sequence, handling nested brackets.
    Returns the parsed Mod and the position after the modification.

    :param proforma_sequence: The proforma sequence to parse.
    :type proforma_sequence: str
    :param opening_bracket: The opening bracket character. Default is '['.
    :type opening_bracket: str
    :param closing_bracket: The closing bracket character. Default is ']'.
    :type closing_bracket: str

    :return: The parsed Mod and the position after the modification.
    :rtype: Tuple[Mod, int]

    .. code-block:: python

        >>> _parse_modification('[1]@C')
        (Mod(1, 1), 3)

        >>> _parse_modification('[3.1415]@S,D')
        (Mod(3.1415, 1), 8)

        >>> _parse_modification('[3.1415]^+2@S,D')
        (Mod(3.1415, 2), 11)

    """
    bracket_depth = 1  # Start with a depth of 1 for the initial opening bracket.
    position = 1  # Start from the character after the opening bracket.
    mod_start = 1  # Modification starts after the opening bracket.

    # Find the matching closing bracket, accounting for nesting.
    while position < len(proforma_sequence) and bracket_depth > 0:
        if proforma_sequence[position] == opening_bracket:
            bracket_depth += 1
        elif proforma_sequence[position] == closing_bracket:
            bracket_depth -= 1
        position += 1

    # Extract the modification string.
    mod_end = position - 1  # Exclude the closing bracket.
    mod_str = proforma_sequence[mod_start:mod_end]

    multiplier = 1
    # Check for a multiplier immediately following the modification.
    if position < len(proforma_sequence) and proforma_sequence[position] == '^':
        multiplier_start = position + 1
        multiplier, parsed_length = _parse_integer(proforma_sequence[multiplier_start:])
        position += parsed_length + 1  # Account for the '^' and the length of the integer parsed.

    return Mod(mod_str, multiplier), position


def _parse_integer(proforma_sequence: str) -> Tuple[int, int]:
    """
    Parse an integer from the proforma sequence, returning the integer and the number of characters parsed.

    :param proforma_sequence: The proforma sequence to parse.
    :type proforma_sequence: str

    :return: The parsed integer and the number of characters parsed.
    :rtype: Tuple[int, int]

    .. code-block:: python

        >>> _parse_integer('123')
        (123, 3)

        >>> _parse_integer('+123')
        (123, 4)

        >>> _parse_integer('-123')
        (-123, 4)

        >>> _parse_integer('-123+PEPTIDE')
        (-123, 4)

    """
    i = 0
    digit_count = 0
    while i < len(proforma_sequence):
        if proforma_sequence[i].isdigit():
            digit_count += 1
            i += 1
        elif digit_count == 0 and proforma_sequence[i] in ['+', '-']:
            i += 1
        else:
            break

    return int(proforma_sequence[:i]), i


def parse_charge_adducts(mod: ModValue) -> ChemComposition:
    """
    Parse charge adducts into a dictionary, mapping the ion to its count
    :param mod: The charge adducts to parse.
    :type mod: ModValue

    :raises TypeError: If the mod is not a string or Mod instance.

    :return: Dictionary of charge adducts.
    :rtype: ChemComposition

    .. code-block:: python

        >>> parse_charge_adducts('+2Na+,+H+')
        {'Na+': 2, 'H+': 1}

        >>> parse_charge_adducts('+2Na+,-H+')
        {'Na+': 2, 'H+': -1}

        >>> parse_charge_adducts('2I-')
        {'I-': 2}

        >>> parse_charge_adducts('+e-')
        {'e-': 1}

        >>> parse_charge_adducts(Mod('-2Na+,+H+', 1))
        {'Na+': -2, 'H+': 1}

        >>> parse_charge_adducts(Mod('+2Mg2+,+H+', 1))
        {'Mg2+': 2, 'H+': 1}


    """

    if isinstance(mod, Mod):
        return parse_charge_adducts(mod.val)

    if not isinstance(mod, str):
        raise TypeError(f"Invalid type for charge adducts: {type(mod)}! Mod or Mod.val must be of type string.")

    charge_adducts = {}

    mods = mod.split(',')
    for m in mods:
        if isinstance(m, Mod):
            mod_str = m.val
        else:
            mod_str = m

        for sign, count, ion in ADDUCT_PATTERN.findall(mod_str):
            if count == '':
                count = 1
            else:
                count = int(count)

            if sign == '-':
                count = -count

            if ion in charge_adducts:
                charge_adducts[ion] += count
            else:
                charge_adducts[ion] = count

    return charge_adducts


def write_charge_adducts(charge_adducts: ChemComposition) -> Mod:
    """
    Converts the dictionary of charge adducts into a list of Mod instances.

    :param charge_adducts: Dictionary of charge adducts.
    :type: ChemComposition

    :return: The charge adducts as a Mod instance.
    :rtype: Mod

    .. code-block:: python

        >>> write_charge_adducts({'Na+': 2, 'H+': 1})
        Mod('+2Na+,+H+', 1)

        >>> write_charge_adducts({'Na+': 2, 'H+': -1})
        Mod('+2Na+,-H+', 1)

        >>> write_charge_adducts({'I-': 2})
        Mod('+2I-', 1)

        >>> write_charge_adducts({'e-': 1})
        Mod('+e-', 1)

        >>> write_charge_adducts({'Mg2+': 2, 'H+': 1})
        Mod('+2Mg2+,+H+', 1)

    """

    adducts_list = []
    for ion, count in charge_adducts.items():
        sign = '+' if count >= 0 else ''
        sign = '-' if count < 0 else sign
        count_abs = abs(count) if abs(count) != 1 else ''
        adducts_list.append(f"{sign}{count_abs}{ion}")

    return Mod(','.join(adducts_list), 1)


def _pop_ion_count(ion: str) -> Tuple[int, str]:
    """
    Parse the charge of an ion from a string.

    :param ion: The ion string to parse.
    :type ion: str

    :return: A tuple containing the charge and the remaining ion string.
    :rtype: Tuple[int, str]

    .. code-block:: python

        >>> _pop_ion_count('+H+')
        (1, 'H+')

        >>> _pop_ion_count('+2Na+')
        (2, 'Na+')

        >>> _pop_ion_count('2I-')
        (2, 'I-')

        >>> _pop_ion_count('+e-')
        (1, 'e-')

        >>> _pop_ion_count('-2Na+')
        (-2, 'Na+')

        >>> _pop_ion_count('+2Mg2+')
        (2, 'Mg2+')
        
    """

    count_str = ''

    charge = 1
    for i, c in enumerate(ion):
        if c == '-':
            charge = -1
        elif c == '+':
            pass
        elif c.isdigit():
            count_str += c
        else:
            cnt = int(count_str) if count_str else 1
            return cnt * charge, ion[i:]


def _pop_ion_symbol(ion: str) -> Tuple[str, str]:
    """
    Parse the symbol of an ion from a string.

    :param ion: The ion string to parse.
    :type ion: str

    :return: A tuple containing the symbol and the remaining ion string.
    :rtype: Tuple[str, str]

    .. code-block:: python

        >>> _pop_ion_symbol('H+')
        ('H', '+')

        >>> _pop_ion_symbol('Na+')
        ('Na', '+')

        >>> _pop_ion_symbol('I-')
        ('I', '-')

        >>> _pop_ion_symbol('e-')
        ('e', '-')

        >>> _pop_ion_symbol('Na+')
        ('Na', '+')

        >>> _pop_ion_symbol('Mg2+')
        ('Mg', '2+')

    """

    symbol = ''
    for i, c in enumerate(ion):
        if c.isdigit() or c in ['+', '-']:
            return symbol, ion[i:]
        symbol += c
    return symbol, ''


def _pop_ion_charge(ion: str) -> Tuple[int, str]:
    """
    Parse the charge of an ion from a string.

    :param ion: The ion string to parse.
    :type ion: str

    :return: A tuple containing the charge and the remaining ion string.
    :rtype: Tuple[int, str]

    .. code-block:: python

        >>> _pop_ion_charge('+')
        (1, '')

        >>> _pop_ion_charge('+')
        (1, '')

        >>> _pop_ion_charge('-')
        (-1, '')

        >>> _pop_ion_charge('-')
        (-1, '')

        >>> _pop_ion_charge('+')
        (1, '')

        >>> _pop_ion_charge('2+')
        (2, '')

        >>> _pop_ion_charge('2-')
        (-2, '')

    """
    count_str, sign = '', 1
    for c in ion:
        if c == '-':
            sign = -1
        elif c == '+':
            continue
        else:
            count_str += c
    cnt = int(count_str) if count_str else 1
    return cnt * sign, ''


def parse_ion_elements(ion: str) -> Tuple[int, str, int]:
    """
    Parse the count, element, and charge of an ion from a string.

    :param ion: The ion string to parse.
    :type ion: str

    :return: A tuple containing the count, element, and charge of the ion.
    :rtype: Tuple[int, str, int]

    Examples:
        >>> parse_ion_elements('+H+')
        (1, 'H', 1)

        >>> parse_ion_elements('+2Na+')
        (2, 'Na', 1)

        >>> parse_ion_elements('2I-')
        (2, 'I', -1)

        >>> parse_ion_elements('+e-')
        (1, 'e', -1)

        >>> parse_ion_elements('-2Na+')
        (-2, 'Na', 1)

        >>> parse_ion_elements('+2Mg2+')
        (2, 'Mg', 2)

        >>> parse_ion_elements('+2Mg2-')
        (2, 'Mg', -2)

    """
    count, ion = _pop_ion_count(ion)
    symbol, charge = _pop_ion_symbol(ion)
    charge, _ = _pop_ion_charge(charge)
    return count, symbol, charge


def parse_static_mods(mods: Optional[List[ModValue]]) -> Dict[str, List[Mod]]:
    """
    Parse static modifications into a dictionary, mapping the location to the modifications.

    :param mods: List of static modifications, where each modification can be a Mod object or a string representation.
    :type mods: List[ModValue]

    :raises TypeError: If the mod is not a string or Mod instance.

    :return: A dictionary with locations as keys and lists of Mod objects as values.
    :rtype: Dict[str, List[Mod]]

    .. code-block:: python

        >>> parse_static_mods([Mod('[1]@C', 1), Mod('[3.1415]@S,D',1)])
        {'C': [Mod(1, 1)], 'S': [Mod(3.1415, 1)], 'D': [Mod(3.1415, 1)]}

        >>> parse_static_mods([Mod('[1]^2@C', 1), Mod('[3.1415]@S,D', 1)])
        {'C': [Mod(1, 2)], 'S': [Mod(3.1415, 1)], 'D': [Mod(3.1415, 1)]}

        >>> parse_static_mods([Mod('[Formula:[13C]H20]@C', 1)])
        {'C': [Mod('Formula:[13C]H20', 1)]}

        >>> parse_static_mods([Mod('[100]@P', 1)])
        {'P': [Mod(100, 1)]}

        >>> parse_static_mods([])
        {}

    """

    static_mod_dict = {}

    if mods is None:
        return static_mod_dict

    for mod in mods:

        if isinstance(mod, Mod):
            mod = mod.val

        if not isinstance(mod, str):
            raise TypeError(f"Invalid type for static mods: {type(mod)}! Mod or Mod.val must be of type string.")

        mod_info, residues = mod.split('@')
        residues = residues.split(',')
        static_mods = _parse_modifications(mod_info, '[', ']')

        for residue in residues:  # Split on ',' for multiple residues
            static_mod_dict.setdefault(residue, []).extend(static_mods)

    return static_mod_dict


def write_static_mods(mods: Dict[str, List[Mod]]) -> List[Mod]:
    """
    Converts the dictionary of static modifications into a list of Mod instances.

    :param mods: Dictionary of static modifications.
    :type: Dict[str, List[Mod]

    :return: List of static modifications.
    :rtype: List[Mod]

    .. code-block:: python

        >>> write_static_mods({'C': [Mod(val=1,mult=1)], 'S': [Mod(val=3.1415,mult=1)], 'D': [Mod(val=3.1415,mult=1)]})
        [Mod('[1]@C', 1), Mod('[3.1415]@S,D', 1)]

        >>> write_static_mods({'C': [Mod(val='Formula:[13C]H20', mult=1)]})
        [Mod('[Formula:[13C]H20]@C', 1)]

        >>> write_static_mods({})
        []

    """

    reverse_mods = {}

    for residue, mod_list in mods.items():
        mod_tuple = tuple(mod_list)
        reverse_mods.setdefault(mod_tuple, []).append(residue)

    static_mods = []
    for mod_tuple, residues in reverse_mods.items():
        residue_str = ",".join(residues)
        mod_str = ''.join(mod.serialize('[]') for mod in mod_tuple)
        static_mods.append(Mod(f'{mod_str}@{residue_str}', 1))

    return static_mods


def parse_isotope_mods(mods: List[ModValue]) -> Dict[str, str]:
    """
    Parse isotope modifications into a dictionary, mapping the elemental symbol to its isotope.

    :param mods: List of isotope modifications.
    :type mods: List[ModValue]

    :raises TypeError: If the mod is not a string or Mod instance.

    :return: Dictionary of isotope modifications.
    :rtype: Dict[str, str]

    .. code-block:: python

        >>> parse_isotope_mods([Mod('13C', 1), Mod('15N', 1), Mod('D', 1)])
        {'C': '13C', 'N': '15N', 'H': 'D'}

        >>> parse_isotope_mods([Mod('13C', 1), Mod('15N', 1), Mod('T', 1)])
        {'C': '13C', 'N': '15N', 'H': 'T'}

        >>> parse_isotope_mods([Mod('13C', 1), Mod('15N', 1), Mod('H', 1)])
        {'C': '13C', 'N': '15N', 'H': 'H'}

        >>> parse_isotope_mods([])
        {}

    """

    isotope_map = {}
    for mod in mods:

        if isinstance(mod, Mod):
            mod = mod.val

        if not isinstance(mod, str):
            raise TypeError(f"Invalid type for isotope mods: {type(mod)}! Mod or Mod.val must be of type string.")

        # remove digits
        base_aa = re.sub(ISOTOPE_NUM_PATTERN, '', mod)
        isotope_map[base_aa] = mod

    # If any keys are D or T, then replace them with H
    if 'D' in isotope_map:
        isotope_map['H'] = isotope_map.pop('D')

    if 'T' in isotope_map:
        isotope_map['H'] = isotope_map.pop('T')

    return isotope_map


def write_isotope_mods(mods: Dict[str, str]) -> List[Mod]:
    """
    Converts the dictionary of isotope modifications into a list of Mod instances.

    :param mods: Dictionary of isotope modifications.
    :type: Dict[str, str]

    :return: List of isotope modifications.
    :rtype: List[Mod]

    .. code-block:: python

        >>> write_isotope_mods({'C': '13C', 'N': '15N', 'H': 'D'})
        [Mod('13C', 1), Mod('15N', 1), Mod('D', 1)]

        >>> write_isotope_mods({'C': '13C', 'N': '15N', 'H': 'T'})
        [Mod('13C', 1), Mod('15N', 1), Mod('T', 1)]

        >>> write_isotope_mods({})
        []

    """

    return [Mod(val=v, mult=1) for v in mods.values()]


def _validate_single_mod_multiplier(func: Callable) -> Callable:
    """
    A decorator that validates the multiplier of a Mod instance (or each Mod in a list of Mods)
    before adding it through the decorated method. It raises a ValueError if any Mod has a multiplier greater than 1.

    :param func: The function to decorate.
    :type func: function

    :return: The decorated function.
    :rtype: function
    """

    @wraps(func)
    def wrapper(self, mod: Union[Mod, List[Mod]]):
        if isinstance(mod, Mod):
            if mod.mult > 1:
                raise ValueError(
                    f"Invalid multiplier {mod.mult} for mod {mod.val}. Multiplier must not be greater than 1.", None,
                    None)
        else:
            for m in mod:
                if m.mult > 1:
                    raise ValueError(
                        f"Invalid multiplier {m.mult} for mod {m.val}. Multiplier must not be greater than 1.",
                        None, None)

        return func(self, mod)

    return wrapper


@dataclass
class ProFormaAnnotation:
    """
    Represents a ProForma sequence annotation.
    """
    _sequence: str
    _isotope_mods: Optional[List[Mod]] = None
    _static_mods: Optional[List[Mod]] = None
    _labile_mods: Optional[List[Mod]] = None
    _unknown_mods: Optional[List[Mod]] = None
    _nterm_mods: Optional[List[Mod]] = None
    _cterm_mods: Optional[List[Mod]] = None
    _internal_mods: Optional[Dict[int, List[Mod]]] = None
    _intervals: Optional[List[Interval]] = None
    _charge: Optional[int] = None
    _charge_adducts: Optional[Union[List[Mod]]] = None

    def has_sequence(self) -> bool:
        """
        Returns True if a sequence is present, otherwise False.
        """
        return self._sequence is not None

    def has_isotope_mods(self) -> bool:
        """
        Returns True if any isotope modifications are present, otherwise False.
        """
        return self._isotope_mods is not None

    def has_static_mods(self) -> bool:
        """
        Returns True if any static modifications are present, otherwise False.
        """
        return self._static_mods is not None

    def has_labile_mods(self) -> bool:
        """
        Returns True if any labile modifications are present, otherwise False.
        """
        return self._labile_mods is not None

    def has_unknown_mods(self) -> bool:
        """
        Returns True if any unknown modifications are present, otherwise False.
        """
        return self._unknown_mods is not None

    def has_nterm_mods(self) -> bool:
        """
        Returns True if any N-terminal modifications are present, otherwise False.
        """
        return self._nterm_mods is not None

    def has_cterm_mods(self) -> bool:
        """
        Returns True if any C-terminal modifications are present, otherwise False.
        """
        return self._cterm_mods is not None

    def has_internal_mods(self) -> bool:
        """
        Returns True if any internal modifications are present, otherwise False.
        """
        return self._internal_mods is not None

    def has_intervals(self) -> bool:
        """
        Returns True if any intervals are present, otherwise False.
        """
        return self._intervals is not None

    def has_charge(self) -> bool:
        """
        Returns True if a charge is present, otherwise False.
        """
        return self._charge is not None

    def has_charge_adducts(self) -> bool:
        """
        Returns True if any charge adducts are present, otherwise False.
        """
        return self._charge_adducts is not None

    def has_mods(self) -> bool:
        """
        Returns True if any modifications are present, otherwise False.
        """
        return any([self.has_isotope_mods(), self.has_static_mods(), self.has_labile_mods(), self.has_unknown_mods(),
                    self.has_nterm_mods(), self.has_cterm_mods(), self.has_internal_mods(), self.has_intervals(),
                    self.has_charge(), self.has_charge_adducts()])

    @property
    def sequence(self) -> str:
        """
        Returns the sequence.
        """
        return self._sequence

    @sequence.setter
    def sequence(self, value: Optional[str]) -> None:
        """
        Sets the sequence.
        """
        self._sequence = value

    @property
    def isotope_mods(self) -> Union[List[Mod], None]:
        """
        Returns the isotope modifications.
        """
        return self._isotope_mods

    @isotope_mods.setter
    def isotope_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the isotope modifications.
        """
        if value is None:
            self._isotope_mods = None
        else:
            value = fix_list_of_mods(value)
            self._isotope_mods = copy.deepcopy(value)

    @property
    def static_mods(self) -> Union[List[Mod], None]:
        """
        Returns the static modifications.
        """
        return self._static_mods

    @static_mods.setter
    def static_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the static modifications.
        """
        if value is None:
            self._static_mods = None
        else:
            value = fix_list_of_mods(value)
            self._static_mods = copy.deepcopy(value)

    @property
    def labile_mods(self) -> Union[List[Mod], None]:
        """
        Returns the labile modifications.
        """
        return self._labile_mods

    @labile_mods.setter
    def labile_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the labile modifications.
        """
        if value is None:
            self._labile_mods = None
        else:
            value = fix_list_of_mods(value)
            self._labile_mods = copy.deepcopy(value)

    @property
    def unknown_mods(self) -> Union[List[Mod], None]:
        """
        Returns the unknown modifications.
        """
        return self._unknown_mods

    @unknown_mods.setter
    def unknown_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the unknown modifications.
        """
        if value is None:
            self._unknown_mods = None
        else:
            value = fix_list_of_mods(value)
            self._unknown_mods = copy.deepcopy(value)

    @property
    def nterm_mods(self) -> Union[List[Mod], None]:
        """
        Returns the N-terminal modifications.
        """
        return self._nterm_mods

    @nterm_mods.setter
    def nterm_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the N-terminal modifications.
        """
        if value is None:
            self._nterm_mods = None
        else:
            value = fix_list_of_mods(value)
            self._nterm_mods = copy.deepcopy(value)

    @property
    def cterm_mods(self) -> Union[List[Mod], None]:
        """
        Returns the C-terminal modifications.
        """
        return self._cterm_mods

    @cterm_mods.setter
    def cterm_mods(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        """
        Sets the C-terminal modifications.
        """
        if value is None:
            self._cterm_mods = None
        else:
            value = fix_list_of_mods(value)
            self._cterm_mods = copy.deepcopy(value)

    @property
    def internal_mods(self) -> Union[Dict[int, List[Mod]], None]:
        """
        Returns the internal modifications.
        """
        return self._internal_mods

    @internal_mods.setter
    def internal_mods(self, value: Optional[Dict[int, Union[List[ModValue], ModValue]]]) -> None:
        if value is None:
            self._internal_mods = None
        else:
            value = fix_dict_of_mods(value)
            self._internal_mods = copy.deepcopy(value)

    @property
    def intervals(self) -> Union[List[Interval], None]:
        """
        Returns the intervals.
        """
        return self._intervals

    @intervals.setter
    def intervals(self, value: Optional[Union[List[INTERVAL_VALUE], INTERVAL_VALUE]]) -> None:
        """
        Sets the intervals.
        """
        if value is None:
            self._intervals = None
        else:
            value = fix_intervals_input(value)
            self._intervals = copy.deepcopy(value)

    @property
    def charge(self) -> Union[int, None]:
        """
        Returns the charge.
        """
        return self._charge

    @charge.setter
    def charge(self, value: Optional[int]):
        """
        Sets the charge.
        """
        self._charge = value

    @property
    def charge_adducts(self) -> Union[List[Mod], None]:
        """
        Returns the charge adducts.
        """
        return self._charge_adducts

    @charge_adducts.setter
    def charge_adducts(self, value: Optional[Union[List[ModValue], ModValue]]) -> None:
        if value is None:
            self._charge_adducts = None
        else:
            value = fix_list_of_mods(value)
            self._charge_adducts = copy.deepcopy(value)

    def __repr__(self) -> str:
        """
        Only shows non None types
        """
        seq = 'ProFormaAnnotation('
        for k, v in self.dict().items():
            if v is not None:
                seq += f"{k}={v}, "
        seq = seq.rstrip(', ')
        seq += ')'

        return seq

    def __len__(self) -> int:
        return len(self.sequence)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, ProFormaAnnotation):
            return NotImplemented

        # Direct sequence comparison
        if self._sequence != other._sequence:
            return False

        if not are_mods_equal(self.labile_mods, other.labile_mods):
            return False

        if not are_mods_equal(self.unknown_mods, other.unknown_mods):
            return False

        if not are_mods_equal(self.nterm_mods, other.nterm_mods):
            return False

        if not are_mods_equal(self.cterm_mods, other.cterm_mods):
            return False

        if not are_mods_equal(self.charge_adducts, other.charge_adducts):
            return False

        if not are_mods_equal(self.isotope_mods, other.isotope_mods):
            return False

        if not are_mods_equal(self.static_mods, other.static_mods):
            return False

        # Compare internal modifications (more complex due to being a dictionary)
        if self.has_internal_mods() or other.has_internal_mods():
            self_keys = set(self.internal_mods.keys()) if self.has_internal_mods() else set()
            other_keys = set(other.internal_mods.keys()) if other.has_internal_mods() else set()

            for k in self_keys.union(other_keys):
                if not are_mods_equal(self.get_internal_mods_by_index(k), other.get_internal_mods_by_index(k)):
                    return False

        # Compare intervals (if applicable and need to ensure a deep comparison)
        if not are_intervals_equal(self.intervals, other.intervals):
            return False

        # Compare charges
        if self.charge != other.charge:
            return False

        return True

    def copy(self) -> 'ProFormaAnnotation':
        """
        Returns a deep copy of the ProFormaAnnotation instance.
        """
        return deepcopy(self)

    def clear_empty_mods(self) -> None:
        """
        Sets the mods to None if they are an empty list
        """
        if self.has_isotope_mods() and len(self.isotope_mods) == 0:
            self._isotope_mods = None
        if self.has_static_mods() and len(self.static_mods) == 0:
            self._static_mods = None
        if self.has_labile_mods() and len(self.labile_mods) == 0:
            self._labile_mods = None
        if self.has_unknown_mods() and len(self.unknown_mods) == 0:
            self._unknown_mods = None
        if self.has_nterm_mods() and len(self.nterm_mods) == 0:
            self._nterm_mods = None
        if self.has_cterm_mods() and len(self.cterm_mods) == 0:
            self._cterm_mods = None
        if self.has_charge_adducts() and len(self.charge_adducts) == 0:
            self._charge_adducts = None

        # Internal mods
        if self.has_internal_mods():
            keys_to_remove = []
            for k, mods in self.internal_mods.items():
                if mods is not None and len(mods) == 0:
                    keys_to_remove.append(k)

            for k in keys_to_remove:
                self.internal_mods.pop(k)

            if len(self.internal_mods) == 0:
                self.internal_mods = None

        # Intervals
        if self.has_intervals():
            for interval in self.intervals:
                if interval.mods is not None and len(interval.mods) == 0:
                    interval.mods = None

            if len(self.intervals) == 0:
                self.intervals = None

    def dict(self) -> Dict[str, Any]:
        """
        Returns a dictionary representation of the ProFormaAnnotation instance.
        """
        return {
            "sequence": copy.deepcopy(self.sequence),
            "isotope_mods": copy.deepcopy(self.isotope_mods),
            "static_mods": copy.deepcopy(self.static_mods),
            "labile_mods": copy.deepcopy(self.labile_mods),
            "unknown_mods": copy.deepcopy(self.unknown_mods),
            "nterm_mods": copy.deepcopy(self.nterm_mods),
            "cterm_mods": copy.deepcopy(self.cterm_mods),
            "intervals": copy.deepcopy(self.intervals),
            "internal_mods": copy.deepcopy(self.internal_mods),
            "charge": copy.deepcopy(self.charge),
            "charge_adducts": copy.deepcopy(self.charge_adducts)
        }

    def mod_dict(self) -> Dict[str, Any]:
        """
        Returns a dictionary representation of the ProFormaAnnotation instance, including only the
        modification-related attributes.
        """
        result = {
            # Only include lists and dicts that are not None, maintaining a clean output
            **{k: v for k, v in {
                "isotope": self.isotope_mods,
                "static": self.static_mods,
                "labile": self.labile_mods,
                "unknown": self.unknown_mods,
                "nterm": self.nterm_mods,
                "cterm": self.cterm_mods,
                "intervals": self.intervals,
                "charge": self.charge,
                "charge_adducts": self.charge_adducts,
            }.items() if v is not None}
        }

        # add internal mods as numbered keys
        if self.internal_mods:
            for index, mods in self.internal_mods.items():
                result[index] = mods

        return copy.deepcopy(result)

    def add_mod_dict(self, mod_dict: Dict, append: bool = False) -> None:
        """
        Add mods from a dictionary
        """
        if 'isotope' in mod_dict:
            self.add_isotope_mods(mod_dict['isotope'], append)
        if 'static' in mod_dict:
            self.add_static_mods(mod_dict['static'], append)
        if 'labile' in mod_dict:
            self.add_labile_mods(mod_dict['labile'], append)
        if 'unknown' in mod_dict:
            self.add_unknown_mods(mod_dict['unknown'], append)
        if 'nterm' in mod_dict:
            self.add_nterm_mods(mod_dict['nterm'], append)
        if 'cterm' in mod_dict:
            self.add_cterm_mods(mod_dict['cterm'], append)
        if 'intervals' in mod_dict:
            self.add_intervals(mod_dict['intervals'], append)
        if 'charge' in mod_dict:
            self.charge = mod_dict['charge']
        if 'charge_adducts' in mod_dict:
            self.add_charge_adducts(mod_dict['charge_adducts'], append)

        internal_mods = {}
        for k, v in mod_dict.items():
            if isinstance(k, int):
                internal_mods[k] = v

        if len(internal_mods) > 0:
            self.add_internal_mods(internal_mods, append)

    def condense_static_mods(self, inplace: bool = False) -> Union['ProFormaAnnotation', None]:
        """
        Condense static mods into internal mods
        """
        if inplace is False:
            new_annotation = deepcopy(self)
        else:
            new_annotation = self

        static_mods = new_annotation.pop_static_mods()
        if static_mods is None:
            if inplace is False:
                return new_annotation

        static_mod_dict = parse_static_mods(static_mods)

        n_term_mod = static_mod_dict.get('N-Term')
        if n_term_mod is not None:
            new_annotation.add_nterm_mods(n_term_mod, append=True)

        c_term_mod = static_mod_dict.get('C-Term')
        if c_term_mod is not None:
            new_annotation.add_cterm_mods(c_term_mod, append=True)

        for aa, mod in static_mod_dict.items():

            if aa in ['N-Term', 'C-Term']:
                continue

            # get indexes of the amino acid
            indexes = [m.start() for m in re.finditer(aa, new_annotation.sequence)]
            for index in indexes:
                new_annotation.add_internal_mods({index: mod}, append=True)

        if inplace is False:
            return new_annotation

    def contains_sequence_ambiguity(self) -> bool:
        """
        Check if the sequence contains any ambiguity (modifications or intervals).
        """
        return self.intervals is not None or self.unknown_mods is not None or self.labile_mods is not None

    def contains_residue_ambiguity(self) -> bool:
        """
        Check if the sequence contains any residue ambiguous amino acids.
        """
        return any(aa in AMBIGUOUS_AMINO_ACIDS for aa in self.sequence)

    def contains_mass_ambiguity(self) -> bool:
        """
        Check if the sequence contains any mass ambiguous amino acids.
        """
        return any(aa in MASS_AMBIGUOUS_AMINO_ACIDS for aa in self.sequence)

    def pop_labile_mods(self) -> Union[List[Mod], None]:
        """
        Pop all labile mods and return them in a list
        """
        m = self.labile_mods
        self.labile_mods = None
        return m

    def pop_unknown_mods(self) -> Union[List[Mod], None]:
        """
        Pop all unknown mods and return them in a list
        """
        m = self.unknown_mods
        self.unknown_mods = None
        return m

    def pop_nterm_mods(self) -> Union[List[Mod], None]:
        """
        Pop all nterm mods and return them in a list
        """
        m = self.nterm_mods
        self.nterm_mods = None
        return m

    def pop_cterm_mods(self) -> Union[List[Mod], None]:
        """
        Pop all cterm mods and return them in a list
        """
        m = self.cterm_mods
        self.cterm_mods = None
        return m

    def pop_internal_mods(self) -> Union[Dict[int, List[Mod]], None]:
        """
        Pop all internal mods and return them in a dictionary
        """
        m = self.internal_mods
        self.internal_mods = None
        return m

    def pop_intervals(self) -> Union[List[Interval], None]:
        """
        Pop all intervals and return them in a list
        """
        m = self.intervals
        self.intervals = None
        return m

    def pop_charge(self) -> Union[int, None]:
        """
        Pop the charge and return it
        """
        m = self.charge
        self.charge = None
        return m

    def pop_charge_adducts(self) -> Union[str, None]:
        """
        Pop all charge adducts and return them in a list
        """
        m = self.charge_adducts
        self.charge_adducts = None
        return m

    def pop_isotope_mods(self) -> Union[List[Mod], None]:
        """
        Pop all isotope mods and return them in a list
        """
        m = self.isotope_mods
        self.isotope_mods = None
        return m

    def pop_static_mods(self) -> Union[List[Mod], None]:
        """
        Pop all static mods and return them in a list
        """
        m = self.static_mods
        self.static_mods = None
        return m

    def pop_mods(self) -> Dict[str, Any]:
        """
        Pop all mods and return them in a dictionary
        """
        d = {}
        if self.has_isotope_mods():
            d['isotope'] = self.pop_isotope_mods()
        if self.has_static_mods():
            d['static'] = self.pop_static_mods()
        if self.has_labile_mods():
            d['labile'] = self.pop_labile_mods()
        if self.has_unknown_mods():
            d['unknown'] = self.pop_unknown_mods()
        if self.has_nterm_mods():
            d['nterm'] = self.pop_nterm_mods()
        if self.has_cterm_mods():
            d['cterm'] = self.pop_cterm_mods()
        if self.has_charge_adducts():
            d['charge_adducts'] = self.pop_charge_adducts()
        if self.has_charge():
            d['charge'] = self.pop_charge()
        if self.has_internal_mods():
            d['internal'] = self.pop_internal_mods()
        if self.has_intervals():
            d['intervals'] = self.pop_intervals()

        return d

    def add_labile_mods(self, mods: Optional[Union[List[Mod], Mod]], append: bool = False) -> None:
        """
        Add labile mods to the annotation. If not append, existing mods will be replaced.
        """
        if mods is None:
            if not append:
                self.labile_mods = None
            return

        if isinstance(mods, Mod):
            mods = [mods]

        if not append:
            self.labile_mods = mods
        else:
            if self.has_labile_mods():
                self.labile_mods.extend(copy.deepcopy(mods))
            else:
                self.labile_mods = mods

    def add_unknown_mods(self, mods: Optional[Union[List[ModValue], ModValue]], append: bool = False) -> None:
        """
        Add unknown mods to the annotation. If not append, existing mods will be replaced.
        """
        if mods is None:
            if not append:
                self.unknown_mods = None
            return

        mods = fix_list_of_mods(mods)

        if not append:
            self.unknown_mods = mods
        else:
            if self.has_unknown_mods():
                self.unknown_mods.extend(copy.deepcopy(mods))
            else:
                self.unknown_mods = mods

    def add_nterm_mods(self, mods: Optional[Union[List[ModValue], ModValue]], append: bool = False) -> None:
        """
        Add nterm mods to the annotation. If not append, existing mods will be replaced.
        """
        if mods is None:
            if not append:
                self.nterm_mods = None
            return

        mods = fix_list_of_mods(mods)

        if not append:
            self.nterm_mods = mods
        else:
            if self.has_nterm_mods():
                self.nterm_mods.extend(copy.deepcopy(mods))
            else:
                self.nterm_mods = mods

    def add_cterm_mods(self, mods: Optional[Union[List[ModValue], ModValue]], append: bool = False) -> None:
        """
        Add cterm mods to the annotation. If not append, existing mods will be replaced.
        """
        if mods is None:
            if not append:
                self.cterm_mods = None
            return

        mods = fix_list_of_mods(mods)

        if not append:
            self.cterm_mods = mods  # Uses the setter to ensure proper copying
        else:
            if self.has_cterm_mods():
                self.cterm_mods.extend(copy.deepcopy(mods))
            else:
                self.cterm_mods = mods

    def count_internal_mods(self) -> int:
        """
        Count the number of internal mods in the annotation
        """
        if not self.has_internal_mods():
            return 0
        return sum(len(v) for v in self.internal_mods.values())

    def count_modified_residues(self) -> int:
        """
        Count the number of modified residues in the annotation
        """
        if not self.has_internal_mods():
            return 0
        return len(self.internal_mods)

    def has_internal_mods_at_index(self, index: int) -> bool:
        """
        Check if there are internal mods at a specific index
        """
        if not self.has_internal_mods():
            return False
        return index in self.internal_mods

    def get_internal_mods_by_index(self, index: int) -> Union[List[Mod], None]:
        """
        Get internal mods at a specific index
        """
        if not self.has_internal_mods():
            return None

        if index not in self.internal_mods:
            return None

        return self.internal_mods[index]

    def pop_internal_mod(self, index: int) -> Union[List[Mod], None]:
        """
        Pop internal mods at a specific index
        """
        if not self.has_internal_mods():
            return None

        if index not in self.internal_mods:
            return None

        return self.internal_mods.pop(index)

    def add_internal_mod(self, index: int, mods: Optional[Union[List[ModValue], ModValue]],
                         append: bool = False) -> None:
        """
        Add internal mods to the annotation. If not append, existing mods will be replaced.
        """

        if mods is None:
            if not append:
                if not self.has_internal_mods():
                    self._internal_mods = {}
                self._internal_mods.pop(index, None)
            return

        mods = fix_list_of_mods(mods)

        if not self.has_internal_mods():
            self._internal_mods = {}

        if not append:
            self._internal_mods[index] = copy.deepcopy(mods)
        else:
            if index in self.internal_mods:
                self._internal_mods[index].extend(copy.deepcopy(mods))
            else:
                self._internal_mods[index] = copy.deepcopy(mods)

    def add_internal_mods(self, mods: Optional[Dict[int, Union[List[ModValue], ModValue]]],
                          append: bool = False) -> None:
        """
        Add internal mods to the annotation. If not append, existing mods will be replaced.
        """

        if mods is None:
            if not append:
                self.internal_mods = None
            return

        mods = fix_dict_of_mods(mods)

        if not append:
            self.internal_mods = mods
        else:
            if not self.has_internal_mods():
                self.internal_mods = copy.deepcopy(mods)
            else:
                for k, v in mods.items():
                    if k in self.internal_mods:
                        self.internal_mods[k].extend(copy.deepcopy(v))
                    else:
                        self.internal_mods[k] = copy.deepcopy(v)

    def add_intervals(self, intervals: Optional[Union[List[INTERVAL_VALUE], INTERVAL_VALUE]],
                      append: bool = False) -> None:
        """
        Add intervals to the annotation. If not append, existing mods will be replaced.
        """

        if intervals is None:
            if append is False:
                self.intervals = None
            return

        intervals = fix_intervals_input(intervals)

        if append is False:
            self.intervals = intervals  # Uses the setter to ensure proper copying
        else:
            if self.has_intervals():
                self.intervals.extend(copy.deepcopy(intervals))
            else:
                self.intervals = intervals

    def add_charge(self, charge: Optional[int]) -> None:
        """
        Add charge to the annotation
        """
        self.charge = charge

    def add_charge_adducts(self, charge_adducts: Optional[Union[List[ModValue], ModValue]],
                           append: bool = False) -> None:
        """
        Add charge adducts to the annotation
        """
        if charge_adducts is None:
            if not append:
                self.charge_adducts = None
            return

        charge_adducts = fix_list_of_mods(charge_adducts)

        if not append:
            self.charge_adducts = charge_adducts  # Uses the setter to ensure proper copying
        else:
            if self.has_charge_adducts():
                self.charge_adducts.extend(copy.deepcopy(charge_adducts))
            else:
                self.charge_adducts = charge_adducts

    def add_isotope_mods(self, mods: Optional[Union[List[ModValue], ModValue]], append: bool = False) -> None:
        """
        Add isotope mods to the annotation. If not append, existing mods will be replaced.
        """

        if mods is None:
            if not append:
                self.isotope_mods = None
            return

        mods = fix_list_of_mods(mods)

        if not append:
            self.isotope_mods = mods  # Uses the setter to ensure proper copying
        else:
            if self.has_isotope_mods():
                self.isotope_mods.extend(copy.deepcopy(mods))
            else:
                self.isotope_mods = mods

    def add_static_mods(self, mods: Optional[Union[List[ModValue], ModValue]], append: bool = False) -> None:
        """
        Add static mods to the annotation. If not append, existing mods will be replaced.
        """

        if mods is None:
            if not append:
                self.static_mods = None
            return

        mods = fix_list_of_mods(mods)

        if not append:
            self.static_mods = mods  # Uses the setter to ensure proper copying
        else:
            if self.has_static_mods():
                self.static_mods.extend(copy.deepcopy(mods))
            else:
                self.static_mods = mods

    def strip(self, inplace: bool = False) -> Union['ProFormaAnnotation', None]:
        """
        Remove all modifications from the annotation and return a new annotation with the stripped sequence.
        """

        if inplace:
            self.isotope_mods = None
            self.static_mods = None
            self.labile_mods = None
            self.unknown_mods = None
            self.nterm_mods = None
            self.cterm_mods = None
            self.internal_mods = None
            self.intervals = None
            self.charge = None
            self.charge_adducts = None
            return None

        return ProFormaAnnotation(_sequence=self.sequence)

    def slice(self, start: Optional[int], stop: Optional[int], inplace: bool = False) -> (
            Union)['ProFormaAnnotation', None]:
        """
        Slice the annotation sequence and return a new annotation with the sliced sequence and modifications.
        """

        if start is None:
            start = 0

        if stop is None:
            stop = len(self.sequence)

        new_sequence = self.sequence[start:stop]

        if not self.has_mods():
            if inplace is True:
                self._sequence = new_sequence
                return None
            return ProFormaAnnotation(_sequence=new_sequence)

        # Adjust internal modifications based on new sequence indices
        new_internal_mods = None
        if self.has_internal_mods():
            new_internal_mods = {}
            for k, mods in self.internal_mods.items():
                if start <= k < stop:
                    new_internal_mods[k - start] = copy.deepcopy(mods)

        new_intervals = None
        if self.has_intervals():
            new_intervals = []
            for interval in self.intervals:
                if interval.start < stop and interval.end >= start:
                    new_start = max(0, interval.start - start)
                    new_end = max(0, interval.end - start)
                    new_intervals.append(Interval(start=new_start,
                                                  end=new_end,
                                                  ambiguous=interval.ambiguous,
                                                  mods=copy.deepcopy(interval.mods)))

        if inplace is True:
            self._sequence = new_sequence
            self._internal_mods = new_internal_mods  # already a copy
            self._intervals = new_intervals  # already a copy
            if start > 0:
                self._nterm_mods = None
            if stop < len(self.sequence):
                self._cterm_mods = None
            return None

        # Create a new annotation with the sliced sequence and modifications
        new_annotation = copy.deepcopy(self)
        new_annotation._sequence = new_sequence
        new_annotation._internal_mods = new_internal_mods  # already a copy
        new_annotation._intervals = new_intervals  # already a copy
        if start > 0:
            new_annotation._nterm_mods = None
        if stop < len(self.sequence):
            new_annotation._cterm_mods = None

        return new_annotation

    def shift(self, n: int, inplace: bool = False) -> Union['ProFormaAnnotation', None]:
        """
        Shift the annotation by n positions in a cyclic manner.
        """
        seq_len = len(self.sequence)
        effective_shift = n % seq_len
        shifted_sequence = self.sequence[effective_shift:] + self.sequence[:effective_shift]

        new_internal_mods = None
        if self.has_internal_mods():
            new_internal_mods = {}
            for mod_index, mods in self.internal_mods.items():
                shifted_index = (mod_index - effective_shift) % seq_len
                new_internal_mods[shifted_index] = copy.deepcopy(mods)

            if len(new_internal_mods) == 0:
                new_internal_mods = None

        # Adjust intervals considering the effective shift and sequence length
        new_intervals = None
        if self.has_intervals():
            new_intervals = []
            for interval in self.intervals:
                new_start = (interval.start - effective_shift) % seq_len
                new_end = (interval.end - effective_shift) % seq_len if interval.end is not None else None
                # Ensure the start is always less than the end for non-ambiguous intervals
                if new_end is not None and new_start > new_end:
                    new_end, new_start = new_start, new_end
                new_intervals.append(
                    Interval(start=new_start,
                             end=new_end,
                             ambiguous=interval.ambiguous,
                             mods=copy.deepcopy(interval.mods)))

            if len(new_intervals) == 0:
                new_intervals = None

        if inplace is True:
            self._sequence = shifted_sequence
            self._internal_mods = new_internal_mods  # already a copy
            self._intervals = new_intervals  # already a copy
            return None

        # Create a new annotation with the shifted sequence and modifications
        new_annotation = copy.deepcopy(self)
        new_annotation._sequence = shifted_sequence
        new_annotation._internal_mods = new_internal_mods  # already a copy
        new_annotation._intervals = new_intervals  # already a copy

        return new_annotation

    def shuffle(self, seed: Optional[Any] = None, inplace: bool = False) -> Union['ProFormaAnnotation', None]:
        """
        Shuffle the annotation sequence and return a new annotation with the shuffled sequence.
        """
        if seed is not None:
            random.seed(seed)

        # Convert sequence to a list of characters for shuffling
        sequence_list = list(self.sequence)
        # Track original positions
        original_positions = list(range(len(sequence_list)))
        # Shuffle the sequence list
        combined = list(zip(sequence_list, original_positions))
        random.shuffle(combined)
        shuffled_sequence, shuffled_positions = zip(*combined)

        # Shuffle internal modifications based on new positions
        new_internal_mods = None
        if self.internal_mods:
            new_internal_mods = {}
            # Create a mapping from original to new positions
            position_mapping = {original: new for new, original in enumerate(shuffled_positions)}
            for original_pos, mods in self.internal_mods.items():
                # Map each original position to its new position
                new_pos = position_mapping[original_pos]
                new_internal_mods[new_pos] = copy.deepcopy(mods)

            if len(new_internal_mods) == 0:
                new_internal_mods = None

        new_sequence = ''.join(shuffled_sequence)
        if inplace:
            self._sequence = new_sequence
            self._internal_mods = new_internal_mods
            return None

        new_annotation = copy.deepcopy(self)
        new_annotation._sequence = new_sequence
        new_annotation._internal_mods = new_internal_mods
        return new_annotation

    def reverse(self, inplace: bool = False, swap_terms: bool = False) -> Union['ProFormaAnnotation', None]:
        """
        Reverse the annotation sequence and return a new annotation with the reversed sequence.
        """
        reversed_sequence = self.sequence[::-1]

        # Reverse internal modifications based on new positions
        new_internal_mods = {}
        if self.internal_mods:
            for original_pos, mods in self.internal_mods.items():
                new_pos = len(self.sequence) - original_pos - 1
                new_internal_mods[new_pos] = copy.deepcopy(mods)

        if len(new_internal_mods) == 0:
            new_internal_mods = None

        # reverse intervals
        new_intervals = None
        if self.intervals is not None:
            new_intervals = []
            for interval in self.intervals:
                new_start = len(self.sequence) - interval.start - 1
                new_end = len(self.sequence) - interval.end - 1 if interval.end is not None else None
                # Ensure the start is always less than the end for non-ambiguous intervals
                if new_end is not None and new_start > new_end:
                    new_end, new_start = new_start, new_end
                new_intervals.append(
                    Interval(start=new_start,
                             end=new_end,
                             ambiguous=interval.ambiguous,
                             mods=copy.deepcopy(interval.mods)))

        if swap_terms:
            nterm_mods = self.cterm_mods
            cterm_mods = self.nterm_mods
        else:
            nterm_mods = self.nterm_mods
            cterm_mods = self.cterm_mods

        if inplace:
            self._sequence = reversed_sequence  # already a copy
            self._internal_mods = new_internal_mods  # already a copy
            self._nterm_mods = nterm_mods
            self._cterm_mods = cterm_mods
            self._intervals = new_intervals
            return None

        # Create a new annotation with the reversed sequence and modifications
        new_annotation = copy.deepcopy(self)
        new_annotation._sequence = reversed_sequence  # already a copy
        new_annotation._internal_mods = new_internal_mods  # already a copy
        new_annotation._nterm_mods = copy.deepcopy(nterm_mods)
        new_annotation._cterm_mods = copy.deepcopy(cterm_mods)
        new_annotation._intervals = new_intervals  # already a copy
        return new_annotation

    def split(self) -> Generator['ProFormaAnnotation', None, None]:
        """
        Split each amino acid in the sequence into a separate ProFormaAnnotation
        """

        # for labile mods only include on first amino acid
        labile_mods = self.pop_labile_mods()

        for i, _ in enumerate(self.sequence):
            s = self.slice(i, i + 1)
            if i == 0 and labile_mods:
                s.add_labile_mods(labile_mods, append=True)
            yield s

    def count_residues(self) -> CounterType:
        """
        Count the occurrences of each residue in the sequence.
        """
        return Counter([a.serialize() for a in self.split()])

    def sort_residues(self, inplace: bool = False) -> Union['ProFormaAnnotation', None]:
        """
        Sort the residues in the annotation sequence and return a new annotation with the sorted sequence.
        """
        # Mapping original positions to their new positions after sorting
        original_to_new_positions = {old: new for new, old in
                                     enumerate(sorted(range(len(self.sequence)), key=lambda x: self.sequence[x]))}

        # Creating new internal mods with adjusted positions
        new_internal_mods = {}
        if self.internal_mods:
            for pos, mods in self.internal_mods.items():
                new_pos = original_to_new_positions[pos]
                new_internal_mods[new_pos] = copy.deepcopy(mods)

        if len(new_internal_mods) == 0:
            new_internal_mods = None

        # Generating sorted sequence
        sorted_sequence = ''.join(sorted(self.sequence))

        if inplace:
            self.sequence = sorted_sequence
            self.internal_mods = new_internal_mods
            return None

        new_annotation = copy.deepcopy(self)
        new_annotation.sequence = sorted_sequence
        new_annotation._internal_mods = new_internal_mods  # already a copy
        return new_annotation

    def serialize(self, include_plus: bool = False) -> str:
        """
        Serialize the entire annotation
        """
        return _serialize_annotation(self, include_plus)

    def serialize_start(self, include_plus: bool = False) -> str:
        """
        Serialize the start of the annotation
        """
        return _serialize_annotation_start(self, include_plus)

    def serialize_middle(self, include_plus: bool = False) -> str:
        """
        Serialize the middle of the annotation
        """
        return _serialize_annotation_middle(self, include_plus)

    def serialize_end(self, include_plus: bool = False) -> str:
        """
        Serialize the end of the annotation
        """
        return _serialize_annotation_end(self, include_plus)

    def is_subsequence(self, other: 'ProFormaAnnotation') -> bool:
        """
        Check if the annotation is a subsequence of another annotation.
        """

        # check if unmodified sequence is a subsequence
        if self.sequence in other.sequence:

            # loop over all starting indexes where the sequence is found
            for start in [m.start() for m in re.finditer(self.sequence, other.sequence)]:
                # check if all modifications are also a subsequence
                sliced_annot = other.slice(start, start + len(self.sequence))
                if sliced_annot == self:
                    return True

        return False

    def find_indices(self, other: 'ProFormaAnnotation') -> List[int]:
        """
        Find all occurrences of the annotation in another annotation.
        """

        # find all starting indexes where the sequence is found and is a subsequence
        return [m.start() for m in re.finditer(self.sequence, other.sequence) if
                self.is_subsequence(other.slice(m.start(), m.start() + len(self.sequence)))]

    def permutations(self, size: Optional[int] = None) -> List['ProFormaAnnotation']:
        """
        Generate all permutations of the annotation sequence.
        """

        if size is None:
            size = len(self)

        start = self.serialize_start()
        end = self.serialize_end()

        mods = self.pop_mods()
        self._internal_mods = mods.get('internal')

        components = [a.serialize() for a in self.split()]

        return [parse(start + ''.join(i) + end) for i in itertools.permutations(components, size)]

    def product(self, repeat: Optional[int] = None) -> List['ProFormaAnnotation']:
        """
        Generate the product of the annotation sequence with itself.
        """

        if repeat is None:
            repeat = len(self)

        start = self.serialize_start()
        end = self.serialize_end()

        mods = self.pop_mods()
        self._internal_mods = mods.get('internal')

        components = [a.serialize() for a in self.split()]

        return [parse(start + ''.join(i) + end) for i in itertools.product(components, repeat=repeat)]

    def combinations(self, size: Optional[int] = None) -> List['ProFormaAnnotation']:
        """
        Generate all combinations of the annotation sequence.
        """

        if size is None:
            size = len(self)

        start = self.serialize_start()
        end = self.serialize_end()

        mods = self.pop_mods()
        self._internal_mods = mods.get('internal')

        components = [a.serialize() for a in self.split()]

        return [parse(start + ''.join(i) + end) for i in itertools.combinations(components, size)]

    def combinations_with_replacement(self, size: Optional[int] = None) -> List['ProFormaAnnotation']:
        """
        Generate all combinations of the annotation sequence with replacement.
        """

        if size is None:
            size = len(self)

        start = self.serialize_start()
        end = self.serialize_end()

        mods = self.pop_mods()
        self._internal_mods = mods.get('internal')

        components = [a.serialize() for a in self.split()]

        return [parse(start + ''.join(i) + end) for i in itertools.combinations_with_replacement(components, size)]


def create_annotation(sequence: str,
                      isotope_mods: Optional[ACCEPTED_MOD_INPUT] = None,
                      static_mods: Optional[ACCEPTED_MOD_INPUT] = None,
                      labile_mods: Optional[ACCEPTED_MOD_INPUT] = None,
                      unknown_mods: Optional[ACCEPTED_MOD_INPUT] = None,
                      nterm_mods: Optional[ACCEPTED_MOD_INPUT] = None,
                      cterm_mods: Optional[ACCEPTED_MOD_INPUT] = None,
                      internal_mods: Optional[Dict[int, ACCEPTED_MOD_INPUT]] = None,
                      intervals: Optional[ACCEPTED_INTERVAL_INPUT] = None,
                      charge: Optional[int] = None,
                      charge_adducts: Optional[ACCEPTED_MOD_INPUT] = None) -> ProFormaAnnotation:
    """
    Create a ProFormaAnnotation from a sequence and modifications

    .. code-block:: python

        >>> create_annotation('PEPTIDE', static_mods=['Carbamidomethyl'])
        ProFormaAnnotation(sequence=PEPTIDE, static_mods=[Mod('Carbamidomethyl', 1)])

    """

    isotope_mods = fix_list_of_mods(isotope_mods) if isotope_mods is not None else None
    static_mods = fix_list_of_mods(static_mods) if static_mods is not None else None
    labile_mods = fix_list_of_mods(labile_mods) if labile_mods is not None else None
    unknown_mods = fix_list_of_mods(unknown_mods) if unknown_mods is not None else None
    nterm_mods = fix_list_of_mods(nterm_mods) if nterm_mods is not None else None
    cterm_mods = fix_list_of_mods(cterm_mods) if cterm_mods is not None else None
    internal_mods = fix_dict_of_mods(internal_mods) if internal_mods is not None else None
    intervals = fix_intervals_input(intervals) if intervals is not None else None
    charge_adducts = fix_list_of_mods(charge_adducts) if charge_adducts is not None else None

    return ProFormaAnnotation(
        _sequence=sequence,
        _isotope_mods=isotope_mods,
        _static_mods=static_mods,
        _labile_mods=labile_mods,
        _unknown_mods=unknown_mods,
        _nterm_mods=nterm_mods,
        _cterm_mods=cterm_mods,
        _internal_mods=internal_mods,
        _intervals=intervals,
        _charge=charge,
        _charge_adducts=charge_adducts
    )


@dataclass
class MultiProFormaAnnotation:
    """
    A multi proforma annotation
    """
    annotations: List[ProFormaAnnotation]
    connections: List[bool]

    def serialize(self, include_plus: bool = False) -> str:
        """
        Convert the multi annotation to a proforma string.

        :return: The serialized multi annotation
        :rtype: str
        """
        seq = ''
        for i, annotation in enumerate(self.annotations):
            seq += annotation.serialize(include_plus=include_plus)
            if i != len(self.annotations) - 1:
                connection = self.connections[i]
                if connection is True:
                    seq += r'\\'
                else:
                    seq += r'+'

        return seq


def create_multi_annotation(annotations: List[ProFormaAnnotation], connections: List[bool]) -> MultiProFormaAnnotation:
    """
    Create a MultiProFormaAnnotation from a list of annotations and connections

    :param annotations: The list of annotations
    :type annotations: List[ProFormaAnnotation]

    :raises ValueError: The number of connections should be one less than the number of annotations

    :param connections: The list of connections
    :type connections: List[bool]

    .. code-block:: python

        >>> annotation = create_multi_annotation([create_annotation('PEP'), create_annotation('TIDE')], [False])
        >>> annotation.annotations[0]
        ProFormaAnnotation(sequence=PEP)
        >>> annotation.annotations[1]
        ProFormaAnnotation(sequence=TIDE)
        >>> annotation.connections
        [False]
    """

    if len(annotations) != len(connections) + 1:
        raise ValueError('The number of connections should be one less than the number of annotations')

    return MultiProFormaAnnotation(annotations=annotations, connections=connections)


def _is_unmodified(proforma_sequence: str) -> bool:
    """
    Check if a proforma sequence is unmodified

    :param proforma_sequence: The proforma sequence
    :type proforma_sequence: str

    :return: True if the sequence is unmodified
    :rtype: bool
    """
    return all(c in AMINO_ACIDS for c in proforma_sequence)


class _ProFormaParser:
    """
    A proforma sequence parser
    """

    def __init__(self, proforma_sequence: str):
        self.sequence = proforma_sequence
        self.position = 0
        self.length = len(proforma_sequence)
        self._amino_acids: List[str] = []
        self._isotope_mods = None
        self._static_mods = None
        self._labile_mods = None
        self._unknown_mods = None
        self._nterm_mods = None
        self._cterm_mods = None
        self._internal_mods = None
        self._charge = None
        self._charge_adducts = None
        self._intervals = None
        self._current_connection = None

    def parse(self) -> Generator[Tuple[ProFormaAnnotation, bool], None, None]:
        """
        Parse the proforma sequence, yielding annotations and connections

        :return: A generator of annotations and connections
        :rtype: Generator[(ProFormaAnnotation, bool), None, None]
        """
        while not self._end_of_sequence():
            self._parse_sequence_start()
            self._parse_sequence_middle()
            self._parse_sequence_end()

            yield self._get_result(), self._current_connection

            if not self._end_of_sequence():
                self._reset_sequence()

    @property
    def _unmod_sequence(self) -> str:
        return ''.join(self._amino_acids)

    def _get_result(self) -> ProFormaAnnotation:
        return ProFormaAnnotation(
            _sequence=self._unmod_sequence,
            _isotope_mods=self._isotope_mods,
            _static_mods=self._static_mods,
            _labile_mods=self._labile_mods,
            _unknown_mods=self._unknown_mods,
            _nterm_mods=self._nterm_mods,
            _cterm_mods=self._cterm_mods,
            _internal_mods=self._internal_mods,
            _intervals=self._intervals,
            _charge=self._charge,
            _charge_adducts=self._charge_adducts
        )

    def _reset_sequence(self) -> None:
        self._amino_acids = []
        self._isotope_mods = None
        self._static_mods = None
        self._labile_mods = None
        self._unknown_mods = None
        self._nterm_mods = None
        self._cterm_mods = None
        self._internal_mods = None
        self._charge = None
        self._charge_adducts = None
        self._intervals = None

    @_validate_single_mod_multiplier
    def _add_static_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._static_mods is None:
            self._static_mods = []
        if isinstance(mod, list):
            self._static_mods.extend(mod)
        else:
            self._static_mods.append(mod)

    @_validate_single_mod_multiplier
    def _add_isotope_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._isotope_mods is None:
            self._isotope_mods = []
        if isinstance(mod, list):
            self._isotope_mods.extend(mod)
        else:
            self._isotope_mods.append(mod)

    def _add_labile_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._labile_mods is None:
            self._labile_mods = []
        if isinstance(mod, list):
            self._labile_mods.extend(mod)
        else:
            self._labile_mods.append(mod)

    def _add_unknown_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._unknown_mods is None:
            self._unknown_mods = []
        if isinstance(mod, list):
            self._unknown_mods.extend(mod)
        else:
            self._unknown_mods.append(mod)

    def _add_nterm_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._nterm_mods is None:
            self._nterm_mods = []
        if isinstance(mod, list):
            self._nterm_mods.extend(mod)
        else:
            self._nterm_mods.append(mod)

    def _add_cterm_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._cterm_mods is None:
            self._cterm_mods = []
        if isinstance(mod, list):
            self._cterm_mods.extend(mod)
        else:
            self._cterm_mods.append(mod)

    def _add_internal_mod(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._internal_mods is None:
            self._internal_mods = {}
        position = len(self._amino_acids) - 1
        if position not in self._internal_mods:
            self._internal_mods[position] = []
        if isinstance(mod, list):
            self._internal_mods[position].extend(mod)
        else:
            self._internal_mods[position].append(mod)

    def _add_interval(self, interval: Union[Interval, List[Interval]]) -> None:
        if self._intervals is None:
            self._intervals = []
        if isinstance(interval, list):
            self._intervals.extend(interval)
        else:
            self._intervals.append(interval)

    @_validate_single_mod_multiplier
    def _add_charge_adducts(self, mod: Union[Mod, List[Mod]]) -> None:
        if self._charge_adducts is None:
            self._charge_adducts = []
        if isinstance(mod, list):
            self._charge_adducts.extend(mod)
        else:
            self._charge_adducts.append(mod)

    def _parse_sequence_start(self) -> None:
        """
        Parse the start of the sequence, up until the first amino acid or interval
        """

        while not self._end_of_sequence():

            cur = self._current()
            if cur in AMINO_ACIDS or cur == '(':  # End of start sequence
                return None
            if cur == '[':  # N-term or unknown mods
                mods = self._parse_modifications('[', ']')
                next_char = self._parse_char()
                if next_char == '-':
                    self._add_nterm_mod(mods)
                elif next_char == '?':
                    self._add_unknown_mod(mods)
                else:
                    raise ProFormaFormatError(f"Expected '-' or '?, but got {cur}", self.position, self.sequence)
            elif cur == '<':  # Global mods
                for mod in self._parse_modifications('<', '>'):
                    if '@' in mod.val:  # Static mod
                        try:
                            self._add_static_mod(mod)
                        except ValueError as err:  # re-raise error with position, and sequence
                            raise ProFormaFormatError(err, self.position, self.sequence) from err
                    else:  # Isotope mod
                        try:
                            self._add_isotope_mod(mod)
                        except ValueError as err:  # re-raise error with position, and sequence
                            raise ProFormaFormatError(err, self.position, self.sequence) from err
            elif cur == '{':  # Labile mods
                self._add_labile_mod(self._parse_modification('{', '}'))
            else:
                raise ProFormaFormatError(r"Expected amino acid, '[', '{', or '<' but got " + cur, self.position,
                                          self.sequence)

    def _parse_sequence_middle(self) -> None:
        """
        Parse the middle of the sequence, up until the end of the sequence defined by '/' or '+'
        """
        dummy_interval = None
        while not self._end_of_sequence():
            cur = self._current()
            if cur in AMINO_ACIDS:  # Amino acid
                self._amino_acids.append(self._parse_char())
            elif cur == '[':  # mods for the previous amino acid
                self._add_internal_mod(self._parse_modifications('[', ']'))
            elif cur == '-':  # cterm mods (end of sequence)
                self._skip(1)
                self._add_cterm_mod(self._parse_modifications('[', ']'))
                return
            elif cur in ('/', '+'):  # charge ( end of sequence)
                return
            elif cur == '(':  # Interval start
                if dummy_interval is not None:
                    raise ProFormaFormatError("Overlapping intervals!", self.position, self.sequence)
                dummy_interval = [len(self._amino_acids), None, False, None]
                self._skip(1)
            elif cur == ')':  # Interval end
                if dummy_interval is None:
                    raise ProFormaFormatError("Interval ended without starting!", self.position, self.sequence)
                dummy_interval[1] = len(self._amino_acids)

                self._skip(1)
                if not self._end_of_sequence() and self._current() == '[':
                    dummy_interval[3] = self._parse_modifications('[', ']')

                self._add_interval(Interval(start=dummy_interval[0],
                                            end=dummy_interval[1],
                                            ambiguous=dummy_interval[2],
                                            mods=dummy_interval[3]))
                dummy_interval = None

            elif cur == '?':  # unknown mods
                if dummy_interval is None:
                    raise ProFormaFormatError("Unknown mod outside of interval", self.position, self.sequence)

                # interval is ambiguous
                dummy_interval[2] = True
                self._skip(1)
            else:
                raise ProFormaFormatError(
                    f"Expected either '[', '(', '?', '-', or '/' but got: {cur}", self.position, self.sequence)

    def _parse_sequence_end(self) -> None:
        """
        Parse the end of the sequence, up until the end of the sequence or start of the next sequence
        """
        while not self._end_of_sequence():
            cur = self._current()
            if cur == '/':  # charge
                self._skip(1)

                # check for // (crosslink)
                if not self._end_of_sequence() and self._current() == '/':
                    self._skip(1)
                    self._current_connection = True
                    return

                self._charge = self._parse_integer()

                # check for charge adducts
                if not self._end_of_sequence() and self._current() == '[':
                    self._add_charge_adducts(self._parse_modifications('[', ']'))

            elif cur == '+':  # next sequence
                self._skip(1)
                self._current_connection = False
                return
            else:
                raise ProFormaFormatError(
                    f"Invalid sequence: expected '/' or '+' but got {cur}", self.position, self.sequence)

    def _parse_char(self) -> str:
        # Assuming any character not a '[' or ']' is an amino acid for simplicity
        aa = self._current()
        self.position += 1
        return aa

    def _parse_modifications(self, opening_bracket='[', closing_bracket=']') -> List[Mod]:
        """
        Parses modifications from the sequence starting with the current position. The function will continue parsing
        until it reaches the end of the sequence or the there are no more sequential modifications.
        """
        mods = []
        while not self._end_of_sequence():
            if self._current() == opening_bracket:
                mod = self._parse_modification(opening_bracket, closing_bracket)
                mods.append(mod)
            else:
                break

        return mods

    def _parse_modification(self, opening_bracket='[', closing_bracket=']') -> Mod:
        """
        Parses a single modification from the sequence starting with the current position.
        """
        self.position += 1
        start = self.position
        bracket_depth = 1
        while not self._end_of_sequence() and bracket_depth > 0:
            if self.sequence[self.position] == opening_bracket:
                bracket_depth += 1
            elif self.sequence[self.position] == closing_bracket:
                bracket_depth -= 1
            self.position += 1

        if bracket_depth != 0:
            msg = f"Unmatched {opening_bracket} at position {self.position}"
            raise ProFormaFormatError(msg, self.position, self.sequence)

        mod = self.sequence[start:self.position - 1]

        multiplier = 1
        if not self._end_of_sequence() and self._peek() == '^':
            self.position += 1
            multiplier_start = self.position
            while not self._end_of_sequence() and self._peek().isdigit():
                self.position += 1
            multiplier = int(self.sequence[multiplier_start:self.position])

        return Mod(mod, multiplier)

    def _parse_integer(self) -> int:
        start = self.position
        digit_count = 0
        while not self._end_of_sequence():
            if self._peek().isdigit():
                digit_count += 1
                self.position += 1
            elif digit_count == 0 and self._peek() in ['+', '-']:
                self.position += 1
            else:
                break
        return int(self.sequence[start:self.position])

    def _current(self) -> str:
        return self.sequence[self.position]

    def _peek(self) -> str:
        return self.sequence[self.position] if not self._end_of_sequence() else None

    def _skip(self, n) -> None:
        self.position += n

    def _end_of_sequence(self) -> bool:
        return self.position >= self.length


def parse(sequence: str) -> Union[ProFormaAnnotation, MultiProFormaAnnotation]:
    """
    Parses a ProForma sequence string and returns its corresponding annotation object.

    Note that the function's behavior and the type of object returned depend on the structure of the input sequence.
    Single sequences result in ProFormaAnnotation objects, while multi-sequences result in MultiProFormaAnnotation
    objects.

    :param sequence: The sequence to parse.
    :type sequence: str

    :raises ProFormaFormatError: If the sequence is not valid.

    :return: Either a ProFormaAnnotation or a MultiProFormaAnnotation, based on the input
    :rtype: Union[ProFormaAnnotation, MultiProFormaAnnotation]

    .. python::

        Parsing a simple peptide sequence:
        >>> isinstance(parse('PEPTIDE'), ProFormaAnnotation)
        True

        Parsing a sequence with multiple peptides or complex modifications:
        >>> isinstance(parse('PEPTIDE+PEPTIDE'), MultiProFormaAnnotation)
        True


    """
    if _is_unmodified(sequence) is True:
        return ProFormaAnnotation(_sequence=sequence)

    annotations_connections = list(_ProFormaParser(sequence).parse())
    annotations = [annotation for annotation, _ in annotations_connections]
    annotations_connections = [connection for _, connection in annotations_connections]

    if len(annotations) == 1:
        return annotations[0]

    return MultiProFormaAnnotation(annotations, annotations_connections[:-1])


def serialize(annotation: Union[ProFormaAnnotation, MultiProFormaAnnotation], include_plus: bool = False) -> str:
    """
    Serializes a ProForma annotation or multiple ProForma annotations into a single string representation.

    :param annotation: Either a ProFormaAnnotation or a MultiProFormaAnnotation.
    :type annotation: Union[ProFormaAnnotation, MultiProFormaAnnotation]

    :return: A string representation of the ProForma annotation.
    :rtype: str

    . python::

        Serializing a simple ProForma annotation:
        >>> serialize(ProFormaAnnotation(_sequence='PEPTIDE'))
        'PEPTIDE'

        >>> pfa1 = ProFormaAnnotation(_sequence='PEPTIDE')
        >>> pfa2 = ProFormaAnnotation(_sequence='PEPTIDE')

        Serializing a MultiProFormaAnnotation with chimeric connections:
        >>> multi_annotation = MultiProFormaAnnotation([pfa1, pfa2], [False])
        >>> serialize(multi_annotation)
        'PEPTIDE+PEPTIDE'

        Serializing a MultiProFormaAnnotation with crosslink connections:
        >>> multi_annotation = MultiProFormaAnnotation([pfa1, pfa2], [True])
        >>> p = serialize(multi_annotation)
        >>> p == r'PEPTIDE\\\PEPTIDE'
        True

    """

    return annotation.serialize(include_plus)


def _serialize_annotation_start(annotation: ProFormaAnnotation, include_plus: bool) -> str:
    comps = []

    # add labile mods
    if annotation.has_labile_mods():
        for mod in annotation.labile_mods:
            comps.append(mod.serialize('{}', include_plus))

    if annotation.has_static_mods():
        for mod in annotation.static_mods:
            comps.append(mod.serialize('<>', include_plus))

    # Add global mods
    if annotation.has_isotope_mods():
        for mod in annotation.isotope_mods:
            comps.append(mod.serialize('<>',include_plus))

    # Unknown mods
    if annotation.has_unknown_mods():
        for mod in annotation.unknown_mods:
            comps.append(mod.serialize('[]',include_plus))
        comps.append("?")

    # N-term mods
    if annotation.has_nterm_mods():
        for mod in annotation.nterm_mods:
            comps.append(mod.serialize('[]', include_plus))
        comps.append("-")

    return ''.join(comps)


def _serialize_annotation_middle(annotation: ProFormaAnnotation, include_plus: bool) -> str:
    comps = []
    # Sequence
    for i, aa in enumerate(annotation.sequence):

        if annotation.intervals:
            for interval in annotation.intervals:
                if interval.start == i:
                    comps.append("(")
                    if interval.ambiguous:
                        comps.append("?")
                if interval.end == i:
                    comps.append(")")

                    if interval.mods:
                        for mod in interval.mods:
                            comps.append(mod.serialize('[]', include_plus))

        comps.append(aa)

        # Internal mods
        if annotation.internal_mods and i in annotation.internal_mods:
            for mod in annotation.internal_mods[i]:
                comps.append(mod.serialize('[]',include_plus))

    # add end interval
    i = len(annotation.sequence)
    if annotation.intervals:
        for interval in annotation.intervals:
            if interval.end == i:
                comps.append(")")
                if interval.mods:
                    for mod in interval.mods:
                        comps.append(mod.serialize('[]',include_plus))

    return ''.join(comps)


def _serialize_annotation_end(annotation: ProFormaAnnotation, include_plus: bool) -> str:
    comps = []
    # C-term mods
    if annotation.cterm_mods:
        comps.append('-')
        for mod in annotation.cterm_mods:
            comps.append(mod.serialize('[]', include_plus))

    # Charge
    if annotation.charge:
        comps.append(f"/{annotation.charge}")

    if annotation.charge_adducts:
        for mod in annotation.charge_adducts:
            comps.append(mod.serialize('[]', include_plus))

    return ''.join(comps)


def _serialize_annotation(annotation: ProFormaAnnotation, include_plus: bool) -> str:
    return _serialize_annotation_start(annotation, include_plus) + \
        _serialize_annotation_middle(annotation, include_plus) + \
        _serialize_annotation_end(annotation, include_plus)
