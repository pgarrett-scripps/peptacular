from functools import wraps
from typing import *

import regex as re

from ..constants import ADDUCT_PATTERN, ISOTOPE_NUM_PATTERN
from ..proforma_dataclasses import ModValue
from ..proforma_dataclasses import ChemComposition
from ..proforma_dataclasses import (
    Mod,
)


def _parse_modifications(
    proforma_sequence: str, opening_bracket: str = "[", closing_bracket: str = "]"
) -> List[Mod]:
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
            mod, new_position = _parse_modification(
                proforma_sequence[position:], opening_bracket, closing_bracket
            )
            mods.append(mod)
            position += new_position  # Update position based on the length of the parsed modification.
        else:
            position += 1  # Move to the next character if it's not an opening bracket.
    return mods


def _parse_modification(
    proforma_sequence: str, opening_bracket: str = "[", closing_bracket: str = "]"
) -> Tuple[Mod, int]:
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
    if position < len(proforma_sequence) and proforma_sequence[position] == "^":
        multiplier_start = position + 1
        multiplier, parsed_length = _parse_integer(proforma_sequence[multiplier_start:])
        position += (
            parsed_length + 1
        )  # Account for the '^' and the length of the integer parsed.

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
        elif digit_count == 0 and proforma_sequence[i] in ["+", "-"]:
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
        raise TypeError(
            f"Invalid type for charge adducts: {type(mod)}! Mod or Mod.val must be of type string."
        )

    charge_adducts = {}

    mods = mod.split(",")
    for m in mods:
        if isinstance(m, Mod):
            mod_str = m.val
        else:
            mod_str = m

        for sign, count, ion in ADDUCT_PATTERN.findall(mod_str):
            if count == "":
                count = 1
            else:
                count = int(count)

            if sign == "-":
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
        sign = "+" if count >= 0 else ""
        sign = "-" if count < 0 else sign
        count_abs = abs(count) if abs(count) != 1 else ""
        adducts_list.append(f"{sign}{count_abs}{ion}")

    return Mod(",".join(adducts_list), 1)


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
            raise TypeError(
                f"Invalid type for static mods: {type(mod)}! Mod or Mod.val must be of type string."
            )

        mod_info, residues = mod.split("@")
        residues = residues.split(",")
        static_mods = _parse_modifications(mod_info, "[", "]")

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

        >>> write_static_mods({'C': [Mod(val=1,multiple=1)], 'S': [Mod(val=3.1415,multiple=1)], 'D': [Mod(val=3.1415,multiple=1)]})
        [Mod('[1]@C', 1), Mod('[3.1415]@S,D', 1)]

        >>> write_static_mods({'C': [Mod(val='Formula:[13C]H20', multiple=1)]})
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
        mod_str = "".join(mod.serialize("[]") for mod in mod_tuple)
        static_mods.append(Mod(f"{mod_str}@{residue_str}", 1))

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
            raise TypeError(
                f"Invalid type for isotope mods: {type(mod)}! Mod or Mod.val must be of type string."
            )

        # remove digits
        base_aa = re.sub(ISOTOPE_NUM_PATTERN, "", mod)
        isotope_map[base_aa] = mod

    # If any keys are D or T, then replace them with H
    if "D" in isotope_map:
        isotope_map["H"] = isotope_map.pop("D")

    if "T" in isotope_map:
        isotope_map["H"] = isotope_map.pop("T")

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
                    f"Invalid multiplier {mod.mult} for mod {mod.val}. Multiplier must not be greater than 1.",
                    None,
                    None,
                )
        else:
            for m in mod:
                if m.mult > 1:
                    raise ValueError(
                        f"Invalid multiplier {m.mult} for mod {m.val}. Multiplier must not be greater than 1.",
                        None,
                        None,
                    )

        return func(self, mod)

    return wrapper


"""
This section is largely based on the code from the `biopython` library, specifically the `Bio.SeqUtils` module.
"""

from .property_data import *


def get_aa_value(
    aa: str,
    aa_data: Dict[str, float],
    missing_aa_handling: Literal[
        "zero", "avg", "min", "max", "median", "error", "skip"
    ],
    weighting_scheme: float = 1,
    normalize: bool = False,
) -> Optional[float]:
    if missing_aa_handling == "zero":
        default_value = 0.0
    elif missing_aa_handling == "avg":
        default_value = sum(aa_data.values()) / len(aa_data)
    elif missing_aa_handling == "min":
        default_value = min(aa_data.values())
    elif missing_aa_handling == "max":
        default_value = max(aa_data.values())
    elif missing_aa_handling == "median":
        sorted_values = sorted(aa_data.values())
        mid = len(sorted_values) // 2
        if len(sorted_values) % 2 == 0:
            default_value = (sorted_values[mid - 1] + sorted_values[mid]) / 2
        else:
            default_value = sorted_values[mid]
    elif missing_aa_handling == "error":
        default_value = "error"
    elif missing_aa_handling == "skip":
        return "skip"
    else:
        raise ValueError(
            f"Invalid default value: {missing_aa_handling}. Choose from 'zero', 'avg', 'min', 'max', 'median', 'error', or 'skip'."
        )

    # B -> Aspartic acid or Asparagine
    # J -> Leucine or Isoleucine
    # Z -> Glutamic acid or Glutamine
    # X -> Any amino acid (not a valid single-letter code)

    if aa in aa_data:
        if normalize:
            # Normalize the value to a range of 0-1
            min_value = min(aa_data.values())
            max_value = max(aa_data.values())
            return (
                (aa_data[aa] - min_value) / (max_value - min_value)
            ) * weighting_scheme
        else:
            # Return the raw value multiplied by the weight
            return aa_data[aa] * weighting_scheme
    elif aa == "B":
        # average of Aspartic acid and Asparagine
        return (
            (
                get_aa_value("D", aa_data, missing_aa_handling)
                + get_aa_value("N", aa_data, missing_aa_handling)
            )
            / 2
        ) * weighting_scheme
    elif aa == "J":
        # average of Leucine and Isoleucine
        return (
            (
                get_aa_value("L", aa_data, missing_aa_handling)
                + get_aa_value("I", aa_data, missing_aa_handling)
            )
            / 2
        ) * weighting_scheme
    elif aa == "Z":
        # average of Glutamic acid and Glutamine
        return (
            (
                get_aa_value("E", aa_data, missing_aa_handling)
                + get_aa_value("Q", aa_data, missing_aa_handling)
            )
            / 2
        ) * weighting_scheme
    elif aa == "X":
        # Mean of all
        if normalize:
            # Normalize the average value to a range of 0-1
            min_value = min(aa_data.values())
            max_value = max(aa_data.values())
            return (
                (sum(aa_data.values()) / len(aa_data) - min_value)
                / (max_value - min_value)
            ) * weighting_scheme
        else:
            # Return the average value multiplied by the weight
            return (sum(aa_data.values()) / len(aa_data)) * weighting_scheme
    else:
        if default_value == "error":
            raise ValueError(
                f"Invalid amino acid: {aa}. No hydrophobicity value found."
            )
        elif default_value == "skip":
            return None
        return default_value * weighting_scheme
