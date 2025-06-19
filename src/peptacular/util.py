"""
Utils.py
"""

import warnings
from typing import Union, List, Tuple, Dict, Generator, Optional
import regex
from functools import wraps
import regex as re
from typing import Callable, Dict, List, Optional, Tuple, Union

from .constants import ADDUCT_PATTERN, ISOTOPE_NUM_PATTERN
from .types import ChemComposition, ModValue
from .proforma_dataclasses import (
    Mod,
)


def get_regex_match_indices(
    input_str: str, regex_str: Union[str, regex.Pattern], offset: int = 0
) -> Generator[int, None, None]:
    """
    Identify the starting indexes of occurrences of a given regex pattern within a string.

    :param input_str: The sequence in which to search.
    :type input_str: str
    :param regex_str: The regex pattern to search for.
    :type regex_str: str
    :param offset: An optional offset to add to each identified index. Default is 0.
    :type offset: int
    :return: A list of starting indexes where the regex pattern is found in the sequence.
    :rtype: List[int]

    .. code-block:: python

        >>> list(get_regex_match_indices("PEPTIDE", "P"))
        [1, 3]

        >>> list(get_regex_match_indices("PEPTIDE", regex.compile("E")))
        [2, 7]

        >>> list(get_regex_match_indices("PEPTIDE", 'E'))
        [2, 7]

        # More complex regex
        >>> list(get_regex_match_indices("PEPTIDE", "P[ST]"))
        [3]

        >>> list(get_regex_match_indices("PPPP", "PP"))
        [1, 2, 3]

        >>> list(get_regex_match_indices("PEPCTIDE", "(?=C)"))
        [3]

        >>> list(get_regex_match_indices("PEPCTIDE", "C"))
        [4]

        >>> list(get_regex_match_indices("PEPCTIDE", "(?<=C)"))
        [4]

    """

    if not isinstance(regex_str, regex.Pattern):
        regex_pattern = regex.compile(regex_str)
    else:
        regex_pattern = regex_str

    for match in regex_pattern.finditer(input_str, overlapped=True):
        if match.start() != match.end():
            warnings.warn(
                message="The regex pattern has a match with a none zero length. Using start index + 1 for the match.",
            )
            yield match.start() + offset + 1
        else:
            yield match.start() + offset


def get_regex_match_range(
    input_str: str, regex_str: Union[str, regex.Pattern], offset: int = 0
) -> List[Tuple[int, int]]:
    """
    Identify the starting indexes of occurrences of a given regex pattern within a string.

    :param input_str: The sequence in which to search.
    :type input_str: str
    :param regex_str: The regex pattern to search for.
    :type regex_str: str
    :param offset: An optional offset to add to each identified index. Default is 0.
    :type offset: int
    :return: A list of starting indexes where the regex pattern is found in the sequence.
    :rtype: List[int]

    .. code-block:: python

        >>> get_regex_match_range("PEPTIDE", "P")
        [(0, 1), (2, 3)]

        >>> get_regex_match_range("PEPTIDE", regex.compile("E"))
        [(1, 2), (6, 7)]

        # More complex regex
        >>> get_regex_match_range("PEPTIDE", "P[ST]")
        [(2, 4)]

        >>> get_regex_match_range("PPPP", "PP")
        [(0, 2), (1, 3), (2, 4)]

    """

    # if regex is compiled
    if isinstance(regex_str, regex.Pattern):
        return [
            (i.start() + offset, i.end() + offset)
            for i in regex_str.finditer(input_str)
        ]

    return [
        (match.start() + offset, match.end() + offset)
        for match in regex.finditer(regex_str, input_str, overlapped=True)
    ]


def _validate_span(span: Tuple[int, int, int]) -> None:
    """
    Validates if a given span is valid.

    :param span: A tuple representing the span.
    :type span: Tuple[int, int, int]
    :raises ValueError: If the span is not valid.

    .. code-block:: python

        >>> _validate_span((0, 5, 0))  # No error raised

        >>> _validate_span((0, 0, 0))  # No error raised

        >>> _validate_span((5, 0, 0))  # Raises ValueError
        Traceback (most recent call last):
        ValueError: Start of span: 5, should be less than or equal to end of span: 0.

        >>> _validate_span((-1, 0, 0))  # Raises ValueError
        Traceback (most recent call last):
        ValueError: Start of span should be non-negative, got -1.

        >>> _validate_span((0, -1, 0))  # Raises ValueError
        Traceback (most recent call last):
        ValueError: End of span should be non-negative, got -1.

    """

    start, end, _ = span
    if start < 0:
        raise ValueError(f"Start of span should be non-negative, got {start}.")
    if end < 0:
        raise ValueError(f"End of span should be non-negative, got {end}.")
    if start > end:
        raise ValueError(
            f"Start of span: {start}, should be less than or equal to end of span: {end}."
        )


def _construct_ambiguity_intervals(
    counts: List[int], reverse: bool
) -> List[Tuple[int, int]]:
    """
    Construct intervals for sequences of zeros in the counts list. When reverse is false, start from the left hand side
    and move to the right. When reverse is true, start from the right hand side and move to the left. Intervals start
    at 0 and end on any positive value. Both are inclusive. Returned intervals should be in forwards format, that is
    have a starting value less than the ending value.

    :param counts: List of integers (typically counts)
    :param reverse: If True, reverse the list before processing
    :return: List of intervals [start, end] indicating runs of zeros

    .. code-block:: python

        # [0, 1, 1, 1, 0, 0, 0]
        # [1, 1, 0, 0, 1, 1, 1] # ambiguity
        >>> _construct_ambiguity_intervals([0, 1, 1, 1, 0, 0, 0], reverse=False)
        [(0, 1), (4, 6)]

        # [0, 0, 1, 1, 1, 1, 0]
        # [1, 1, 0, 0, 0, 1, 1] # ambiguity
        >>> _construct_ambiguity_intervals([0, 0, 1, 1, 1, 0, 0], reverse=True)
        [(0, 1), (4, 6)]

        >>> _construct_ambiguity_intervals([0, 1, 1, 1, 0, 0, 1], reverse=False)
        [(0, 1), (4, 6)]
    """

    if reverse:
        ambiguity_intervals = _construct_ambiguity_intervals(
            counts[::-1], reverse=False
        )
        ambiguity_intervals = [
            (len(counts) - 1 - end, len(counts) - 1 - start)
            for start, end in ambiguity_intervals
        ]
        # sort the intervals
        ambiguity_intervals.sort(key=lambda x: x[0])
        return ambiguity_intervals

    ambiguity_intervals = []
    current_interval = None
    for i, cnt in enumerate(counts):
        if cnt == 0:
            if current_interval is not None:
                current_interval = (current_interval[0], i)
            else:
                current_interval = (i, i)

        else:
            if current_interval is not None:
                current_interval = (current_interval[0], i)
                ambiguity_intervals.append(current_interval)
                current_interval = None
            else:
                continue

    if current_interval is not None:
        current_interval = (current_interval[0], len(counts) - 1)
        ambiguity_intervals.append(current_interval)

    return ambiguity_intervals


def _combine_ambiguity_intervals(
    *interval_lists: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """
    Merge multiple lists of ambiguity intervals into a single list of common ambiguity intervals.

    This function identifies positions that are ambiguous across all provided interval lists.
    For a position to be considered ambiguous in the result, it must be contained in at least
    one interval from each input list. The function then constructs optimized intervals
    covering these common ambiguous positions.

    Intervals are represented as tuples (start, end) where:
    - start is inclusive
    - end is exclusive

    Intervals with identical start and end values (zero-length intervals) are removed.

    :param interval_lists: Variable number of lists containing ambiguity intervals
    :type interval_lists: List[Tuple[int, int]]

    :return: A list of merged intervals representing positions that are ambiguous across all input lists
    :rtype: List[Tuple[int, int]]

    .. code-block:: python

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0,1)])
        [(0, 1)]

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0,1), (4,5)])
        [(0, 1), (4, 5)]

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0, 4), (5, 6)])
        [(0, 1), (5, 6)]

        >>> _combine_ambiguity_intervals([(2, 5)], [(3, 6)])
        [(3, 5)]

        >>> _combine_ambiguity_intervals([(0, 1)], [(4, 6)])
        []
    """

    # First, collect all unique intervals from input lists
    all_intervals = set()
    for interval_list in interval_lists:
        for interval in interval_list:
            all_intervals.add(interval)

    # Remove intervals where start == end (these are not ambiguous)
    filtered_intervals = {(start, end) for start, end in all_intervals if start != end}

    # Find all possible indices that are covered by any interval
    all_indices = set()
    for start, end in filtered_intervals:
        for i in range(start, end):
            all_indices.add(i)

    # For each index, check if it's contained in at least one interval from each input list
    common_indices = set()
    for idx in all_indices:
        is_common = True
        for interval_list in interval_lists:
            if not any(start <= idx < end for start, end in interval_list):
                is_common = False
                break
        if is_common:
            common_indices.add(idx)

    # If no common indices found, return empty list
    if not common_indices:
        return []

    # Construct new intervals from the common indices
    result = []
    if common_indices:
        sorted_indices = sorted(common_indices)
        start = sorted_indices[0]
        for i in range(1, len(sorted_indices)):
            if sorted_indices[i] > sorted_indices[i - 1] + 1:
                # Gap found, close the current interval and start a new one
                result.append((start, sorted_indices[i - 1] + 1))
                start = sorted_indices[i]
        # Add the last interval
        result.append((start, sorted_indices[-1] + 1))

    return result


def _get_mass_shift_interval(
    forward_coverage: List[int], reverse_coverage: List[int]
) -> Optional[Tuple[int, int]]:
    """
        Determine the interval where a mass shift should be placed based on fragment ion coverage.

    This function examines the forward and reverse ion coverage to identify the region where
    a mass shift (such as a modification) should be positioned. It returns the start and end
    indices (inclusive) of this region, or None if no suitable region is found.

    The mass shift interval is determined by:
    1. Finding the highest position with forward ion coverage
    2. Finding the lowest position with reverse ion coverage
    3. The mass shift belongs between these two positions

    :param forward_coverage: Binary list indicating forward ion coverage (1) or no coverage (0)
    :type forward_coverage: List[int]
    :param reverse_coverage: Binary list indicating reverse ion coverage (1) or no coverage (0)
    :type reverse_coverage: List[int]

    :return: A tuple containing the start and end indices (inclusive) for the mass shift,
             or None if there is no valid interval
    :rtype: Optional[Tuple[int, int]]

    .. code-block:: python

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,0,1,1,1])
        (3, 3)

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,1,1,1,1])
        (3, 3)

        >>> _get_mass_shift_interval([1,1,0,0,0,0,0], [0,0,0,0,1,1,1])
        (2, 3)

        >>> _get_mass_shift_interval([0,0,0,0,0,0,0], [0,0,0,0,1,1,1])
        (0, 3)

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,0,0,0,0])
        (3, 6)

        >>> _get_mass_shift_interval([1,1,1,1,1,0,0], [0,0,0,0,1,1,1]) # None

    """

    highest_forward_fragment = [i for i, cnt in enumerate(forward_coverage) if cnt > 0]
    if len(highest_forward_fragment) == 0:
        highest_forward_fragment = -1
    else:
        highest_forward_fragment = max(highest_forward_fragment)

    highest_reverse_fragment = [i for i, cnt in enumerate(reverse_coverage) if cnt > 0]
    if len(highest_reverse_fragment) == 0:
        highest_reverse_fragment = len(reverse_coverage)
    else:
        highest_reverse_fragment = min(highest_reverse_fragment)

    if highest_forward_fragment >= highest_reverse_fragment:
        return None

    if highest_forward_fragment == highest_reverse_fragment - 1:
        return (highest_forward_fragment + 1, highest_forward_fragment + 1)

    # if there is a gap between the two, return the gap
    return (highest_forward_fragment + 1, highest_reverse_fragment - 1)


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

    count_str = ""

    charge = 1
    for i, c in enumerate(ion):
        if c == "-":
            charge = -1
        elif c == "+":
            pass
        elif c.isdigit():
            count_str += c
        else:
            cnt = int(count_str) if count_str else 1
            return cnt * charge, ion[i:]

    raise ValueError(f"Bad Ion Count: {ion}")


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

    symbol = ""
    for i, c in enumerate(ion):
        if c.isdigit() or c in ["+", "-"]:
            return symbol, ion[i:]
        symbol += c
    return symbol, ""


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
    count_str, sign = "", 1
    for c in ion:
        if c == "-":
            sign = -1
        elif c == "+":
            continue
        else:
            count_str += c
    cnt = int(count_str) if count_str else 1
    return cnt * sign, ""


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
