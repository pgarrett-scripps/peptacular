import re
from copy import deepcopy
from typing import Dict, List, Any

import regex as reg

from .constants import PROTEASES
from .spans import span_to_sequence, build_spans, build_non_enzymatic_spans, build_left_semi_spans, \
    build_right_semi_spans


def _are_parentheses_balanced(text):
    """
    Check if parentheses in the given text are balanced or not.

    The function uses a stack data structure to ensure each opening parenthesis '('
    has a corresponding closing parenthesis ')' and vice versa.

    Args:
        text (str): The input string to check for balanced parentheses.

    Returns:
        bool: True if parentheses are balanced, False otherwise.
    """
    stack = []
    for i, char in enumerate(text):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if not stack:
                return False
            else:
                stack.pop()
    return len(stack) == 0


def convert_to_ms2_pip_style(sequence: str) -> str:
    """
    Converts a peptide sequence to the format used by MS2PIP prediction software.

    This function parses a peptide sequence and transforms it into a string representation
    where each modification's location is paired with the type of modification, separated
    by a pipe ('|') symbol.

    Args:
        sequence (str): The peptide sequence, potentially including modifications.

    Returns:
        str: The sequence transformed into MS2PIP format.

    Example:
        >>> convert_to_ms2_pip_style('PEP(Phospho)TIDE(Phospho)K')
        '2|Phospho|6|Phospho'
    """

    # Parse modifications from the sequence
    modification_dict = parse_modified_sequence(sequence)

    ms2_pip_pairs = []
    for location, modification in modification_dict.items():
        # Adjust location index for MS2PIP
        ms2_pip_location = location + 1
        if ms2_pip_location == len(strip_modifications(sequence)):
            ms2_pip_location = -1

        # Append the new location and modification pair
        ms2_pip_pairs.append(f'{ms2_pip_location}|{modification}')

    # Combine pairs into a single string, separated by '|'
    return '|'.join(ms2_pip_pairs)


def parse_modified_sequence(sequence: str) -> Dict[int, str]:
    """
    Parses a peptide sequence with modifications and returns a dictionary where keys
    represent the position of the modified amino acid and values are the respective modifications.

    Args:
        sequence (str): The peptide sequence, potentially including modifications.

    Returns:
        Dict[int, str]: A dictionary with the modification values indexed by the position
                        of the modified amino acid.

    Raises:
        ValueError: If the sequence contains unmatched parentheses indicating incorrect
                    modification notation.

    Example:
        >>> parse_modified_sequence('PEP(Phospho)TIDEK')
        {3: 'Phospho'}
    """

    # Check parentheses balance
    if not _are_parentheses_balanced(sequence):
        raise ValueError(f'Incorrect modification notation in peptide sequence : {sequence}!')

    mod_pattern = re.compile(r'\(([^)]*)\)')
    matches = mod_pattern.finditer(sequence)

    modifications = {}
    mod_offset = 0

    for match in matches:
        modification_position = match.start() - mod_offset - 1
        modification_value = match.group(1)
        modifications[modification_position] = modification_value

        # Compute offset for the next modification
        mod_offset += len(match.group())

    return modifications


def create_modified_sequence(sequence: str, modifications: Dict[int, Any]) -> str:
    """
    Generates a peptide sequence with modifications based on an unmodified sequence and a dictionary of modifications.

    The modifications are provided as a dictionary where keys correspond to the indices of the amino acids
    in the sequence to be modified, and values are the modifications to be applied. Modifications are inserted
    into the sequence in descending order of their indices to ensure their validity after each modification.

    Args:
        sequence (str): Unmodified peptide sequence.
        modifications (Dict[int, str]): Dictionary with indices of the modified amino acids and corresponding
                                        modifications.

    Returns:
        str: The peptide sequence with the specified modifications.

    Raises:
        ValueError: If an index in modifications is invalid for the sequence.

    Example:
        >>> create_modified_sequence('PEPTIDEK', {3: 5.0})
        'PEP(5.0)TIDEK'
    """

    modified_sequence = []
    prev_index = 0

    # Sort the modifications by index in descending order
    for i, mod in sorted(modifications.items()):

        if i < -1 or i >= len(sequence):
            raise ValueError(f'Index of modification: {i} is invalid for peptide sequence: {sequence}')

        # Insert the modification into the modified peptide sequence
        modified_sequence.append(sequence[prev_index: i + 1])
        modified_sequence.append(f"({mod})")
        prev_index = i + 1

    modified_sequence.append(sequence[prev_index:])
    return ''.join(modified_sequence)


def strip_modifications(sequence: str) -> str:
    """
    Remove any non-amino-acid characters from the given peptide sequence. This function specifically removes
    the modifications which are generally annotated within parentheses.

    Args:
        sequence (str): The peptide sequence to be processed.

    Raises:
        ValueError: If the parentheses in the given sequence are unbalanced.

    Returns:
        str: The peptide sequence without any modification annotations.
    """
    if _are_parentheses_balanced(sequence) is False:
        raise ValueError(f'Incorrect modification notation in peptide sequence : {sequence}!')

    unmodified_sequence = re.sub(r'\([^)]*\)', '', sequence)
    return unmodified_sequence


def get_left_semi_sequences(sequence: str, min_len: int = 1, max_len: int = None) -> List[str]:
    """
    Returns a set of left substrings of string `sequence` that have lengths between `min_len` and `max_len`.

    Args:
        sequence (str): The string to generate subsequences from.
        min_len (int, optional): The minimum length of the subsequence. Defaults to 1.
        max_len (int, optional): The maximum length of the subsequence. If not provided, it defaults to the length of
                                 the sequence minus one.

    Returns:
        set: A list of left subsequences of the input sequence.
    """
    spans = build_left_semi_spans((0, len(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def get_right_semi_sequences(sequence: str, min_len: int = 1, max_len: int = None) -> List[str]:
    """
    Returns a set of right substrings of string `sequence` that have lengths between `min_len` and `max_len`.

    Args:
        sequence (str): The string to generate subsequences from.
        min_len (int, optional): The minimum length of the subsequence. Defaults to 1.
        max_len (int, optional): The maximum length of the subsequence. If not provided, it defaults to the length of
                                 the sequence minus one.

    Returns:
        set: A list of right subsequences of the input sequence.
    """

    spans = build_right_semi_spans((0, len(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


from typing import List

def get_semi_sequences(sequence: str, min_len: int = None, max_len: int = None) -> List[str]:
    """
    Generate a list of all semi-enzymatic amino acid sequences of a given string sequence.

    Args:
        sequence (str): A string representing the protein sequence.
        min_len (int, optional): The minimum length of the semi-enzymatic sequence. Default is None.
        max_len (int, optional): The maximum length of the semi-enzymatic sequence. Default is None.

    Returns:
        List[str]: A list of semi-enzymatic sequences within the defined length bounds.
    """

    return get_left_semi_sequences(sequence, min_len, max_len) + get_right_semi_sequences(sequence, min_len, max_len)


def get_non_enzymatic_sequences(sequence: str, min_len: int = None, max_len: int = None) -> List[str]:
    """
    Generate a list of all non-enzymatic amino acid sequences of a given string sequence.

    Args:
        sequence (str): A string representing the protein sequence.
        min_len (int, optional): The minimum length of the non-enzymatic sequence. Default is None.
        max_len (int, optional): The maximum length of the non-enzymatic sequence. Default is None.

    Returns:
        List[str]: A list of non-enzymatic sequences within the defined length bounds.
    """

    spans = build_non_enzymatic_spans((0, len(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def get_enzymatic_sequences(sequence: str, enzyme_regex: str, missed_cleavages: int, semi: bool, min_len: int = None,
                            max_len: int = None) -> List[str]:
    """
    Generate a list of all enzymatic amino acid sequences of a given string sequence based on provided enzyme_regex.

    Args:
        sequence (str): A string representing the protein sequence.
        enzyme_regex (str): The enzyme regular expression for identifying cleavage sites.
        missed_cleavages (int): The maximum number of missed cleavages.
        semi (bool): Whether the sequence includes semi enzymatic peptides.
        min_len (int, optional): The minimum length of the enzymatic sequence. Default is None.
        max_len (int, optional): The maximum length of the enzymatic sequence. Default is None.

    Returns:
        List[str]: A list of enzymatic sequences within the defined length bounds and satisfying the enzyme cleavage rules.
    """

    cleavage_sites = identify_cleavage_sites(sequence, enzyme_regex)
    spans = build_spans(len(sequence), cleavage_sites, missed_cleavages, min_len, max_len, semi)
    return [span_to_sequence(sequence, span) for span in spans]


def add_static_mods(sequence: str, mod_map: Dict[str, Any]):
    """
    Add static modifications to an amino acid sequence.

    Args:
        sequence (str): Original amino acid sequence.
        mod_map (Dict[str, float]): Dictionary mapping amino acids to the mass of their modifications.

    Returns:
        str: Modified amino acid sequence.
    """

    stripped_sequence = strip_modifications(sequence)
    sequence_mod_map = parse_modified_sequence(sequence)
    for aa, mass in mod_map.items():
        indexes = [i for i, x in enumerate(stripped_sequence) if x == aa]
        for index in indexes:
            sequence_mod_map[index] = mass

    return create_modified_sequence(stripped_sequence, sequence_mod_map)


def rec_mod_builder(mods: Dict[str, float], sequence: str, i: int, mod_dict: Dict[int, float], max_mod_count: int):
    """
    Recursive function to generate all possible combinations of modifications on an amino acid sequence.

    Args:
        mods (Dict[str, float]): Dictionary mapping amino acids to the mass of their modifications.
        sequence (str): amino acid sequence to modify.
        i (int): Current index on the amino acid sequence.
        mod_dict (Dict[int, float]): Dictionary mapping the indices of the amino acids to the mass of their
        modifications.
        max_mod_count (int): Maximum number of modifications allowed on the amino acid sequence.

    Yields:
        Dict[int, float]: Modified dictionary mapping the indices of the amino acid sequence to the mass of their
        modifications.
    """

    if i == len(sequence):
        yield mod_dict
        return

    original_mod_dict = deepcopy(mod_dict)
    for aa, mass in mods.items():

        if len(mod_dict) >= max_mod_count:
            break

        if sequence[i] == aa and i not in mod_dict:
            mod_dict[i] = mass
            yield from rec_mod_builder(mods, sequence, i + 1, mod_dict, max_mod_count)

    yield from rec_mod_builder(mods, sequence, i + 1, original_mod_dict, max_mod_count)


def add_variable_mods(sequence: str, mod_map: Dict[str, float], max_mods: int) -> List[str]:
    """
    Add variable modifications to an amino acid sequence.

    Args:
        sequence (str): Original amino acid sequence.
        mod_map (Dict[str, float]): Dictionary mapping amino acids to the mass of their modifications.
        max_mods (int): Maximum number of modifications allowed on the peptide.

    Returns:
        List[str]: List of all possible modified peptide sequences, each represented as a string.
    """

    original_mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)
    mods = rec_mod_builder(mod_map, stripped_sequence, 0, original_mods, max_mods + len(original_mods))
    return [create_modified_sequence(stripped_sequence, mod) for mod in mods]


def sequence_generator(sequence: str, forward: bool):
    """
    Generates amino acid sequences either from the front or the back.

    Args:
        sequence (str): The initial amino acid sequence.
        forward (bool): If True, start generating from the front; else, start from the back.

    Returns:
        generator: A generator yielding amino acid sequences.
    """
    if forward:
        return _sequence_generator_front(sequence)
    else:
        return _sequence_generator_back(sequence)


def _sequence_generator_back(sequence: str):
    """
    Generates amino acid sequences starting from the back of the sequence.

    Args:
        sequence (str): The initial amino acid sequence.

    Yields:
        str: The sequence with the back amino acid (and modification) removed.

    Returns:
        None: If the modification starts at the beginning of the sequence.
    """
    while sequence:
        index = 0
        if sequence[-1] == ')':
            start_of_modification = sequence.rfind('(')

            if start_of_modification == 0:
                return None

            index = start_of_modification

        yield sequence
        sequence = sequence[:index - 1]


def _sequence_generator_front(sequence: str):
    """
    Generates sub-sequences starting from the front of the sequence.

    Args:
        sequence (str): The initial sequence.

    Yields:
        str: The sequence with the front amino acid (and modification) removed.

    Returns:
        None: If the modification ends at the end of the sequence.
    """
    if sequence[0] == '(':
        yield sequence

        end_of_modification = sequence.find(')')
        sequence = sequence[end_of_modification + 2:]

        if sequence[0] == '(':
            end_of_modification = sequence.find(')')
            sequence = sequence[end_of_modification + 1:]

    while sequence:
        index = 0
        if len(sequence) > 1 and sequence[1] == '(':
            end_of_modification = sequence.find(')')

            if end_of_modification == len(sequence):
                return None

            index = end_of_modification

        yield sequence
        sequence = sequence[index + 1:]


def reverse_sequence(sequence: str):
    """
    Reverses the sequence, while preserving the position of any modifications.

    Args:
        sequence (str): The amino acid sequence to be reversed.

    Returns:
        str: The reversed sequence with modifications preserved.
    """
    mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)[::-1]
    mod_reverse = {len(stripped_sequence) - 1 - k if k != -1 else -1: v for k, v in mods.items()}
    return create_modified_sequence(stripped_sequence, mod_reverse)


def shift_sequence_left(sequence: str, shift: int):
    """
    Shifts the sequence to the left by a given number of positions, while preserving the position of any modifications.

    Args:
        sequence (str): The sequence to be shifted.
        shift (int): The number of positions to shift the sequence to the left.

    Returns:
        str: The shifted sequence with modifications preserved.
    """
    mods = parse_modified_sequence(sequence)
    stripped_sequence = strip_modifications(sequence)
    seq_len = len(stripped_sequence)

    # Take modulus to ensure proper wrapping around
    effective_shift = shift % seq_len

    shifted_stripped_sequence = stripped_sequence[effective_shift:] + stripped_sequence[:effective_shift]

    # Shift mod positions, taking care to avoid negative keys
    shifted_sequence_mods = {(k - effective_shift) % seq_len if k != -1 else -1: v for k, v in mods.items()}

    return create_modified_sequence(shifted_stripped_sequence, shifted_sequence_mods)


def identify_cleavage_sites(protein_sequence: str, enzyme_regex: str):
    """
    Identifies cleavage sites in a protein sequence, based on enzyme_regex pattern.
    If enzyme_regex is a key in PROTEASES, the corresponding regex pattern will be used.
    Cleavage sites are identified as the positions after matching residues.

    Args:
        protein_sequence (str): The protein sequence to be processed.
        enzyme_regex (str): Regular expression pattern or key in PROTEASES dictionary that defines the cleavage
                            rule of the enzyme.

    Returns:
        List[int]: A list of positions (1-indexed) in the protein sequence where cleavage occurs.
    """

    if enzyme_regex in PROTEASES:
        enzyme_regex = PROTEASES[enzyme_regex]

    enzyme_sites = []
    for site in reg.finditer(enzyme_regex, protein_sequence, overlapped=True):
        enzyme_sites.append(site.span(0))
    return [site[0] + 1 for site in enzyme_sites]


def digest_sequence(sequence: str, enzyme_regex: str, missed_cleavages: int, min_len: int,
                    max_len: int, semi: bool) -> List[str]:
    """
    Digests a given amino acid sequence using specified enzyme rules and parameters, returning a list of peptides.

    The sequence is broken down into peptides based on the enzyme's cleavage rules, with the number of
    missed cleavages taken into account. The returned peptides are filtered based on minimum and maximum length.

    Args:
        sequence (str): The amino acid sequence to be digested.
        enzyme_regex (str): Regular expression representing the enzyme's cleavage rules.
        missed_cleavages (int): The number of allowed missed cleavages.
        min_len (int): Minimum length of the returned peptides.
        max_len (int): Maximum length of the returned peptides.
        semi (bool): Whether to include semi-digested peptides.

    Returns:
        List[str]: The list of digested peptides.
    """
    spans = build_spans(len(sequence), identify_cleavage_sites(sequence, enzyme_regex), missed_cleavages, min_len, max_len, semi)
    return [span_to_sequence(sequence, span) for span in spans]
