import re
from typing import Dict

from .util import check_parentheses


def convert_ip2_mod_to_uniprot_mod(peptide_sequence: str, mod_dict: Dict[str, str]):
    """
    Convert a peptide sequence with IP2-style modifications to UniProt-style modifications.

    This function takes a peptide sequence in string format and a dictionary of IP2-style
    modifications and their corresponding UniProt-style modifications, and returns a
    modified peptide sequence with the IP2-style modifications replaced with their
    corresponding UniProt-style modifications.
    """

    for ip2_mod, uniprot_mod in mod_dict.items():
        ip2_mod = f'({ip2_mod})'
        uniprot_mod = f'({uniprot_mod})'
        peptide_sequence = peptide_sequence.replace(ip2_mod, uniprot_mod)

    return peptide_sequence


def convert_to_ms2_pip_style(peptide_sequence: str) -> str:
    """
    Convert a peptide sequence to MS2PIP-style format.

    This function takes a peptide sequence in string format and converts it to a modified
    format that is used by the MS2PIP prediction software. The modified format consists of
    a series of location and modification pairs, separated by a pipe symbol ('|').

    Example:
    >>> convert_to_ms2_pip_style('PEPTIDEK')
    '2|Phospho|6|Lysine'
    """

    mod_dict = parse_modified_peptide(peptide_sequence)
    locs, mods = [], []
    for loc in mod_dict:
        mod = mod_dict[loc]
        ms2_pip_style_loc = loc + 1
        if ms2_pip_style_loc == len(strip_modifications(peptide_sequence)):
            ms2_pip_style_loc = -1
        locs.append(ms2_pip_style_loc)
        mods.append(mod)

    return '|'.join([f'{loc}|{mod}' for loc, mod in zip(locs, mods)])


def parse_modified_peptide(peptide_sequence: str) -> Dict[int, str]:
    """
    This function reads a peptide sequence with modifications and returns a dictionary
    with the modification values indexed by the position of the modified amino acid.
    """

    if check_parentheses(peptide_sequence) is False:
        raise ValueError(f'Incorrect modification notation in peptide sequence : {peptide_sequence}!')

    matches = re.finditer(r'\(([^)]*)\)', peptide_sequence)

    modifications = {}
    mod_offset = 0

    # Loop through the substrings
    for match in matches:
        modifications[match.start() - mod_offset - 1] = match.group()[1:-1]
        mod_offset += match.end() - match.start()
    return modifications


def create_modified_peptide(unmodified_sequence: str, modifications: Dict[int, str]) -> str:
    """
    Creates a modified peptide sequence from an unmodified peptide sequence and a dictionary of modifications.

    The modifications are specified as a dictionary where the keys are the indices of the modified amino acids
    in the peptide sequence, and the values are the modifications to apply at those indices. The modifications are
    added to the peptide sequence in descending order of index, so that the indices remain valid after each
    modification is applied.
    """

    modified_sequence = []
    prev_index = 0

    # Sort the modifications by index in descending order
    for i, mod in sorted(modifications.items()):

        if i < 0 or i >= len(unmodified_sequence):
            raise ValueError(f'Index of modification: {i} is invalid for peptide sequence: {unmodified_sequence}')

        # Insert the modification into the modified peptide sequence
        modified_sequence.append(unmodified_sequence[prev_index: i + 1])
        modified_sequence.append(f"({mod})")
        prev_index = i + 1

    modified_sequence.append(unmodified_sequence[prev_index:])
    return ''.join(modified_sequence)


def strip_modifications(peptide_sequence: str) -> str:
    """
    Removes any non-amino-acid characters from the given peptide sequence.
    """

    if check_parentheses(peptide_sequence) is False:
        raise ValueError(f'Incorrect modification notation in peptide sequence : {peptide_sequence}!')

    unmodified_sequence = re.sub(r'\([^)]*\)', '', peptide_sequence)
    return unmodified_sequence


def get_left_sequences(sequence: str, min_len: int = None, max_len: int = None):
    """
    Returns a set of left substrings of string `sequence` that have lengths between `min_len` and `max_len`.

    Example:
    >>> get_left_sequences("abc", 1, 2)
    {'a', 'ab'}
    """

    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = len(sequence)
    return {sequence[0:i] for i in range(min_len, min(max_len, len(sequence)) + 1)}


def get_right_sequences(sequence: str, min_len: int = None, max_len: int = None):
    """
    Returns a set of right substrings of string `sequence` that have lengths between `min_len` and `max_len`.

    Example:
    >>> get_right_sequences("abc", 1, 2)
    {'c', 'bc'}
    """

    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = len(sequence)
    return {sequence[-i:] for i in range(min_len, min(max_len, len(sequence)) + 1)}


def get_semi_sequences(sequence: str, min_len: int = None, max_len: int = None):
    """
    Returns a set of all semi-peptides of string `sequence` that have lengths between `min_len` and `max_len`.

    Example:
    >>> get_semi_sequences("abc", 1, 2)
    {'a', 'ab', 'c', 'bc'}
    """

    return get_left_sequences(sequence, min_len, max_len).union(get_right_sequences(sequence, min_len, max_len))


def get_non_enzymatic_sequences(sequence: str, min_len: int = None, max_len: int = None):
    """
    Returns a set of all non-enzymatic peptides of string `sequence` that have lengths between `min_len` and `max_len`.

    Example:
    >>> get_non_enzymatic_sequences("abc", 1, 2)
    {'a', 'b', 'c', 'ab', 'bc'}
    """

    if min_len is None:
        min_len = 1

    if max_len is None:
        max_len = len(sequence)

    substrings = set()
    for i in range(len(sequence)):
        for j in range(i + min_len, min(i + max_len + 1, len(sequence) + 1)):
            substrings.add(sequence[i:j])
    return substrings
