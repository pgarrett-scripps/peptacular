import re
from copy import deepcopy
from typing import Dict, List, Any, Generator

from peptacular.term import get_c_term_modification, strip_term_modifications, add_n_term_modification, \
    add_c_term_modification, get_n_term_modification


def calculate_sequence_length(sequence: str) -> int:
    """
    Compute the length of the peptide sequence, excluding terminal modification notations.

    :param sequence: The peptide sequence with potential terminal notations.
    :type sequence: str

    :return: The length of the peptide sequence.
    :rtype: int

    :Example:

    .. code-block:: python

        >>> calculate_sequence_length("[Acetyl]PEP(1.2345)TID(3.14)E[Amide]")
        7

        >>> calculate_sequence_length("PEP(1.2345)TID(3.14)E")
        7

        >>> calculate_sequence_length("[Acetyl]PEPTIDE[Amide]")
        7

        >>> calculate_sequence_length("PEPTIDE")
        7
    """

    return len(strip_modifications(sequence))


def parse_modifications(sequence: str) -> Dict[int, str]:
    """
    Parses a peptide sequence with modifications and returns a dictionary where keys
    represent the position of the modified amino acid and values are the respective modifications.

    :param sequence: The peptide sequence, potentially including modifications.
    :type sequence: str

    :return: A dictionary with the modification values indexed by the position
             of the modified amino acid.
    :rtype: Dict[int, str]

    :raises ValueError: If the sequence contains invalid modification notation.

    :Example:

    .. code-block:: python

        >>> parse_modifications('PEP(Phospho)TIDE')
        {2: 'Phospho'}

        >>> parse_modifications('[Acetyl]P(3.14)EPTIDE(1.234)[Amide]')
        {-1: 'Acetyl', 0: '3.14', 6: '1.234', 7: 'Amide'}

        >>> parse_modifications('PEPTIDE')
        {}
    """

    n_term_mod = get_n_term_modification(sequence)
    c_term_mod = get_c_term_modification(sequence)

    sequence = strip_term_modifications(sequence)
    stripped_sequence = strip_modifications(sequence)

    mod_pattern = re.compile(r'\(([^)]*)\)')
    matches = mod_pattern.finditer(sequence)

    modifications = {}

    if n_term_mod:
        modifications[-1] = n_term_mod

    mod_offset = 0

    for match in matches:
        modification_position = match.start() - mod_offset - 1
        modification_value = match.group(1)
        modifications[modification_position] = modification_value

        # Compute offset for the next modification
        mod_offset += len(match.group())

    if c_term_mod:
        modifications[len(stripped_sequence)] = c_term_mod

    return modifications


def _build_modified_sequence(unmodified_sequence: str, modifications: Dict[int, Any]) -> str:
    """
    Builds a modified peptide sequence from an unmodified sequence and a dictionary of modifications.
    """
    n_term_mod = modifications.pop(-1, None)
    c_term_mod = modifications.pop(calculate_sequence_length(unmodified_sequence), None)

    modified_sequence = []
    prev_index = 0

    # Sort the modifications by index in descending order
    for i, mod in sorted(modifications.items()):

        if i < -1 or i >= len(unmodified_sequence):
            raise ValueError(f'Index of modification: {i} is invalid for peptide sequence: {unmodified_sequence}')

        # Insert the modification into the modified peptide sequence
        modified_sequence.append(unmodified_sequence[prev_index: i + 1])
        modified_sequence.append(f"({mod})")
        prev_index = i + 1

    modified_sequence.append(unmodified_sequence[prev_index:])
    modified_sequence = ''.join(modified_sequence)
    modified_sequence = add_n_term_modification(modified_sequence, n_term_mod)
    modified_sequence = add_c_term_modification(modified_sequence, c_term_mod)
    return modified_sequence


def add_modifications(sequence: str, modifications: Dict[int, Any]) -> str:
    """
    Adds modifications to the given peptide sequence, overriding any existing modifications at overlapping indexes.

    The modifications are provided as a dictionary where keys correspond to the indices of the amino acids
    in the sequence to be modified, and values are the modifications to be applied. Modifications are inserted
    into the sequence in descending order of their indices to ensure their validity after each modification.

    :param sequence: Unmodified peptide sequence.
    :type sequence: str
    :param modifications: Dictionary with indices of the modified amino acids and corresponding modifications.
    :type modifications: Dict[int, str]

    :return: The peptide sequence with the specified modifications.
    :rtype: str

    :raises ValueError: If an index in modifications is invalid for the sequence.

    :Example:

    .. code-block:: python

        >>> add_modifications('PEPTIDE', {2: 'phospho'})
        'PEP(phospho)TIDE'

        >>> add_modifications('PEPTIDE', {-1: 'Acetyl', 6: '1.234', 7: 'Amide'})
        '[Acetyl]PEPTIDE(1.234)[Amide]'

        >>> add_modifications('PEPTIDE', {})
        'PEPTIDE'

        >>> add_modifications('[Acetyl]PEPTIDE(1.234)[Amide]', {-1: 'Amide', 7: 'Acetyl'})
        '[Amide]PEPTIDE(1.234)[Acetyl]'
    """

    stripped_sequence = strip_modifications(sequence)
    original_mods = parse_modifications(sequence)

    # Merge the original modifications with the new ones, overriding the original ones
    original_mods.update(modifications)

    return _build_modified_sequence(stripped_sequence, original_mods)


def strip_modifications(sequence: str) -> str:
    """
    Removes modifications from the given sequence.

    :param sequence: modified sequence to be stripped
    :type sequence: str

    :return: unmodified sequence
    :rtype: str

    :raises ValueError: If an index in modifications is invalid for the sequence.

    :Example:

    .. code-block:: python

        >>> strip_modifications('PEP(phospho)TIDE')
        'PEPTIDE'

        >>> strip_modifications('[Acetyl](3.14)PEPTIDE(1.234)[Amide]')
        'PEPTIDE'

        >>> strip_modifications('PEPTIDE')
        'PEPTIDE'

    """

    # Removing modifications in parentheses
    sequence = re.sub(r'\([^)]*\)', '', sequence)

    # Removing modifications in square brackets
    sequence = re.sub(r'\[[^\]]*\]', '', sequence)

    return sequence


def apply_static_modifications(sequence: str, mod_map: Dict[str, Any]) -> str:
    """
    Add static modifications to an amino acid sequence. If a modification is already present in the sequence,
    it will be replaced by the new one.

    :param sequence: Original amino acid sequence.
    :type sequence: str
    :param mod_map: Dictionary mapping amino acids to the mass of their modifications.
    :type mod_map: Dict[str, Any]

    :return: Modified amino acid sequence.
    :rtype: str

    :raises ValueError: If an index in modifications is invalid for the sequence.

    :Example:

    .. code-block:: python

        >>> apply_static_modifications('PEPTIDE', {'P': 'phospho'})
        'P(phospho)EP(phospho)TIDE'

        >>> apply_static_modifications('[3.14]PEPTIDE(1.234)', {'P': 'phospho'})
        '[3.14]P(phospho)EP(phospho)TIDE(1.234)'

        >>> apply_static_modifications('PEP(1.234)TIDE', {'P': 'phospho'})
        'P(phospho)EP(phospho)TIDE'

        >>> apply_static_modifications('PEPTIDE', {})
        'PEPTIDE'

        >>> apply_static_modifications('PEPTIDE', {'S': 'phospho'})
        'PEPTIDE'

    """

    stripped_sequence = strip_modifications(sequence)
    sequence_mod_map = parse_modifications(sequence)
    for aa, mass in mod_map.items():
        indexes = [i for i, x in enumerate(stripped_sequence) if x == aa]
        for index in indexes:
            sequence_mod_map[index] = mass

    return add_modifications(stripped_sequence, sequence_mod_map)


def _variable_modification_builder(mods: Dict[str, Any], sequence: str, curr_idx: int, mod_dict: Dict[int, Any],
                                   max_mod_count: int) -> Generator[Dict[int, str], None, None]:
    """
    Recursive function to generate all possible combinations of modifications on an amino acid sequence.

    :param mods: Dictionary mapping amino acids to the mass of their modifications.
    :type mods: Dict[str, Any]
    :param sequence: amino acid sequence to modify.
    :type sequence: str
    :param curr_idx: Current index on the amino acid sequence.
    :type curr_idx: int
    :param mod_dict: Dictionary mapping the indices of the amino acids to the mass of their modifications.
    :type mod_dict: Dict[int, str]
    :param max_mod_count: Maximum number of modifications allowed on the amino acid sequence.
    :type max_mod_count: int

    :yield: Modified dictionary mapping the indices of the amino acid sequence to the mass of their modifications.
    :rtype: Generator[Dict[int, float], None, None]
    """

    if curr_idx == len(sequence):
        yield mod_dict
        return

    original_mod_dict = deepcopy(mod_dict)
    for aa, mass in mods.items():

        if len(mod_dict) >= max_mod_count:
            break

        if sequence[curr_idx] == aa and curr_idx not in mod_dict:
            mod_dict[curr_idx] = mass
            yield from _variable_modification_builder(mods, sequence, curr_idx + 1, mod_dict, max_mod_count)

    yield from _variable_modification_builder(mods, sequence, curr_idx + 1, original_mod_dict, max_mod_count)


def apply_variable_modifications(sequence: str, mod_map: Dict[str, Any], max_mods: int) -> List[str]:
    """
    Add variable modifications to an amino acid sequence.

    :param sequence: Original amino acid sequence.
    :type sequence: str
    :param mod_map: Dictionary mapping amino acids to the mass of their modifications.
    :type mod_map: Dict[str, float]
    :param max_mods: Maximum number of modifications allowed on the peptide.
    :type max_mods: int

    :return: List of all possible modified peptide sequences, each represented as a string.
    :rtype: List[str]

    :raises ValueError: If an index in modifications is invalid for the sequence.

    :Example:

    .. code-block:: python

        >>> apply_variable_modifications('PEPTIDE', {'P': 'phospho'}, 2)
        ['P(phospho)EP(phospho)TIDE', 'P(phospho)EPTIDE', 'PEP(phospho)TIDE', 'PEPTIDE']

        >>> apply_variable_modifications('PEPTIDE', {'P': 'phospho'}, 0)
        ['PEPTIDE']

        >>> apply_variable_modifications('[Acetyl]P(3.14)EPTIDE[Amide]', {'P': 'phospho'}, 2)
        ['[Acetyl]P(3.14)EP(phospho)TIDE[Amide]', '[Acetyl]P(3.14)EPTIDE[Amide]']

    """

    original_mods = parse_modifications(sequence)
    stripped_sequence = strip_modifications(sequence)
    mods = _variable_modification_builder(mod_map, stripped_sequence, 0, original_mods, max_mods + len(original_mods))
    return [add_modifications(stripped_sequence, mod) for mod in mods]


def reverse_sequence(sequence: str, swap_terms: bool = False) -> str:
    """
    Reverses the sequence, while preserving the position of any modifications.

    :param sequence: The amino acid sequence to be reversed.
    :type sequence: str

    :param swap_terms: If True, the N- and C-terminal modifications will be swapped.
    :type swap_terms: bool

    :return: The reversed sequence with modifications preserved.
    :rtype: str

    :raises ValueError: If the sequence contains unmatched parentheses indicating incorrect modification notation.

    :Example:

    .. code-block:: python

        >>> reverse_sequence('PEPTIDE')
        'EDITPEP'

        >>> reverse_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]')
        '[Acetyl]EDITP(phospho)EP(phospho)[Amide]'

        >>> reverse_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', True)
        '[Amide]EDITP(phospho)EP(phospho)[Acetyl]'

    """

    mods = parse_modifications(sequence)
    stripped_sequence = strip_modifications(sequence)[::-1]

    n_term = mods.pop(-1, None)
    c_term = mods.pop(len(stripped_sequence), None)

    mod_reverse = {len(stripped_sequence) - 1 - k: v for k, v in mods.items()}

    modified_sequence = add_modifications(stripped_sequence, mod_reverse)

    if swap_terms:
        modified_sequence = add_n_term_modification(modified_sequence, c_term)
        modified_sequence = add_c_term_modification(modified_sequence, n_term)

    else:
        modified_sequence = add_n_term_modification(modified_sequence, n_term)
        modified_sequence = add_c_term_modification(modified_sequence, c_term)

    return modified_sequence


def shift_sequence_left(sequence: str, shift: int) -> str:
    """
    Shifts the sequence to the left by a given number of positions, while preserving the position of any modifications.

    :param sequence: The sequence to be shifted.
    :type sequence: str
    :param shift: The number of positions to shift the sequence to the left.
    :type shift: int

    :return: The shifted sequence with modifications preserved.
    :rtype: str

    :raises ValueError: If the sequence contains unmatched parentheses indicating incorrect modification notation.

    :Example:

    .. code-block:: python

        >>> shift_sequence_left('PEPTIDE', 2)
        'PTIDEPE'

        >>> shift_sequence_left('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', 2)
        '[Acetyl]P(phospho)TIDEP(phospho)E[Amide]'

        >>> shift_sequence_left('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', 7)
        '[Acetyl]P(phospho)EP(phospho)TIDE[Amide]'

        >>> shift_sequence_left('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', -2)
        '[Acetyl]DEP(phospho)EP(phospho)TI[Amide]'

        >>> shift_sequence_left('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', 0)
        '[Acetyl]P(phospho)EP(phospho)TIDE[Amide]'

    """

    mods = parse_modifications(sequence)
    stripped_sequence = strip_modifications(sequence)

    n_term = mods.pop(-1, None)
    c_term = mods.pop(len(stripped_sequence), None)

    seq_len = len(stripped_sequence)

    # Take modulus to ensure proper wrapping around
    effective_shift = shift % seq_len

    shifted_stripped_sequence = stripped_sequence[effective_shift:] + stripped_sequence[:effective_shift]

    # Shift mod positions, taking care to avoid negative keys
    shifted_sequence_mods = {(k - effective_shift) % seq_len if k != -1 else -1: v for k, v in mods.items()}

    modified_sequence = add_modifications(shifted_stripped_sequence, shifted_sequence_mods)
    modified_sequence = add_n_term_modification(modified_sequence, n_term)
    modified_sequence = add_c_term_modification(modified_sequence, c_term)

    return modified_sequence


def is_sequence_valid(sequence: str) -> bool:
    """
    Checks if the sequence is valid. A sequence is valid if it does not contain doubly enclosed modifications.
    :param sequence: modified sequence
    :type sequence: str
    :return: True, if the sequence is valid, False otherwise
    :rtype: bool

    :Example:

    .. code-block:: python

        >>> is_sequence_valid('PEPTIDE')
        True

        >>> is_sequence_valid('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]')
        True

        >>> is_sequence_valid('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]]')
        False

        >>> is_sequence_valid('[Acetyl]P(phospho)EP((phospho))TIDE[Amide]')
        False
    """

    stripped_sequence = strip_modifications(sequence)
    residues = set(stripped_sequence)

    # check if parenthesis are in the stripped sequence
    if '(' in residues or ')' in residues:
        return False

    # check if square brackets are in the stripped sequence
    if '[' in residues or ']' in residues:
        return False

    # check modifications
    mods = parse_modifications(sequence)
    for index, mod in mods.items():
        if '(' in mod or ')' in mod:
            return False
        if '[' in mod or ']' in mod:
            return False

    return True
