"""
Sequence.py - Amino Acid Sequence Manipulation with Modifications

This module provides a set of utility functions to work with amino acid sequences that may contain modifications.
Such modifications can occur at any position within the sequence and can also be present at the N-terminus or
C-terminus of a peptide sequence.

Key Features:
- Parsing sequences with modifications and extracting modification details.
- Calculating the length of sequences excluding any modification notations.
- Adding, shifting, reversing, and stripping modifications from sequences.
- Applying both static and variable modifications to sequences.
- Validating the structure and integrity of modified sequences.

Modification Notations:
- Term modifications (N-terminus and C-terminus) are specified using square brackets: `[]`.
- Residue modifications within the sequence are specified using parentheses: `()`.

Modification Types and Representation:
- Modifications can be represented as strings (e.g., "Phospho", "Acetyl"), integers, or floats.
- String-based modifications typically represent the name or type of the modification.
- Numerical modifications (integers or floats) represent specific mass changes associated with a modification.
- During parsing, the module automatically identifies the modification type based on its representation.

Example Sequences:
    [Acetyl]P(Oxydation)EPTI(Phospho)DE[Amide] or [1.234]PE(1.0)PTI(2.0)DE[3.1415]

The module ensures that sequences and their modifications are handled correctly, preserving the positional
relationship between amino acids and their modifications during operations.
"""

import re
from copy import deepcopy
from typing import Dict, List, Any, Generator, Union, Set, Tuple

from peptacular.term.modification import get_c_term_modification, strip_term_modifications, add_n_term_modification, \
    add_c_term_modification, get_n_term_modification
from peptacular.term.residue import pop_n_term_residue
from peptacular.util import convert_type, identify_regex_indexes, _validate_span


def calculate_sequence_length(sequence: str) -> int:
    """
    Compute the length of the peptide sequence, excluding any modification notation.

    :param sequence: The amino acid sequence
    :type sequence: str
    :return: Length of the unmodified sequence.
    :rtype: int

    .. code-block:: python

        >>> calculate_sequence_length("[Acetyl]PEP(1.2345)TID(3.14)E[Amide]")
        7
        >>> calculate_sequence_length("PEPTIDE")
        7

    """

    return len(strip_modifications(sequence))


def get_modifications(sequence: str) -> Dict[int, Union[str, float, int]]:
    """
    Parses a peptide sequence with modifications and returns a dictionary where keys
    represent the position of the modified amino acid and values are the respective modifications.

    :param sequence: The peptide sequence, potentially including modifications.
    :type sequence: str
    :return: A dictionary with the modification values indexed by the position
             of the modified amino acid.
    :rtype: Dict[int, Union[str, int, float]]

    .. code-block:: python

        # All modifications will be returned as either strings, ints or floats.
        >>> get_modifications('PEP(Phospho)T(1)IDE(3.14)')
        {2: 'Phospho', 3: 1, 6: 3.14}

        # N-terminus and C-terminal modifications are index by -1, and the length of the sequence, respectively:
        >>> get_modifications('[Acetyl]PEPTIDE[Amide]')
        {-1: 'Acetyl', 7: 'Amide'}

        # A sequecnes without any modifications will return an empty dictionary:
        >>> get_modifications('PEPTIDE')
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
        mod_index = match.start() - mod_offset - 1
        mod = convert_type(match.group(1))
        modifications[mod_index] = mod

        # Compute offset for the next modification
        mod_offset += len(match.group())

    if c_term_mod:
        modifications[len(stripped_sequence)] = c_term_mod

    return modifications


def _construct_sequence_with_modifications(sequence: str, mod_map: Dict[int, Any]) -> str:
    """
    Builds a modified peptide sequence from an unmodified sequence and a dictionary of modifications.
    """

    n_term_mod = mod_map.pop(-1, None)
    c_term_mod = mod_map.pop(calculate_sequence_length(sequence), None)

    sequence_components = []

    if n_term_mod:
        sequence_components.append(f"[{n_term_mod}]")

    prev_index = 0

    # Sort the modifications by index in descending order
    for mod_index, mod in sorted(mod_map.items()):

        if mod_index < -1 or mod_index >= len(sequence):
            raise ValueError(f'Index of modification: {mod_index} is invalid for peptide sequence: {sequence}')

        # Insert the modification into the modified peptide sequence
        sequence_components.append(sequence[prev_index: mod_index + 1])
        sequence_components.append(f"({mod})")
        prev_index = mod_index + 1

    sequence_components.append(sequence[prev_index:])

    if c_term_mod:
        sequence_components.append(f"[{c_term_mod}]")

    return ''.join(sequence_components)


def add_modifications(sequence: str, modifications: Dict[int, Any], overwrite: bool = True) -> str:
    """
    Adds modifications to the given peptide sequence.

    :param sequence: Unmodified peptide sequence.
    :type sequence: str
    :param modifications: Dictionary with indices of the modified amino acids and corresponding modifications.
    :type modifications: Dict[int, Any]
    :param overwrite: If True, existing modifications will be overwritten. If False, existing modifications will be
                        preserved.
    :type overwrite: bool

    :return: The peptide sequence with the specified modifications.
    :rtype: str

    .. code-block:: python

        # Add internal modifications to an unmodified peptide
        >>> add_modifications('PEPTIDE', {2: 'phospho'})
        'PEP(phospho)TIDE'

        # Can also add N and C terminal modifications
        >>> add_modifications('PEPTIDE', {-1: 'Acetyl', 6: '1.234', 7: 'Amide'})
        '[Acetyl]PEPTIDE(1.234)[Amide]'

        # Empty modification dicts will simply return the input sequence
        >>> add_modifications('PEPTIDE', {})
        'PEPTIDE'

        # Existing modifications will be overwritten by default
        >>> add_modifications('PEP(phospho)TIDE', {2: 'acetyl'})
        'PEP(acetyl)TIDE'

        # Can also preserve existing modifications
        >>> add_modifications('PEP(phospho)TIDE', {2: 'acetyl'}, overwrite=False)
        'PEP(phospho)TIDE'

    """

    stripped_sequence = strip_modifications(sequence)
    original_mods = get_modifications(sequence)

    for mod_index, mod in sorted(modifications.items()):
        if mod is None:
            continue
        if mod_index in original_mods:
            if overwrite is True:
                original_mods[mod_index] = mod
        else:
            original_mods[mod_index] = mod

    return _construct_sequence_with_modifications(stripped_sequence, original_mods)


def pop_modifications(sequence: str) -> Tuple[str, Dict[int, Any]]:
    """
    Removes all modifications from the given sequence, returning the unmodified sequence and a dictionary of the
    removed modifications.

    :param sequence: The sequence to be stripped of modifications.
    :type sequence: str

    :return: A tuple containing the unmodified sequence and a dictionary of the removed modifications.
    :rtype: Tuple[str, Dict[int, Any]]

    .. code-block:: python

        # Simply combines the functionality of strip_modifications and get_modifications
        >>> pop_modifications('PEP(phospho)TIDE')
        ('PEPTIDE', {2: 'phospho'})

    """

    return strip_modifications(sequence), get_modifications(sequence)


def strip_modifications(sequence: str) -> str:
    """
    Strips all modifications from the given sequence, returning the unmodified sequence.

    :param sequence: The sequence to be stripped of modifications.
    :type sequence: str

    :return: The stripped sequence
    :rtype: str

    .. code-block:: python

        # Removes internal modifications:
        >>> strip_modifications('PEP(phospho)TIDE')
        'PEPTIDE'

        # Also removes N and C terminal modifications:
        >>> strip_modifications('[Acetyl](3.14)PEPTIDE(1.234)[Amide]')
        'PEPTIDE'

        # Using a sequence without modifications will return the same sequence:
        >>> strip_modifications('PEPTIDE')
        'PEPTIDE'

    """

    # Removing modifications in parentheses
    sequence = re.sub(r'\([^)]*\)', '', sequence)

    # Removing modifications in square brackets
    sequence = re.sub(r'\[[^]]*]', '', sequence)

    return sequence


def apply_static_modifications(sequence: str, mod_map: Dict[str, Any], overwrite: bool = True) -> str:
    """
    Add static modifications to an amino acid sequence. If a modification is already present in the sequence,
    it will be replaced by the new one.

    :param sequence: Original amino acid sequence.
    :type sequence: str
    :param mod_map: Dictionary mapping amino acids to the mass of their modifications.
    :type mod_map: Dict[str, Any]
    :param overwrite: If True, existing modifications will be replaced by the new ones.
    :type overwrite: bool

    :return: Modified amino acid sequence.
    :rtype: str

    .. code-block:: python

        # Applies static modifcation to all matching residues in the sequence:
        >>> apply_static_modifications('PEPTIDE', {'P': 'phospho'})
        'P(phospho)EP(phospho)TIDE'

        # Can specify multiple static modifications:
        >>> apply_static_modifications('PEPTIDE', {'P': 'phospho', 'E': '3.0'})
        'P(phospho)E(3.0)P(phospho)TIDE(3.0)'

        # Works with already modified sequences:
        >>> apply_static_modifications('[3.14]PEPTIDE(1.234)', {'P': 'phospho'})
        '[3.14]P(phospho)EP(phospho)TIDE(1.234)'

        # By default, any overlapping modifications will be replaced:
        >>> apply_static_modifications('PEP(1.234)TIDE', {'P': 'phospho'})
        'P(phospho)EP(phospho)TIDE'

        # To not reaplce existing modifications, set override to False:
        >>> apply_static_modifications('PEP(1.234)TIDE', {'P': 'phospho'}, overwrite=False)
        'P(phospho)EP(1.234)TIDE'

        # Can also use regular expressions to match multiple residues:
        >>> apply_static_modifications('PEPTIDE', {'(?<=P)E': 'phospho'})
        'PE(phospho)PTIDE'

        >>> apply_static_modifications('PEPTIDE', {'PE': 'phospho'})
        'P(phospho)EPTIDE'

    """

    stripped_sequence = strip_modifications(sequence)
    original_mod_map = get_modifications(sequence)

    for regex_str, mod_mass in mod_map.items():
        for mod_index in identify_regex_indexes(stripped_sequence, regex_str):
            if overwrite is True or mod_index not in original_mod_map:
                original_mod_map[mod_index] = mod_mass

    return add_modifications(stripped_sequence, original_mod_map)


def _apply_variable_modifications_rec(mods: Dict[int, Set], sequence: str, index: int, current_mods: Dict[int, Any],
                                      max_mod_count: int) -> Generator[Dict[int, Any], None, None]:
    """
    Recursively applies variable modifications to the amino acid sequence.

    :param mods: Dictionary mapping amino acids to the mass of their modifications.
    :type mods: Dict[int, Any]
    :param sequence: the unmodified sequence.
    :type sequence: str
    :param index: Current index on the amino acid sequence.
    :type index: int
    :param current_mods: Dictionary mapping the indices of the amino acids to the mass of their modifications.
    :type current_mods: Dict[int, str]
    :param max_mod_count: Maximum number of modifications allowed on the amino acid sequence.
    :type max_mod_count: int

    :yield: Modified dictionary mapping the indices of the amino acid sequence to the mass of their modifications.
    :rtype: Generator[Dict[int, float], None, None]
    """

    if index == len(sequence) or len(current_mods) == max_mod_count:
        yield current_mods
        return

    original_mod_dict = deepcopy(current_mods)

    curr_mods = mods.get(index, None)

    if curr_mods is not None:
        for curr_mod in curr_mods:
            updated_mod_dict = deepcopy(current_mods)
            if index not in current_mods:
                updated_mod_dict[index] = curr_mod
                yield from _apply_variable_modifications_rec(mods, sequence, index + 1, updated_mod_dict,
                                                             max_mod_count)
    yield from _apply_variable_modifications_rec(mods, sequence, index + 1, original_mod_dict, max_mod_count)


def apply_variable_modifications(sequence: str, mod_map: Dict[str, Any], max_mods: int) -> List[str]:
    """
    Apply variable modifications to a sequence.

    :param sequence: Original amino acid sequence.
    :type sequence: str
    :param mod_map: Dictionary mapping amino acids to the mass of their modifications.
    :type mod_map: Dict[str, Any]
    :param max_mods: Maximum number of modifications allowed on the peptide.
    :type max_mods: int

    :return: List of all possible modified peptide sequences.
    :rtype: List[str]

    .. code-block:: python

        >>> apply_variable_modifications('PEPTIDE', {'P': 'phospho'}, 2)
        ['P(phospho)EP(phospho)TIDE', 'P(phospho)EPTIDE', 'PEP(phospho)TIDE', 'PEPTIDE']

        # Setting max_mods to 0 will return the original sequence:
        >>> apply_variable_modifications('PEPTIDE', {'P': 'phospho'}, 0)
        ['PEPTIDE']

        # Works with already modified sequences, and will not replace existing modifications:
        >>> apply_variable_modifications('[Acetyl]P(3.14)EPTIDE[Amide]', {'P': 'phospho'}, 2)
        ['[Acetyl]P(3.14)EP(phospho)TIDE[Amide]', '[Acetyl]P(3.14)EPTIDE[Amide]']

        # Can also use regular expressions to match multiple residues:
        >>> apply_variable_modifications('PEPTIDE', {'(?<=P)E': 'phospho'}, 2)
        ['PE(phospho)PTIDE', 'PEPTIDE']

    """

    original_mods = get_modifications(sequence)
    stripped_sequence = strip_modifications(sequence)

    new_mod_map: Dict[int, Set] = {}
    for regex_str, mod_mass in mod_map.items():
        for mod_index in identify_regex_indexes(stripped_sequence, regex_str):
            new_mod_map.setdefault(mod_index, set()).add(mod_mass)

    mods = _apply_variable_modifications_rec(new_mod_map, stripped_sequence, 0, original_mods,
                                             max_mods + len(original_mods))
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

    .. code-block:: python

        # For unmodified sequences, the result is the same as the built-in string reversal:
        >>> reverse_sequence('PEPTIDE')
        'EDITPEP'

        # For modified sequences, the modifications are preserved on the associated residues:
        >>> reverse_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]')
        '[Acetyl]EDITP(phospho)EP(phospho)[Amide]'

        # If swap_terms is True, the N- and C-terminal modifications will be swapped too:
        >>> reverse_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', swap_terms=True)
        '[Amide]EDITP(phospho)EP(phospho)[Acetyl]'

    """

    mods = get_modifications(sequence)
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


def shift_sequence(sequence: str, shift: int) -> str:
    """
    Shifts the sequence to the left by a given number of positions, while preserving the position of any modifications.

    :param sequence: The sequence to be shifted.
    :type sequence: str
    :param shift: The number of positions to shift the sequence to the left.
    :type shift: int

    :return: The shifted sequence with modifications preserved.
    :rtype: str

    .. code-block:: python

        >>> shift_sequence('PEPTIDE', 2)
        'PTIDEPE'

        >>> shift_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', 2)
        '[Acetyl]P(phospho)TIDEP(phospho)E[Amide]'

        # Shifting by 0 or length of sequence positions returns the original sequence:
        >>> shift_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', 7)
        '[Acetyl]P(phospho)EP(phospho)TIDE[Amide]'
        >>> shift_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', 0)
        '[Acetyl]P(phospho)EP(phospho)TIDE[Amide]'

        # Shifting by a negative number shifts the sequence to the right:
        >>> shift_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]', -2)
        '[Acetyl]DEP(phospho)EP(phospho)TI[Amide]'

    """

    mods = get_modifications(sequence)
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
    mods = get_modifications(sequence)
    for _, mod in mods.items():
        if isinstance(mod, str):
            if '(' in mod or ')' in mod:
                return False
            if '[' in mod or ']' in mod:
                return False

    return True


def span_to_sequence(sequence: str, span: Tuple[int, int, int]) -> str:
    """
    Extracts a subsequence from the input sequence based on the provided span.

    :param sequence: The original sequence.
    :type sequence: str
    :param span: A tuple representing the span of the subsequence to be extracted.
    :type span: Tuple[int, int, int]

    :return: The subsequence of the input sequence defined by the span.
    :rtype: str

    .. code-block:: python

        # Works with unmodified sequences
        >>> span_to_sequence('PEPTIDE', (0, 4, 0))
        'PEPT'

        >>> span_to_sequence('PEPTIDE', (1, 6, 0))
        'EPTID'

        >>> span_to_sequence('PEPTIDE', (4, 7, 0))
        'IDE'

        # but will also preserve modifications, including terminal modifications
        >>> span_to_sequence('[Acetyl]P(1.0)EPTIDE', (0, 4, 0))
        '[Acetyl]P(1.0)EPT'

        >>> span_to_sequence('PEPTIDE(1.0)[Amide]', (4, 7, 0))
        'IDE(1.0)[Amide]'

    """

    _validate_span(span)

    stripped_sequence = strip_modifications(sequence)

    base_sequence = stripped_sequence[span[0]:span[1]]

    # if peptide is unmodified skip the modification parsing
    if stripped_sequence == sequence:
        return base_sequence

    mods = get_modifications(sequence)
    n_term = get_n_term_modification(sequence)
    c_term = get_c_term_modification(sequence)

    mods = {k - span[0]: mods[k] for k in mods if span[0] <= k < span[1]}

    if span[0] == 0 and n_term is not None:
        mods[-1] = n_term

    if span[1] == len(stripped_sequence) and c_term is not None:
        mods[len(base_sequence)] = c_term

    return add_modifications(base_sequence, mods)


def split_sequence(sequence: str) -> List[str]:
    """
    Splits sequence into a list of amino acids, preserving modifications.

    :param sequence: The sequence to be split.
    :type sequence: str

    :return: A list of amino acids, preserving modifications.
    :rtype: List[str]

    .. code-block:: python

        >>> split_sequence('PEPTIDE')
        ['P', 'E', 'P', 'T', 'I', 'D', 'E']

        >>> split_sequence('[Acetyl]P(phospho)EP(phospho)TIDE[Amide]')
        ['[Acetyl]P(phospho)', 'E', 'P(phospho)', 'T', 'I', 'D', 'E[Amide]']

    """

    split_sequences = []
    while sequence:
        n_term, sequence = pop_n_term_residue(sequence)
        split_sequences.append(n_term)

    return split_sequences
