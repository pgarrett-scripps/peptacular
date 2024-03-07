"""
Sequence.py - Amino Acid Sequence Manipulation with Modifications

Modifications can be a str, int, or float

Modifications are mapped to the unmodified sequence by their index.

Special modifications are mapped to the sequence by the following:
    - N-terminal modifications are mapped to the index 'n'
    - C-terminal modifications are mapped to the index 'c'
    - Isotopic modifications are mapped to the index 'i'
    - Static modifications are mapped to the index 's'
    - Labile modifications are mapped to the index 'l'
    - Unknown modifications are mapped to the index 'u'


Global Mods and labile mods are popped from the sequence and returned as a tuple with the stripped sequence. As
a result they can be positioned anywhere in the sequence. Maybe add a check to see if they are at the beginning or
end of the sequence and if not raise an error.

"""

import collections
import random
import regex as re
from typing import Dict, List, Union, Tuple, Counter, Callable

from peptacular.sequence.global_mods import pop_global_mods, parse_static_mods
from peptacular.sequence.labile_mods import pop_labile_mods
from peptacular.spans import Span
from peptacular.util import convert_type, map_bracket_content_to_index

from peptacular.types import ModDict


def sequence_length(sequence: str) -> int:
    """
    Compute the length of the peptide sequence, excluding any modification notation.

    :param sequence: The amino acid sequence
    :type sequence: str
    :return: Length of the unmodified sequence.
    :rtype: int

    .. code-block:: python

        >>> sequence_length("[Acetyl]-PEP[1.2345]TID[3.14]E-[Amide]")
        7
        >>> sequence_length("PEPTIDE")
        7
        >>> sequence_length("<C13>PEPTIDE[1.2345]")
        7

        # Skips ambiguous sequence notation: (?)
        >>> sequence_length("(?PE)PTIDE[1.2345]")
        7

    """

    stripped_sequence = strip_mods(sequence)
    cnt = 0
    for c in stripped_sequence:
        if c.isalpha():
            cnt += 1

    return cnt


def is_ambiguous(sequence: str) -> bool:
    """
    Check if the sequence contains ambiguous amino acids.

    :param sequence: The amino acid sequence
    :type sequence: str

    :return: True if the sequence contains ambiguous amino acids, False otherwise.
    :rtype: bool

    .. code-block:: python

        >>> is_ambiguous("PEPTIDE")
        False
        >>> is_ambiguous("PEPTIDE")
        False
        >>> is_ambiguous("PEPTIDE[1.2345]")
        False
        >>> is_ambiguous("<13C>PEPTIDE")
        False
        >>> is_ambiguous("(?PE)PEPTIDE")
        True

    """

    stripped_sequence = strip_mods(sequence)

    if '?' in stripped_sequence:
        return True

    if '(' in stripped_sequence:
        return True

    return False


def get_mods(sequence: str) -> ModDict:
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
        >>> get_mods('PEP[Phospho]T[1]IDE[3.14]')
        {2: ['Phospho'], 3: [1], 6: [3.14]}

        >>> get_mods('PEP[Phospho][1.0]TIDE')
        {2: ['Phospho', 1.0]}

        # N-terminus and C-terminal modifications are index by -1, and the length of the sequence, respectively:
        >>> get_mods('[Acetyl]-PEPTIDE-[Amide]')
        {'n': ['Acetyl'], 'c': ['Amide']}

        >>> get_mods('[Acetyl]-PEPTIDE')
        {'n': ['Acetyl']}

        >>> get_mods('PEPTIDE-[Amide]')
        {'c': ['Amide']}

        # A sequecnes without any modifications will return an empty dictionary:
        >>> get_mods('PEPTIDE')
        {}

        # gets labile modifications
        >>> get_mods('{1.0}PEPTIDE')
        {'l': [1.0]}

        # gets isotopic and static modifications
        >>> get_mods('<13C><15N><[Acetyl]@C>PEPTIDE')
        {'i': ['13C', '15N'], 's': ['[Acetyl]@C']}

        # get multipe sequence
        >>> get_mods('PEP[1][2]TIDE')
        {2: [1, 2]}

        # get multipe sequence
        >>> get_mods('PEP[1]^2TIDE')
        {2: [1, 1]}

        >>> get_mods('[1]?PEP[1]^2TIDE')
        {'u': [1], 2: [1, 1]}

        >>> get_mods('[1][2]^2?[100]^3-PEP[1]^2TIDE')
        {'u': [1, 2, 2], 'n': [100, 100, 100], 2: [1, 1]}

        >>> get_mods('(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K')
        {10: ['Oxidation'], 20: ['Oxidation']}

        >>> get_mods('P[2]EP[Formula:[13C]H12]TIDE')
        {0: [2], 2: ['Formula:[13C]H12']}

        >>> get_mods('P[2]EPTIDE')
        {0: [2]}

        >>> get_mods('[P]PEP[Formula:[13C]H12]TIDE')
        Traceback (most recent call last):
        ValueError: Invalid sequence

        >>> get_mods('?[P]PEP[Formula:[13C]H12]TIDE')
        Traceback (most recent call last):
        ValueError: Invalid sequence

        >>> get_mods('??[P]??-[P]PEP[Formula:[13C]H12]TIDE')
        Traceback (most recent call last):
        ValueError: Invalid sequence

    """

    if '<' not in sequence and '[' not in sequence and '{' not in sequence:
        return {}

    # get labile modifications
    sequence, labile_mods = pop_labile_mods(sequence)
    sequence, isotope_mods, static_mods = pop_global_mods(sequence)
    sequence, bracket_mods = map_bracket_content_to_index(sequence, False)

    seq_len = len(sequence)
    mod_dict = {}
    term_offset = 0
    for i, bracket_mods in bracket_mods.items():

        bracket_mods = [convert_type(mod) for mod in bracket_mods]

        if sequence[i] == '-':  # n or c term
            if i == seq_len - 1:
                if 'c' in mod_dict:
                    raise ValueError('Invalid sequence: Multiple locations for C-terminal modification')
                mod_dict['c'] = bracket_mods
            else:
                if 'n' in mod_dict:
                    raise ValueError('Invalid sequence: Multiple locations for N-terminal modification')
                mod_dict['n'] = bracket_mods
                term_offset += 1

        elif sequence[i] == '?':
            if 'u' in mod_dict:
                raise ValueError('Invalid sequence: Multiple locations for unknown modification')
            mod_dict['u'] = bracket_mods
            term_offset += 1
        else:
            index = i - term_offset

            mod_dict[index] = bracket_mods

    if labile_mods:
        mod_dict['l'] = labile_mods

    if isotope_mods:
        mod_dict['i'] = isotope_mods

    if static_mods:
        mod_dict['s'] = static_mods

    return mod_dict


def _construct_sequence(sequence: str, mods: ModDict) -> str:
    """
    Builds a modified peptide sequence from an unmodified sequence and a dictionary of modifications.
    """

    n_term_mods = mods.get('n', [])
    c_term_mods = mods.get('c', [])
    isotope_mods = mods.get('i', [])
    static_mods = mods.get('s', [])
    labile_mods = mods.get('l', [])
    unknown_mods = mods.get('u', [])

    sequence_components = []

    prev_index = 0

    # Sort the modifications by index in descending order
    for mod_index in sorted([k for k in mods.keys() if isinstance(k, int)]):

        if mod_index < 0 or mod_index > len(sequence):
            raise ValueError(f'Index of modification: {mod_index} is invalid for peptide sequence: {sequence}')

        mod = ''.join([f'[{m}]' for m in mods[mod_index]])

        # Insert the modification into the modified peptide sequence
        sequence_components.append(sequence[prev_index: mod_index + 1])
        sequence_components.append(f"{mod}")
        prev_index = mod_index + 1

    sequence_components.append(sequence[prev_index:])

    if n_term_mods:
        n_term_mod = ''.join([f'[{m}]' for m in n_term_mods]) + '-'
        sequence_components.insert(0, n_term_mod)

    if c_term_mods:
        c_term_mod = '-' + ''.join([f'[{m}]' for m in c_term_mods])
        sequence_components.append(c_term_mod)

    if unknown_mods:
        unknown_mod = ''.join([f'[{m}]' for m in unknown_mods]) + '?'
        sequence_components.insert(0, unknown_mod)

    if isotope_mods:
        isotope_mod = ''.join([f'<{m}>' for m in isotope_mods])
        sequence_components.insert(0, isotope_mod)

    if static_mods:
        static_mod = ''.join([f'<{m}>' for m in static_mods])
        sequence_components.insert(0, static_mod)

    if labile_mods:
        labile_mod = ''.join(['{' + str(m) + '}' for m in labile_mods])
        sequence_components.insert(0, labile_mod)

    return ''.join(sequence_components)


def add_mods(sequence: str, mods: ModDict, overwrite: bool = False) -> str:
    """
    Adds modifications to the given peptide sequence.

    :param sequence: Unmodified peptide sequence.
    :type sequence: str
    :param mods: Dictionary with indices of the modified amino acids and corresponding modifications.
    :type mods: Dict[int, Any]
    :param overwrite: If True, existing modifications will be overwritten. If False, existing modifications will be
                        preserved.
    :type overwrite: bool

    :return: The peptide sequence with the specified modifications.
    :rtype: str

    .. code-block:: python

        # Add internal modifications to an unmodified peptide
        >>> add_mods('PEPTIDE', {2: ['phospho']})
        'PEP[phospho]TIDE'

        # Can also add N and C terminal modifications
        >>> add_mods('PEPTIDE', {'n': ['Acetyl'], 6: ['1.234'], 'c': ['Amide']})
        '[Acetyl]-PEPTIDE[1.234]-[Amide]'

        # Empty modification dicts will simply return the input sequence
        >>> add_mods('PEPTIDE', {})
        'PEPTIDE'

        >>> add_mods('PEP[phospho]TIDE', {2: ['acetyl']}, overwrite=True)
        'PEP[acetyl]TIDE'

        >>> add_mods('PEP[phospho]TIDE', {2: ['acetyl']})
        'PEP[phospho][acetyl]TIDE'

        >>> add_mods('PEP[phospho]TIDE', {'u': ['acetyl']})
        '[acetyl]?PEP[phospho]TIDE'

        >>> add_mods('PEPTIDE', {'s': ['[100][200]@C']}, overwrite=True)
        '<[100][200]@C>PEPTIDE'

    """

    stripped_sequence, original_mods = pop_mods(sequence)

    for mod_index, mod in mods.items():

        if not isinstance(mod, list):
            mod = [mod]

        if mod_index in original_mods:
            if overwrite is True:
                original_mods[mod_index] = mod
            else:
                original_mods[mod_index].extend(mod)
        else:
            original_mods[mod_index] = mod

    return _construct_sequence(stripped_sequence, original_mods)


def condense_static_mods(sequence: str) -> str:
    """
    Condenses static modifications into a single modification.

    :param sequence: The peptide sequence.
    :type sequence: str
    :return: The peptide sequence with condensed static modifications.
    :rtype: str

    .. code-block:: python

        # Condenses static modifications into a single modification
        >>> condense_static_mods('<13C><[100]@P>PEPTIDE')
        '<13C>P[100]EP[100]TIDE'

        >>> condense_static_mods('<13C><[100]@P>P[10]EPTIDE')
        '<13C>P[10][100]EP[100]TIDE'

    """

    stripped_sequence, original_mods = pop_mods(sequence)

    static_mods = original_mods.pop('s', [])

    static_mod_dict = parse_static_mods(static_mods)
    for aa in static_mod_dict:
        # get indexes of the amino acid
        indexes = [m.start() for m in re.finditer(aa, stripped_sequence)]
        for index in indexes:
            original_mods.setdefault(index, []).extend(static_mod_dict[aa])

    return add_mods(stripped_sequence, original_mods)


def pop_mods(sequence: str) -> Tuple[str, ModDict]:
    """
    Removes all modifications from the given sequence, returning the unmodified sequence and a dictionary of the
    removed modifications.

    :param sequence: The sequence to be stripped of modifications.
    :type sequence: str

    :return: A tuple containing the unmodified sequence and a dictionary of the removed modifications.
    :rtype: Tuple[str, Dict[int, Any]]

    .. code-block:: python

        # Simply combines the functionality of strip_modifications and get_modifications
        >>> pop_mods('PEP[phospho]TIDE')
        ('PEPTIDE', {2: ['phospho']})

    """

    return strip_mods(sequence), get_mods(sequence)


def strip_mods(sequence: str) -> str:
    """
    Strips all modifications from the given sequence, returning the unmodified sequence.

    :param sequence: The sequence to be stripped of modifications.
    :type sequence: str

    :return: The stripped sequence
    :rtype: str

    .. code-block:: python

        # Removes internal modifications:
        >>> strip_mods('PEP[phospho]TIDE')
        'PEPTIDE'

        # Also removes N and C terminal modifications:
        >>> strip_mods('[Acetyl]-[3.14]PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Also remove labile modifications:
        >>> strip_mods('{1.0}[Acetyl]-[3.14]PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Also remove isotope notations:
        >>> strip_mods('<C13>[Acetyl]-[3.14]PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Using a sequence without modifications will return the same sequence:
        >>> strip_mods('PEPTIDE')
        'PEPTIDE'

        >>> strip_mods('PEP[Formula:[13C]H12]TIDE')
        'PEPTIDE'

        >>> strip_mods('(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K')
        '(?DQ)NGTWEMESNENFEGYMK'

        >>> strip_mods('[1][2]^2?[100]^3-PEP[1]^2TIDE')
        'PEPTIDE'

    """

    # remove {}
    sequence = re.sub(r'{[^}]*}', '', sequence)

    # remove <>
    sequence = re.sub(r'<[^>]*>', '', sequence)

    # Removing modifications in square brackets
    sequence, bracket_mods = map_bracket_content_to_index(sequence)

    if sequence.startswith('?'):
        sequence = sequence[1:]

    return sequence.replace('-', '')


def reverse(sequence: str, swap_terms: bool = False) -> str:
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
        >>> reverse('PEPTIDE')
        'EDITPEP'

        >>> reverse('<13C>PEPTIDE')
        '<13C>EDITPEP'

        # For modified sequences, the modifications are preserved on the associated residues:
        >>> reverse('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]')
        '[Acetyl]-EDITP[phospho]EP[phospho]-[Amide]'

        # If swap_terms is True, the N- and C-terminal modifications will be swapped too:
        >>> reverse('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]', swap_terms=True)
        '[Amide]-EDITP[phospho]EP[phospho]-[Acetyl]'

    """

    mods = get_mods(sequence)
    stripped_sequence = strip_mods(sequence)[::-1]

    if is_ambiguous(stripped_sequence):
        raise ValueError("Cannot reverse ambiguous sequence")

    mod_reverse = {}
    for mod_index in mods:
        if isinstance(mod_index, int):
            mod_reverse[len(stripped_sequence) - mod_index - 1] = mods[mod_index]
        elif mod_index == 'n' and swap_terms is True:
            mod_reverse['c'] = mods[mod_index]
        elif mod_index == 'c' and swap_terms is True:
            mod_reverse['n'] = mods[mod_index]
        else:
            mod_reverse[mod_index] = mods[mod_index]

    return _construct_sequence(stripped_sequence, mod_reverse)


def shuffle(sequence: str, seed: int = None) -> str:
    """
    Shuffles the sequence, while preserving the position of any modifications.

    :param sequence: The sequence to be shuffled.
    :type sequence: str
    :param seed: Seed for the random number generator.
    :type seed: int

    :return: The shuffled sequence with modifications preserved.
    :rtype: str

    .. code-block:: python

        # For unmodified sequences, the result is a random permutation of the original sequence:
        >>> shuffle('PEPTIDE', seed=0)
        'IPEPDTE'

        # For modified sequences, the modifications are preserved on the associated residues:
        >>> shuffle('[Acetyl]-PEPTIDE', seed=0)
        '[Acetyl]-IPEPDTE'

        >>> shuffle('<13C>PEPTIDE', seed=0)
        '<13C>IPEPDTE'

    """

    if seed is not None:
        random.seed(seed)

    stripped_sequence, mods = pop_mods(sequence)

    if is_ambiguous(stripped_sequence):
        raise ValueError("Cannot shuffle ambiguous sequence")

    n_term_mods = mods.pop('n', [])
    c_term_mods = mods.pop('c', [])
    isotope_mods = mods.pop('i', [])
    static_mods = mods.pop('s', [])
    labile_mods = mods.pop('l', [])
    unknown_mods = mods.pop('u', [])

    internally_modified_sequence = add_mods(stripped_sequence, mods)
    components = split(internally_modified_sequence)
    random.shuffle(components)
    shuffled_sequence = ''.join(components)

    updated_mods = {
        'n': n_term_mods,
        'c': c_term_mods,
        'i': isotope_mods,
        's': static_mods,
        'l': labile_mods,
        'u': unknown_mods
    }

    return add_mods(shuffled_sequence, updated_mods)


def shift(sequence: str, n: int) -> str:
    """
    Shifts the sequence to the left by a given number of positions, while preserving the position of any modifications.

    :param sequence: The sequence to be shifted.
    :type sequence: str
    :param n: The number of positions to shift the sequence to the left.
    :type n: int

    :return: The shifted sequence with modifications preserved.
    :rtype: str

    .. code-block:: python

        >>> shift('PEPTIDE', 2)
        'PTIDEPE'

        >>> shift('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]', 2)
        '[Acetyl]-P[phospho]TIDEP[phospho]E-[Amide]'

        # Shifting by 0 or length of sequence positions returns the original sequence:
        >>> shift('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]', 7)
        '[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]'
        >>> shift('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]', 0)
        '[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]'

        # Shifting by a negative number shifts the sequence to the right:
        >>> shift('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]', -2)
        '[Acetyl]-DEP[phospho]EP[phospho]TI-[Amide]'

        >>> shift('<13C>PEPTIDE', 2)
        '<13C>PTIDEPE'

    """

    stripped_sequence, mods = pop_mods(sequence)

    if is_ambiguous(stripped_sequence):
        raise ValueError("Cannot shift ambiguous sequence")

    seq_len = len(stripped_sequence)
    effective_shift = n % seq_len
    shifted_stripped_sequence = stripped_sequence[effective_shift:] + stripped_sequence[:effective_shift]

    shifted_sequence_mods = {}
    for mod_index in mods:
        if isinstance(mod_index, int):
            shifted_sequence_mods[(mod_index - effective_shift) % seq_len] = mods[mod_index]
        else:
            shifted_sequence_mods[mod_index] = mods[mod_index]

    return _construct_sequence(shifted_stripped_sequence, shifted_sequence_mods)


def span_to_sequence(sequence: str, span: Span) -> str:
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
        >>> span_to_sequence('[Acetyl]-P[1.0]EPTIDE', (0, 4, 0))
        '[Acetyl]-P[1.0]EPT'

        >>> span_to_sequence('PEPTIDE[1.0]-[Amide]', (4, 7, 0))
        'IDE[1.0]-[Amide]'

        >>> span_to_sequence('<13C>PEPTIDE[1.0]-[Amide]', (1, 6, 0))
        '<13C>EPTID'

    """

    stripped_sequence, mods = pop_mods(sequence)
    return span_to_sequence_fast(stripped_sequence, mods, span)


def span_to_sequence_fast(stripped_sequence: str, mods: ModDict, span: Span) -> str:
    """
    Faster version of span_to_sequence

    :param stripped_sequence: The stripped sequence.
    :type stripped_sequence: str
    :param mods: The modifications for sequence
    :type mods: Dict[int, Union[str, float, int]]
    :param span: A tuple representing the span of the subsequence to be extracted.
    :type span: Tuple[int, int, int]

    :return: The subsequence of the input sequence defined by the span.
    :rtype: str

    .. code-block:: python

        >>> span_to_sequence_fast('PEPTIDE', {}, (0, 4, 0))
        'PEPT'

        >>> span_to_sequence_fast('PEPTIDE', {}, (1, 6, 0))
        'EPTID'

        >>> span_to_sequence_fast('PEPTIDE', {'n': [1], 'c': [2], 2: [2]}, (4, 7, 0))
        'IDE-[2]'

        >>> span_to_sequence_fast('PEPTIDE', {'n': [1], 'c': [2], 2: [2]}, (0, 3, 0))
        '[1]-PEP[2]'

        >>> span_to_sequence_fast('TIDERTIDEKTIDE', {9: [1], 'c': [2]}, (5, 14, 1))
        'TIDEK[1]TIDE-[2]'

    """

    base_sequence = stripped_sequence[span[0]:span[1]]

    # if peptide is unmodified skip the modification parsing
    if len(mods) == 0:
        return base_sequence

    new_mods = {}
    for k in mods:
        if isinstance(k, int):
            if span[0] <= k < span[1]:
                new_mods[k - span[0]] = mods[k]
        else:
            new_mods[k] = mods[k]

    if span[0] != 0:
        _ = new_mods.pop('n', None)

    if span[1] != len(stripped_sequence):
        _ = new_mods.pop('c', None)

    # if peptide is unmodified skip the modification parsing
    if len(mods) == 0:
        return base_sequence

    # update key values
    return _construct_sequence(base_sequence, new_mods)


def split(sequence: str) -> List[str]:
    """
    Splits sequence into a list of amino acids, preserving modifications.

    :param sequence: The sequence to be split.
    :type sequence: str

    :return: A list of amino acids, preserving modifications.
    :rtype: List[str]

    .. code-block:: python

        >>> split('PEPTIDE')
        ['P', 'E', 'P', 'T', 'I', 'D', 'E']

        >>> split('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]')
        ['[Acetyl]-P[phospho]', 'E', 'P[phospho]', 'T', 'I', 'D', 'E-[Amide]']

        >>> split('<C13>PEP[1]TIDE')
        ['<C13>P', '<C13>E', '<C13>P[1]', '<C13>T', '<C13>I', '<C13>D', '<C13>E']

    """

    mods = get_mods(sequence)
    stripped_sequence = strip_mods(sequence)

    if is_ambiguous(stripped_sequence):
        raise ValueError("Cannot split ambiguous sequence")

    n_term_mods = mods.pop('n', [])
    c_term_mods = mods.pop('c', [])
    isotope_mods = mods.pop('i', [])
    static_mods = mods.pop('s', [])
    labile_mods = mods.pop('l', [])
    unknown_mods = mods.pop('u', [])

    comps = []
    for i, aa in enumerate(stripped_sequence):
        tmp_mods = {}
        if i == 0 and n_term_mods:
            tmp_mods['n'] = n_term_mods
        if i == len(stripped_sequence) - 1 and c_term_mods:
            tmp_mods['c'] = c_term_mods
        if i in mods:
            tmp_mods[0] = mods[i]

        tmp_mods['i'] = isotope_mods
        tmp_mods['s'] = static_mods
        tmp_mods['l'] = labile_mods
        tmp_mods['u'] = unknown_mods

        comps.append(_construct_sequence(aa, tmp_mods))

    return comps


def count_components(sequence: str) -> Counter:
    """
    Counts the occurrences of each amino acid in the input sequence.

    :param sequence: The sequence to be counted.
    :type sequence: str

    :return: A Counter object containing the occurrences of each amino acid in the input sequence.
    :rtype: Counter

    .. code-block:: python

        >>> count_components('PEPTIDE')
        Counter({'P': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1})

        >>> count_components('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]')
        Counter({'P[phospho]': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1})

        >>> count_components('<13C>PE[3.14]PTIDE')
        Counter({'P': 2, 'E[3.14]': 1, 'T': 1, 'I': 1, 'D': 1, 'E': 1})

    """
    sequence, mods = pop_mods(sequence)

    _ = mods.pop('l', [])
    _ = mods.pop('s', [])
    _ = mods.pop('i', [])
    _ = mods.pop('n', [])
    _ = mods.pop('c', [])
    _ = mods.pop('u', [])

    internal_modified_sequence = _construct_sequence(sequence, mods)
    components = split(internal_modified_sequence)
    return collections.Counter(components)


def is_subsequence(subsequence: str, sequence: str, order: bool = True) -> bool:
    """
    Checks if the input subsequence is a subsequence of the input sequence. If order is True, the subsequence must be in
    the same order as in the sequence. If order is False, the subsequence can be in any order.

    :param subsequence: The subsequence to be checked.
    :type subsequence: str
    :param sequence: The sequence to be checked.
    :type sequence: str
    :param order: If True, the subsequence must be in the same order as in the sequence.
    :type order: bool

    :return: True if the subsequence is a subsequence of the sequence, False otherwise.
    :rtype: bool

    .. code-block:: python

        >>> is_subsequence('PEP', 'PEPTIDE')
        True

        >>> is_subsequence('PET', 'PEPTIDE', order=False)
        True

        >>> is_subsequence('PET', 'PEPTIDE', order=True)
        False

        >>> is_subsequence('<13C>PEP', '<13C>PEPTIDE', order=True)
        True

        >>> is_subsequence('<13C>PEP[1.0]', '<13C>PEP[1.0]TIDE', order=True)
        True

        >>> is_subsequence('<13C>PEP', '<13C>PEP[1.0]TIDE', order=True)
        False

    """

    if order is True:
        return len(find_subsequence_indices(sequence, subsequence)) != 0
    else:

        _, sequence_mods = pop_mods(sequence)
        sort_mods(sequence_mods)
        _, subsequence_mods = pop_mods(subsequence)
        sort_mods(subsequence_mods)

        if sequence_mods.get('n') != subsequence_mods.get('n'):
            return False

        if sequence_mods.get('c') != subsequence_mods.get('c'):
            return False

        if sequence_mods.get('i') != subsequence_mods.get('i'):
            return False

        if sequence_mods.get('s') != subsequence_mods.get('s'):
            return False

        if sequence_mods.get('u') != subsequence_mods.get('u'):
            return False

        subsequence_counts = count_components(subsequence)
        sequence_counts = count_components(sequence)
        return all(subsequence_counts[aa] <= sequence_counts[aa] for aa in subsequence_counts)


def sort_mods(mods: ModDict, sort_function: Callable[[str], str] = lambda x: str(x)) -> None:
    """
    Sorts the modifications in the input dictionary using the provided sort function.

    :param mods: The modifications to be sorted.
    :type mods: Dict[int, List[ModValue]]
    :param sort_function: The sort function to be used. Defaults to identity function.
    :type sort_function: Callable[[str], str]

    :return: None
    :rtype: None

    .. code-block:: python

        >>> mods = {1: ['phospho', 1], 2: ['phospho']}
        >>> sort_mods(mods)
        >>> mods
        {1: [1, 'phospho'], 2: ['phospho']}

    """

    for k in mods:
        mods[k].sort(key=sort_function)


def sort(sequence: str, sort_function: Callable[[str], str] = lambda x: x[0]) -> str:
    """
    Sorts the input sequence using the provided sort function. Terminal sequence are kept in place.

    :param sequence: The sequence to be sorted.
    :type sequence: str

    :param sort_function: The sort function to be used. Defaults to identity function.
    :type sort_function: Callable[[str], str]

    :return: The input sequence sorted using the provided sort function.
    :rtype: str

    .. code-block:: python

        >>> sort('PEPTIDE')
        'DEEIPPT'

        >>> sort('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]')
        '[Acetyl]-DEEIP[phospho]P[phospho]T-[Amide]'

        >>> sort('[Acetyl]-P[1][phospho]EP[phospho]TIDE-[Amide]')
        '[Acetyl]-DEEIP[1][phospho]P[phospho]T-[Amide]'

    """

    stripped_sequence, mods = pop_mods(sequence)

    sort_mods(mods)

    residue_mods = {k: v for k, v in mods.items() if isinstance(k, int)}
    other_mods = {k: v for k, v in mods.items() if isinstance(k, str)}

    internal_modified_sequence = _construct_sequence(stripped_sequence, residue_mods)
    components = split(internal_modified_sequence)
    components.sort(key=sort_function)
    return add_mods(''.join(components), other_mods)


def get_kmers(sequence: str, k: int) -> List[str]:
    """
    Builds kmers from a given sequence.

    :param sequence: The sequence to build kmers from.
    :type sequence: str
    :param k: The length of the kmers.
    :type k: int

    :return: A list of kmers.
    :rtype: list

    .. code-block:: python

        >>> get_kmers('PEPTIDE', 2)
        ['PE', 'EP', 'PT', 'TI', 'ID', 'DE']

        >>> get_kmers('PEPTIDE', 3)
        ['PEP', 'EPT', 'PTI', 'TID', 'IDE']

        >>> get_kmers('<13C>PE[1][2]PTIDE', 3)
        ['<13C>PE[1][2]P', '<13C>E[1][2]PT', '<13C>PTI', '<13C>TID', '<13C>IDE']

        >>> get_kmers('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]', 3)
        ['[Acetyl]-P[phospho]EP[phospho]', 'EP[phospho]T', 'P[phospho]TI', 'TID', 'IDE-[Amide]']

        >>> get_kmers('==PE[2]P', 2)
        ['==', '=P', 'PE[2]', 'E[2]P']


    """

    stripped_sequence, mods = pop_mods(sequence)

    peptide_mods = {k: v for k, v in mods.items() if isinstance(k, int) or k in ['n', 'c']}
    global_mods = {k: v for k, v in mods.items() if isinstance(k, str) and k not in ['n', 'c']}

    internal_modified_sequence = _construct_sequence(stripped_sequence, peptide_mods)
    components = split(internal_modified_sequence)

    kmers = []
    for i in range(len(stripped_sequence) - k + 1):
        kmer = ''.join(components[i:i + k])
        kmers.append(add_mods(kmer, global_mods))

    return kmers


def find_subsequence_indices(sequence: str, subsequence: str, ignore_mods: bool = False) -> List[int]:
    """
    Retrieves all starting indexes of a given subsequence within a sequence.

    :param sequence: The parent sequence.
    :type sequence: str
    :param subsequence: The subsequences to search for.
    :type subsequence: str
    :param ignore_mods: Whether to ignore modifications.
    :type ignore_mods: bool

    :return: A list of starting indexes of the subsequence within the sequence.
    :rtype: List[int]

    .. code-block:: python

        # Find the starting indexes of a subsequence
        >>> find_subsequence_indices("PEPTIDE", "PEP")
        [0]
        >>> find_subsequence_indices("PEPTIDE", "EPT")
        [1]
        >>> find_subsequence_indices("PEPTIDE", "E")
        [1, 6]

        # By default the function will not ignore modifications
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "PEP")
        []
        >>> find_subsequence_indices("<13C>PEPTIDE", "PEP")
        []
        >>> find_subsequence_indices("<13C>PEP[1][Phospho]TIDE", "<13C>PEP[Phospho][1]")
        [0]
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "[Acetyl]-PEP")
        [0]
        >>> find_subsequence_indices("PEPTIDE", "[Acetyl]-PEP")
        []
        >>> find_subsequence_indices("PEPTIDE[1.0]", "IDE")
        []
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE-[Amide]", "[Acetyl]-PEPTIDE-[Amide]")
        [0]

        # If ignore_mods is set to True, the function will ignore modifications
        >>> find_subsequence_indices("[Acetyl]-PEPTIDE", "PEP", ignore_mods=True)
        [0]
        >>> find_subsequence_indices("PEPTIDE", "[Acetyl]-PEP", ignore_mods=True)
        [0]

    """

    unmodified_peptide, peptide_mods = pop_mods(subsequence)
    unmodified_protein, protein_mods = pop_mods(sequence)

    if ignore_mods:
        return [i.start() for i in re.finditer(unmodified_peptide, unmodified_protein, overlapped=True)]

    sort_mods(peptide_mods)
    sort_mods(protein_mods)

    subsequence = span_to_sequence_fast(unmodified_peptide, peptide_mods, (0, len(unmodified_peptide), 0))

    spans = [(i.start(), i.end(), 0) for i in re.finditer(unmodified_peptide, unmodified_protein, overlapped=True)]
    hit_peptides = [span_to_sequence_fast(unmodified_protein, protein_mods, s) for s in spans]

    return [s[0] for s, p in zip(spans, hit_peptides) if p == subsequence]


def coverage(sequence: str, subsequences: List[str], accumulate: bool = False,
             ignore_mods: bool = False) -> List[int]:
    """
    Calculate the sequence coverage given a list of subsequecnes.

    The coverage is represented as a binary list where each position in the protein sequence is marked as 1 if it
    is covered by at least one peptide and 0 otherwise.

    :param sequence: The parent sequence.
    :type sequence: str
    :param subsequences: List of subsequences.
    :type subsequences: List[str]
    :param accumulate: If True, the coverage array will be accumulated, i.e. if a position is covered by more than
                        one subsequence, it will be marked as the sum of the number of subsequences covering it.
                        If False, the coverage array will be binary, i.e. if a position is covered by more than one
                        subsequence, it will be marked as 1.
    :type accumulate: bool
    :param ignore_mods: Whether to ignore modifications when calcualting the coverage.
    :type ignore_mods: bool

    :return: The sequence coverage array.
    :rtype: List[int]

    .. code-block:: python

        >>> coverage("PEPTIDE", ["PEP"])
        [1, 1, 1, 0, 0, 0, 0]

        >>> coverage("PEPTIDE", ["PEP", "EPT"])
        [1, 1, 1, 1, 0, 0, 0]

        # If accumulate is True, overlapping indecies will be accumulated
        >>> coverage("PEPTIDE", ["PEP", "EPT"], accumulate=True)
        [1, 2, 2, 1, 0, 0, 0]

    """

    cov_arr = [0] * len(sequence)
    for peptide in subsequences:
        peptide_indexes = find_subsequence_indices(sequence, peptide, ignore_mods=ignore_mods)
        for peptide_index in peptide_indexes:
            if accumulate:
                cov_arr[peptide_index:peptide_index + len(peptide)] = \
                    [x + 1 for x in cov_arr[peptide_index:peptide_index + len(peptide)]]
            else:
                cov_arr[peptide_index:peptide_index + len(peptide)] = [1] * len(peptide)

    return cov_arr


def percent_coverage(sequence: str, subsequences: List[str], ignore_mods: bool = False) -> float:
    """
    Calculates the coverage given a list of subsequences.

    :param sequence: The parent sequence.
    :type sequence: str
    :param subsequences: List of subsequences.
    :type subsequences: List[str]
    :param ignore_mods: Whether to ignore modifications when calcualting the coverage.
    :type ignore_mods: bool

    :return: The percent coverage.
    :rtype: float

    .. code-block:: python

        >>> percent_coverage("PEPTIDE", ["PEP"])
        0.42857142857142855

        >>> percent_coverage("PEPTIDE", ["PEP", "EPT"])
        0.5714285714285714

    """

    cov_arr = coverage(sequence, subsequences, accumulate=False, ignore_mods=ignore_mods)

    if len(cov_arr) == 0:
        return 0

    return sum(cov_arr) / len(cov_arr)
