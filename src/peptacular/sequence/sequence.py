"""
Sequence.py - Amino Acid Sequence Manipulation with Modifications

Modifications can be a str, int, or float

Residue modifications are mapped to the unmodified sequence by their index.

Special modifications are mapped to the sequence by the following:
    - N-terminal modifications are mapped to the index 'nterm'
    - C-terminal modifications are mapped to the index 'cterm'
    - Isotopic modifications are mapped to the index 'isotope'
    - Static modifications are mapped to the index 'static'
    - Labile modifications are mapped to the index 'labile'
    - Unknown modifications are mapped to the index 'unknown'
    - Internals are mapped to the index 'interval'


Global Mods and labile mods are popped from the sequence and returned as a tuple with the stripped sequence. As
a result they can be positioned anywhere in the sequence. Maybe add a check to see if they are at the beginning or
end of the sequence and if not raise an error.

"""

# TODO: add mode to all mods additions: overwrite, append, skip...

from __future__ import annotations

import typing

import regex as re
from typing import Dict, List, Union, Tuple, Counter, Callable

from peptacular.sequence.proforma import parse, ProFormaAnnotation, serialize, Mod
from peptacular.spans import Span

from peptacular.types import ModDict, ModValue
from peptacular.input_parser import fix_list_of_mods


def parse_single_sequence(sequence: str) -> ProFormaAnnotation:
    annotations = parse(sequence)
    if len(annotations) != 1:
        raise ValueError(f"Invalid sequence: {sequence}")

    return annotations[0]


def sequence_length(sequence: str | ProFormaAnnotation) -> int:
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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    return len(annotation)


def is_ambiguous(sequence: str | ProFormaAnnotation) -> bool:
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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    return annotation.is_ambiguous()


def get_mods(sequence: str | ProFormaAnnotation) -> ModDict:
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
        {2: [Mod('Phospho', 1)], 3: [Mod(1, 1)], 6: [Mod(3.14, 1)]}

        >>> get_mods('PEP[Phospho][1.0]TIDE')
        {2: [Mod('Phospho', 1), Mod(1.0, 1)]}

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    return annotation.mod_dict()


def add_mods(sequence: str | ProFormaAnnotation,
             mods: Dict[str, List[ModValue] | ModValue],
             append: bool = True) -> str:
    """
    Adds modifications to the given peptide sequence.

    :param sequence: Unmodified peptide sequence.
    :type sequence: str
    :param mods: Dictionary with indices of the modified amino acids and corresponding modifications.
    :type mods: Dict[int, Any]
    :param append: If True, the modifications will be appended to the existing modifications.
                      If False, the existing modifications will be replaced.
    :type append: bool

    :return: The peptide sequence with the specified modifications.
    :rtype: str

    .. code-block:: python

        # Add internal modifications to an unmodified peptide
        >>> add_mods('PEPTIDE', {2: [Mod('phospho', 1)]})
        'PEP[phospho]TIDE'

        # Can also add N and C terminal modifications
        >>> add_mods('PEPTIDE', {'nterm': [Mod('Acetyl', 1)], 6: [Mod(1.234, 1)], 'cterm': [Mod('Amide', 1)]})
        '[Acetyl]-PEPTIDE[1.234]-[Amide]'

        # Empty modification dicts will simply return the input sequence
        >>> add_mods('PEPTIDE', {})
        'PEPTIDE'

        >>> add_mods('PEP[phospho]TIDE', {2: [Mod('acetyl', 1)]}, append=False)
        'PEP[acetyl]TIDE'

        >>> add_mods('PEP[phospho]TIDE', {2: [Mod('acetyl', 1)]})
        'PEP[phospho][acetyl]TIDE'

        >>> add_mods('PEP[phospho]TIDE', {'unknown': [Mod('acetyl', 1)]})
        '[acetyl]?PEP[phospho]TIDE'

        >>> add_mods('PEPTIDE', {'static': [Mod('[100][200]@C', 1)]}, append=False)
        '<[100][200]@C>PEPTIDE'

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    mods = {k: fix_list_of_mods(v) for k, v in mods.items()}
    annotation.add_mod_dict(mods, append=append)
    return serialize([annotation])


def condense_static_mods(sequence: str | ProFormaAnnotation) -> str:
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

        >>> condense_static_mods('PEPTIDE')
        'PEPTIDE'

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    return annotation.condense_static_mods().serialize()


def pop_mods(sequence: str | ProFormaAnnotation) -> Tuple[str, ModDict]:
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
        ('PEPTIDE', {2: [Mod('phospho', 1)]})

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    return annotation.sequence, annotation.mod_dict()


def strip_mods(sequence: str | ProFormaAnnotation) -> str:
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
        >>> strip_mods('[Acetyl]-PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Also remove labile modifications:
        >>> strip_mods('{1.0}[Acetyl]-PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Also remove isotope notations:
        >>> strip_mods('<C13>[Acetyl]-PEPTIDE[1.234]-[Amide]')
        'PEPTIDE'

        # Using a sequence without modifications will return the same sequence:
        >>> strip_mods('PEPTIDE')
        'PEPTIDE'

        >>> strip_mods('PEP[Formula:[13C]H12]TIDE')
        'PEPTIDE'

        >>> strip_mods('(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K')
        'DQNGTWEMESNENFEGYMK'

        >>> strip_mods('[1][2]^2?[100]^3-PEP[1]^2TIDE')
        'PEPTIDE'

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    return annotation.sequence


def reverse(sequence: str | ProFormaAnnotation, swap_terms: bool = False) -> str:
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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)

    else:
        annotation = sequence

    reversed_annotation = annotation.reverse(swap_terms=swap_terms)
    return reversed_annotation.serialize()


def shuffle(sequence: str | ProFormaAnnotation, seed: int = None) -> str:
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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)

    else:
        annotation = sequence

    shifted_annotation = annotation.shuffle(seed)
    return shifted_annotation.serialize()


def shift(sequence: str | ProFormaAnnotation, n: int) -> str:
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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    shifted_annotation = annotation.shift(n)
    return shifted_annotation.serialize()


def span_to_sequence(sequence: str | ProFormaAnnotation, span: Span) -> str:
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

        >>> span_to_sequence('(PEPT)IDE', (1, 6, 0))
        '(EPT)ID'

    """

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    sliced_annotation = annotation.slice(span[0], span[1])
    return sliced_annotation.serialize()


def split(sequence: str | ProFormaAnnotation) -> List[str]:
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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)

    else:
        annotation = sequence

    return [annot.serialize() for annot in annotation.split()]


def count_residues(sequence: str) -> typing.Counter:
    """
    Counts the occurrences of each amino acid in the input sequence.

    :param sequence: The sequence to be counted.
    :type sequence: str

    :return: A Counter object containing the occurrences of each amino acid in the input sequence.
    :rtype: Counter

    .. code-block:: python

        >>> count_residues('PEPTIDE')
        Counter({'P': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1})

        >>> count_residues('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]')
        Counter({'[Acetyl]-P[phospho]': 1, 'E': 1, 'P[phospho]': 1, 'T': 1, 'I': 1, 'D': 1, 'E-[Amide]': 1})

        >>> count_residues('<13C>PE[3.14]PTIDE')
        Counter({'<13C>P': 2, '<13C>E[3.14]': 1, '<13C>T': 1, '<13C>I': 1, '<13C>D': 1, '<13C>E': 1})

    """
    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)
    else:
        annotation = sequence

    return annotation.count_residues()


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
        _sort_mods(sequence_mods)
        _, subsequence_mods = pop_mods(subsequence)
        _sort_mods(subsequence_mods)

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

        subsequence_counts = count_residues(subsequence)
        sequence_counts = count_residues(sequence)
        return all(subsequence_counts[aa] <= sequence_counts[aa] for aa in subsequence_counts)


def _sort_mods(mods: ModDict, sort_function: Callable[[str], str] = lambda x: str(x)) -> None:
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
        >>> _sort_mods(mods)
        >>> mods
        {1: [1, 'phospho'], 2: ['phospho']}

    """

    for k in mods:
        mods[k].sort(key=sort_function)


def sort(sequence: str) -> str:
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

    if isinstance(sequence, str):
        annotation = parse_single_sequence(sequence)

    else:
        annotation = sequence

    sorted_annotation = annotation.sort_residues()
    return sorted_annotation.serialize()


def find_subsequence_indices(sequence: str | ProFormaAnnotation,
                             subsequence: str | ProFormaAnnotation,
                             ignore_mods: bool = False) -> List[int]:
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

    if isinstance(sequence, str):
        sequence = parse_single_sequence(sequence)

    if isinstance(subsequence, str):
        subsequence = parse_single_sequence(subsequence)

    if ignore_mods:
        sequence = sequence.strip()
        subsequence = subsequence.strip()

    return subsequence.find_indices(sequence)


def coverage(sequence: str | ProFormaAnnotation,
             subsequences: List[str | ProFormaAnnotation],
             accumulate: bool = False,
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

    if isinstance(sequence, str):
        sequence = parse_single_sequence(sequence)

    if isinstance(subsequences[0], str):
        subsequences = [parse_single_sequence(sub) for sub in subsequences]

    cov_arr = [0] * sequence_length(sequence)
    for subsequence in subsequences:
        peptide_indexes = find_subsequence_indices(sequence, subsequence, ignore_mods=ignore_mods)
        for peptide_index in peptide_indexes:
            if accumulate:
                cov_arr[peptide_index:peptide_index + sequence_length(subsequence)] = \
                    [x + 1 for x in cov_arr[peptide_index:peptide_index + sequence_length(subsequence)]]
            else:
                cov_arr[peptide_index:peptide_index + sequence_length(subsequence)] = [1] * sequence_length(subsequence)

    return cov_arr


def percent_coverage(sequence: str | ProFormaAnnotation,
                     subsequences: List[str | ProFormaAnnotation],
                     ignore_mods: bool = False) -> float:
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


def convert_ip2_sequence(sequence: str) -> str:
    """
    Converts a IP2-Like sequence to a proforma2.0 compatible sequence.

    :param sequence: The sequence to be converted.
    :type sequence: str

    :return: Proforma2.0 compatable sequence.
    :rtype: str

    .. code-block:: python

        >>> convert_ip2_sequence('K.PEP(phospho)TIDE.K')
        'PEP[phospho]TIDE'

        >>> convert_ip2_sequence('K.(1)PEP(phospho)TIDE.K')
        '[1]-PEP[phospho]TIDE'

        >>> convert_ip2_sequence('K.PEPTIDE(2).K')
        'PEPTIDE[2]'

        >>> convert_ip2_sequence('K.PEPTIDE(2)(3).K')
        'PEPTIDE[2]-[3]'

        >>> convert_ip2_sequence('-.(1)PEP(phospho)TIDE(2)(3).-')
        '[1]-PEP[phospho]TIDE[2]-[3]'

        >>> convert_ip2_sequence('P')
        'P'

        >>> convert_ip2_sequence('')
        ''

        >>> convert_ip2_sequence('PEPTIDE')
        'PEPTIDE'
    """

    # Use regex to check if sequence starts and ends with the specified pattern
    if re.match(r'^([A-Z]|-)\..*\.([A-Z]|-)$', sequence):
        # If it matches, remove the leading and trailing characters (first and last two characters)
        sequence = sequence[2:-2]

    # Step 2: Replace () with []
    sequence = re.sub(r'\(([^)]+)\)', r'[\1]', sequence)

    # Step 3: Handle modifications at the start
    sequence = re.sub(r'^\[(\d+)\]', r'[\1]-', sequence)

    # Step 4: Convert consecutive modifications to use a dash
    sequence = re.sub(r'\]\[', r']-[', sequence)

    return sequence
