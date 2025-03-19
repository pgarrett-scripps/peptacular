"""
Sequence.py - A module for manipulating peptide sequences with modifications.

Modifications can be a str, int, float, or Mod object.

Residue modifications are mapped to the unmodified sequence by their index.

Special modifications are mapped to the sequence by the following:
    - N-terminal modifications are mapped to the index 'nterm'
    - C-terminal modifications are mapped to the index 'cterm'
    - Isotopic modifications are mapped to the index 'isotope'
    - Static modifications are mapped to the index 'static'
    - Labile modifications are mapped to the index 'labile'
    - Unknown modifications are mapped to the index 'unknown'
    - Intervals are mapped to the index 'interval'


Global Mods and labile mods are popped from the sequence and returned as a tuple with the stripped sequence. As
a result they can be positioned anywhere in the sequence. Maybe add a check to see if they are at the beginning or
end of the sequence and if not raise an error.

"""

from typing import Counter as CounterType, Optional
from typing import Dict, List, Tuple, Callable, Union

import regex as re

from peptacular.constants import ORDERED_AMINO_ACIDS
from peptacular.proforma.proforma_parser import parse, ProFormaAnnotation, serialize, MultiProFormaAnnotation
from peptacular.spans import Span
from peptacular.proforma.input_convert import ModDict, fix_list_of_mods, fix_intervals_input


def sequence_to_annotation(sequence: str) -> ProFormaAnnotation:
    """
    Parses a peptide sequence with modifications and returns a ProFormaAnnotation object.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A ProFormaAnnotation object representing the input sequence.
    :rtype: ProFormaAnnotation

    """
    annotation = parse(sequence)

    if isinstance(annotation, MultiProFormaAnnotation):
        raise ValueError(f"Invalid sequence: {sequence}")

    return annotation


def sequence_length(sequence: Union[str, ProFormaAnnotation]) -> int:
    """
    Compute the length of the peptide sequence based on the unmodified sequence.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence length.
    :rtype: int

    .. code-block:: python

        # The length of the unmodified sequence is the same as the length of the string
        >>> sequence_length("PEPTIDE")
        7

        # Modifications are not counted in the sequence length
        >>> sequence_length("[Acetyl]-PEP[1.2345]TID[3.14]E-[Amide]")
        7
        >>> sequence_length("<C13>PEPTIDE[1.2345]")
        7
        >>> sequence_length("(?PE)PTIDE[1.2345]")
        7

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return len(annotation)


def is_ambiguous(sequence: Union[str, ProFormaAnnotation]) -> bool:
    """
    Check if the sequence contains ambiguous amino acids.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: True if the sequence contains ambiguous amino acids, False otherwise.
    :rtype: bool

    .. code-block:: python

        # Any intervals will be considered ambiguous
        >>> is_ambiguous("(?PE)PTIDE")
        True

        # Any unknown modifications will also be considered ambiguous
        >>> is_ambiguous("[Oxidation]?PEPTIDE")
        True

        # Unmodified sequences and explicit modifications are not considered ambiguous
        >>> is_ambiguous("PEPTIDE")
        False

        >>> is_ambiguous("PEPTIDE[Oxidation]")
        False

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return annotation.contains_sequence_ambiguity()


def is_modified(sequence: Union[str, ProFormaAnnotation]) -> bool:
    """
    Check if the sequence contains any modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: True if the sequence contains any modifications, False otherwise.
    :rtype: bool

    .. code-block:: python

        # Any modifications will return True
        >>> is_modified("PEP[Phospho]TIDE")
        True

        # Unmodified sequences will return False
        >>> is_modified("PEPTIDE")
        False

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return annotation.has_mods()


def get_mods(sequence: Union[str, ProFormaAnnotation]) -> ModDict:
    """
    Parses a sequence with modifications and returns a dictionary where keys represent the position/type of the
    modifications.

    Internal modifications are mapped to the unmodified sequence by their index in the unmodified sequence.

    Special modifications are mapped to the sequence by the following:
    - N-terminal modifications are mapped to the index 'nterm'
    - C-terminal modifications are mapped to the index 'cterm'
    - Isotopic modifications are mapped to the index 'isotope'
    - Static modifications are mapped to the index 'static'
    - Labile modifications are mapped to the index 'labile'
    - Unknown modifications are mapped to the index 'unknown'
    - Intervals are mapped to the index 'interval'
    - charge state is mapped to the index 'charge'
    - charge adducts are mapped to the index 'charge_adducts'

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A dictionary with the modifications
    :rtype: ModDict

    .. code-block:: python

        # All modifications will be returned as Mod objects which contain the modification value and multiplier
        >>> get_mods('PEP[Phospho]T[1]IDE[-3.14]')
        {2: [Mod('Phospho', 1)], 3: [Mod(1, 1)], 6: [Mod(-3.14, 1)]}

        >>> get_mods('PEP[Phospho][1.0]TIDE')
        {2: [Mod('Phospho', 1), Mod(1.0, 1)]}

        # N-terminal modifications are mapped to the index 'nterm'
        >>> get_mods('[Acetyl]-PEPTIDE')
        {'nterm': [Mod('Acetyl', 1)]}

        # C-terminal modifications are mapped to the index 'cterm'
        >>> get_mods('PEPTIDE-[Amide]')
        {'cterm': [Mod('Amide', 1)]}

        # Isotopic modifications are mapped to the index 'isotope'
        >>> get_mods('<13C>PEPTIDE')
        {'isotope': [Mod('13C', 1)]}

        # Static modifications are mapped to the index 'static'
        >>> get_mods('<[+1.234]@P>PEPTIDE')
        {'static': [Mod('[+1.234]@P', 1)]}

        # Labile modifications are mapped to the index 'labile'
        >>> get_mods('{Glycan:Hex}PEPTIDE')
        {'labile': [Mod('Glycan:Hex', 1)]}

        # Unknown modifications are mapped to the index 'unknown'
        >>> get_mods('[Phospho]^3?PEPTIDE')
        {'unknown': [Mod('Phospho', 3)]}

        # Intervals are mapped to the index 'interval'
        >>> get_mods('PEP(TI)[Phospho]DE')
        {'intervals': [Interval(3, 5, False, [Mod('Phospho', 1)])]}

        # Charge state is mapped to the index 'charge'
        >>> get_mods('PEPTIDE/+2')
        {'charge': 2}

        # Charge adducts are mapped to the index 'charge_adducts'
        >>> get_mods('PEPTIDE/+2[+2Na+,-H+]')
        {'charge': 2, 'charge_adducts': [Mod('+2Na+,-H+', 1)]}

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return annotation.mod_dict()


def add_mods(sequence: Union[str, ProFormaAnnotation],
             mods: Dict,
             append: bool = True,
             include_plus: bool = False) -> str:
    """
    Adds modifications to the given sequence. The modifications can be of type Mod, str, int, or float, and can be
    a single value or a list of values. The modifications will be added to the sequence in the order they are provided.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param mods: Dictionary representing the modifications to be added to the sequence.
    :type mods: Dict
    :param append: If True, the modifications will be appended to the existing modifications.
                      If False, the existing modifications will be replaced. Defaults to True.
    :type append: bool
    :param include_plus: If True, the modifications will be serialized with a '+' sign for positive values.
    :type include_plus: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The peptide sequence with the specified modifications.
    :rtype: str

    .. code-block:: python

        >>> from peptacular.proforma.proforma_dataclasses import Mod

        # Add internal modifications to an unmodified peptide
        >>> add_mods('PEPTIDE', {2: [Mod('phospho', 1)]})
        'PEP[phospho]TIDE'

        # Can also add N and C terminal modifications
        >>> add_mods('PEPTIDE', {'nterm': 'Acetyl', 6: 1.234, 'cterm': 'Amide'})
        '[Acetyl]-PEPTIDE[1.234]-[Amide]'

        >>> add_mods('PEPTIDE', {'nterm': 'Acetyl', 6: 1.234, 'cterm': 'Amide'}, include_plus=True)
        '[Acetyl]-PEPTIDE[+1.234]-[Amide]'

        # Can also add isotopic modifications
        >>> add_mods('PEPTIDE', {'isotope': ['13C', '15N']})
        '<13C><15N>PEPTIDE'

        # Can also add static modifications
        >>> add_mods('PEPTIDE', {'static': '[+1.234]@P'})
        '<[+1.234]@P>PEPTIDE'

        # Can also add labile modifications
        >>> add_mods('PEPTIDE', {'labile': 'Glycan:Hex'})
        '{Glycan:Hex}PEPTIDE'

        # Can also add unknown modifications
        >>> add_mods('PEPTIDE', {'unknown': Mod('Phospho', 3)})
        '[Phospho]^3?PEPTIDE'

        # Can also add intervals
        >>> add_mods('PEPTIDE', {'intervals': (3, 5, False, 'Phospho')})
        'PEP(TI)[Phospho]DE'

        # Can also add charge state
        >>> add_mods('PEPTIDE', {'charge': 2})
        'PEPTIDE/2'

        # Can also add charge adducts
        >>> add_mods('PEPTIDE', {'charge': 2, 'charge_adducts': '+2Na+,-H+'})
        'PEPTIDE/2[+2Na+,-H+]'

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    for k in mods:
        if k == 'charge':
            continue
        if k == 'intervals':
            mods[k] = fix_intervals_input(mods[k])
        else:
            mods[k] = fix_list_of_mods(mods[k])

    annotation.add_mod_dict(mods, append=append)
    return serialize(annotation, include_plus)


def condense_static_mods(sequence: Union[str, ProFormaAnnotation], include_plus: bool = False) -> str:
    """
    Condenses static modifications into internal modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The peptide sequence with condensed static modifications.
    :rtype: str

    .. code-block:: python

        # Condenses static modifications to specified internal modifications
        >>> condense_static_mods('<13C><[100]@P>PEPTIDE')
        '<13C>P[100]EP[100]TIDE'

        # If residue is already modified, the static modification will be appended
        >>> condense_static_mods('<13C><[100]@P>P[10]EPTIDE')
        '<13C>P[10][100]EP[100]TIDE'

        # Example for unmodified sequences
        >>> condense_static_mods('PEPTIDE')
        'PEPTIDE'

        # Example for N-Term static modifications
        >>> condense_static_mods('<[Oxidation]@N-Term>PEPTIDE')
        '[Oxidation]-PEPTIDE'

        # Example for C-Term static modifications
        >>> condense_static_mods('<[Oxidation]@C-Term>PEPTIDE')
        'PEPTIDE-[Oxidation]'

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return annotation.condense_static_mods().serialize(include_plus)


def pop_mods(sequence: Union[str, ProFormaAnnotation]) -> Tuple[str, ModDict]:
    """
    Removes all modifications from the given sequence, returning the unmodified sequence and a dictionary of the
    removed modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A tuple containing the unmodified sequence and a dictionary of the removed modifications.
    :rtype: Tuple[str, ModDict]

    .. code-block:: python

        # Simply combines the functionality of strip_modifications and get_modifications
        >>> pop_mods('PEP[phospho]TIDE')
        ('PEPTIDE', {2: [Mod('phospho', 1)]})

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return annotation.sequence, annotation.mod_dict()


def strip_mods(sequence: Union[str, ProFormaAnnotation]) -> str:
    """
    Strips all modifications from the given sequence, returning the unmodified sequence.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return annotation.sequence


def reverse(sequence: Union[str, ProFormaAnnotation], swap_terms: bool = False, include_plus: bool = False) -> str:
    """
    Reverses the sequence, while preserving the position of any modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param swap_terms: If True, the N- and C-terminal modifications will be swapped.
    :type swap_terms: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        annotation = sequence_to_annotation(sequence)

    else:
        annotation = sequence

    reversed_annotation = annotation.reverse(swap_terms=swap_terms)
    return reversed_annotation.serialize(include_plus)


def shuffle(sequence: Union[str, ProFormaAnnotation], seed: int = None, include_plus: bool = False) -> str:
    """
    Shuffles the sequence, while preserving the position of any modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param seed: Seed for the random number generator.
    :type seed: int

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        annotation = sequence_to_annotation(sequence)

    else:
        annotation = sequence

    shifted_annotation = annotation.shuffle(seed)
    return shifted_annotation.serialize(include_plus)


def shift(sequence: Union[str, ProFormaAnnotation], n: int, include_plus: bool = False) -> str:
    """
    Shifts the sequence to the left by a given number of positions, while preserving the position of any modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param n: The number of positions to shift the sequence to the left.
    :type n: int

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    shifted_annotation = annotation.shift(n)
    return shifted_annotation.serialize(include_plus)


def span_to_sequence(sequence: Union[str, ProFormaAnnotation], span: Span, include_plus: bool = False) -> str:
    """
    Extracts a subsequence from the input sequence based on the provided span.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param span: A tuple representing the span of the subsequence to be extracted.
    :type span: Tuple[int, int, int]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return annotation.slice(span[0], span[1]).serialize(include_plus)


def split(sequence: Union[str, ProFormaAnnotation], include_plus: bool = False) -> List[str]:
    """
    Splits sequence into a list of amino acids, preserving modifications.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        annotation = sequence_to_annotation(sequence)

    else:
        annotation = sequence

    return [annot.serialize(include_plus) for annot in annotation.split()]


def count_residues(sequence: Union[str, ProFormaAnnotation]) -> CounterType:
    """
    Counts the occurrences of each amino acid in the input sequence.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A Counter object containing the occurrences of each amino acid in the input sequence.
    :rtype: Counter

    .. code-block:: python

        >>> count_residues('PEPTIDE')
        Counter({'P': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1})

        >>> count_residues('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]')
        Counter({'[Acetyl]-P[phospho]': 1, 'E': 1, 'P[phospho]': 1, 'T': 1, 'I': 1, 'D': 1, 'E-[Amide]': 1})

        >>> count_residues('<13C>PE[3.14]PTIDE')
        Counter({'<13C>P': 2, '<13C>E[3.14]': 1, '<13C>T': 1, '<13C>I': 1, '<13C>D': 1, '<13C>E': 1})

        >>> count_residues('<[3.14]@C>PEPTIDE')
        Counter({'P': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1})

        >>> count_residues('<[3.14]@E>PEPTIDE')
        Counter({'P': 2, 'E[3.14]': 2, 'T': 1, 'I': 1, 'D': 1})

    """
    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    new_annotation = annotation.condense_static_mods(inplace=False)
    return new_annotation.count_residues()


def is_subsequence(subsequence: Union[str, ProFormaAnnotation], sequence: Union[str, ProFormaAnnotation],
                   order: bool = True) -> bool:
    """
    Checks if the input subsequence is a subsequence of the input sequence. If order is True, the subsequence must be in
    the same order as in the sequence. If order is False, the subsequence can be in any order.

    :param subsequence: The sequence or ProFormaAnnotation object, representing the subsequence.
    :type subsequence: Union[str, ProFormaAnnotation]
    :param sequence: The sequence or ProFormaAnnotation object, representing the sequence.
    :type sequence: Union[str, ProFormaAnnotation]
    :param order: If True, the subsequence must be in the same order as in the sequence.
    :type order: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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

    subsequence_counts = count_residues(subsequence)
    sequence_counts = count_residues(sequence)
    return all(subsequence_counts[aa] <= sequence_counts[aa] for aa in subsequence_counts)


def _sort_mods(mods: ModDict, sort_function: Optional[Callable[[str], str]] = None) -> None:
    """
    Sorts the modifications in the input dictionary using the provided sort function.

    :param mods: The modifications to be sorted.
    :type mods: Dict[int, List[ModValue]]
    :param sort_function: The sort function to be used. Defaults to identity function.
    :type sort_function: Callable[[str], str]

    :return: None

    .. code-block:: python

        >>> mod_dict = {1: ['phospho', 1], 2: ['phospho']}
        >>> _sort_mods(mod_dict)
        >>> mod_dict
        {1: [1, 'phospho'], 2: ['phospho']}

    """

    if sort_function is None:
        sort_function = lambda x: str(x)

    for k in mods:
        mods[k].sort(key=sort_function)


def sort(sequence: Union[str, ProFormaAnnotation], include_plus: bool = False) -> str:
    """
    Sorts the input sequence using the provided sort function. Terminal sequence are kept in place.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        annotation = sequence_to_annotation(sequence)

    else:
        annotation = sequence

    sorted_annotation = annotation.sort_residues()
    return sorted_annotation.serialize(include_plus)


def find_subsequence_indices(sequence: Union[str, ProFormaAnnotation],
                             subsequence: Union[str, ProFormaAnnotation],
                             ignore_mods: bool = False) -> List[int]:
    """
    Retrieves all starting indexes of a given subsequence within a sequence.

    :param sequence: The sequence or ProFormaAnnotation objectto look for the 'subsequence' in.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequence: The sequence or ProFormaAnnotation object to be found within 'sequence'.
    :type subsequence: Union[str, ProFormaAnnotation]
    :param ignore_mods: Whether to ignore modifications.
    :type ignore_mods: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        sequence = sequence_to_annotation(sequence)

    if isinstance(subsequence, str):
        subsequence = sequence_to_annotation(subsequence)

    if not sequence.has_sequence() or sequence.sequence == '':
        return []

    if not subsequence.has_sequence() or subsequence.sequence == '':
        return []

    if ignore_mods:
        sequence = sequence.strip()
        subsequence = subsequence.strip()

    return subsequence.find_indices(sequence)


def coverage(sequence: Union[str, ProFormaAnnotation],
             subsequences: List[Union[str, ProFormaAnnotation]],
             accumulate: bool = False,
             ignore_mods: bool = False) -> List[int]:
    """
    Calculate the sequence coverage given a list of subsequecnes.

    The coverage is represented as a binary list where each position in the protein sequence is marked as 1 if it
    is covered by at least one peptide and 0 otherwise.

    :param sequence: The sequence or ProFormaAnnotation object to be covered.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequences: The sequence's or ProFormaAnnotation object's subsequences to be used for coverage.
    :type subsequences: List[Union[str, ProFormaAnnotation]]
    :param accumulate: If True, the coverage array will be accumulated, i.e. if a position is covered by more than
                        one subsequence, it will be marked as the sum of the number of subsequences covering it.
                        If False, the coverage array will be binary, i.e. if a position is covered by more than one
                        subsequence, it will be marked as 1.
    :type accumulate: bool
    :param ignore_mods: Whether to ignore modifications when calcualting the coverage.
    :type ignore_mods: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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
        sequence = sequence_to_annotation(sequence)

    cov_arr = [0] * sequence_length(sequence)
    for subsequence in subsequences:

        if isinstance(subsequence, str):
            subsequence = sequence_to_annotation(subsequence)

        peptide_indexes = find_subsequence_indices(sequence, subsequence, ignore_mods=ignore_mods)
        for peptide_index in peptide_indexes:
            if accumulate:
                cov_arr[peptide_index:peptide_index + sequence_length(subsequence)] = \
                    [x + 1 for x in cov_arr[peptide_index:peptide_index + sequence_length(subsequence)]]
            else:
                cov_arr[peptide_index:peptide_index + sequence_length(subsequence)] = [1] * sequence_length(subsequence)

    return cov_arr


def percent_coverage(sequence: Union[str, ProFormaAnnotation],
                     subsequences: List[Union[str, ProFormaAnnotation]],
                     ignore_mods: bool = False) -> float:
    """
    Calculates the coverage given a list of subsequences.

    :param sequence: The sequence or ProFormaAnnotation object to be covered.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequences: The sequence's or ProFormaAnnotation object's subsequences to be used for coverage.
    :type subsequences: List[Union[str, ProFormaAnnotation]]
    :param ignore_mods: Whether to ignore modifications when calcualting the coverage.
    :type ignore_mods: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

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


def convert_diann_sequence(sequence: str) -> str:
    """
    Converts a IP2-Like sequence to a proforma2.0 compatible sequence.

    :param sequence: The sequence to be converted.
    :type sequence: str

    :return: Proforma2.0 compatable sequence.
    :rtype: str

    .. code-block:: python

        >>> convert_diann_sequence('_YMGTLRGC[Carbamidomethyl]LLRLYHD_')
        'YMGTLRGC[Carbamidomethyl]LLRLYHD'

        >>> convert_diann_sequence('_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD_')
        '[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD'

        >>> convert_diann_sequence('_[Acytel]YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]_[Methyl]')
        '[Acytel]-YMGTLRGC[Carbamidomethyl]LLRLYHD[1.0]-[Methyl]'

    """

    # Check if sequence starts and ends with underscores and remove them
    if sequence.startswith('_'):
        sequence = sequence[1:]

        # Check for a modification at the start of the sequence
        if re.match(r'^\[[^\]]+\]', sequence):
            sequence = re.sub(r'^\[([^\]]+)\]', r'[\1]-', sequence)

    if sequence.endswith('_'):
        sequence = sequence[:-1]

    elif re.search(r'_\[[^\]]+\]$', sequence):
        sequence = re.sub(r'_\[([^\]]+)\]$', r'-[\1]', sequence)

    return sequence


def convert_casanovo_sequence(sequence: str) -> str:
    """
    Converts a sequence with modifications to a proforma2.0 compatible sequence.
    :param sequence: The sequence to be converted.
    :type sequence: str
    :return: Proforma2.0 compatable sequence.
    :rtype: str

    .. code-block:: python

        >>> convert_casanovo_sequence('+43.006P+100EPTIDE')
        '[+43.006]-P[+100]EPTIDE'

    """
    new_sequence = []
    in_mod = False  # Tracks if we are within a modification
    is_nterm = False  # Tracks if the current modification is at the N-terminus

    for _, char in enumerate(sequence):
        if char in {'+', '-'}:
            # Check if it's at the start (N-terminal)
            is_nterm = len(new_sequence) == 0

            # Start a new modification block
            new_sequence.append('[')
            new_sequence.append(char)
            in_mod = True
        elif in_mod and char.isalpha():
            # End the modification block
            new_sequence.append(']')

            if is_nterm:
                # Add a dash if it's an N-terminal modification
                new_sequence.append('-')
                is_nterm = False

            # Add the current character and close modification
            in_mod = False
            new_sequence.append(char)
        else:
            # Add regular characters
            new_sequence.append(char)

    # Close any unclosed modification at the end of the sequence
    if in_mod:
        new_sequence.append(']')

    return ''.join(new_sequence)


def is_sequence_valid(sequence: Union[str, ProFormaAnnotation]) -> bool:
    """
    Checks if the input sequence is a valid ProForma sequence.

    :param sequence: The sequence or ProFormaAnnotation object to be validated.
    :type sequence: Union[str, ProFormaAnnotation]

    :return: True if the sequence is a valid ProForma sequence, False otherwise.
    :rtype: bool

    """

    if isinstance(sequence, str):
        try:
            _ = sequence_to_annotation(sequence)
        except Exception as err:
            return False
    return True

def count_aa(sequence: Union[str, ProFormaAnnotation]) -> Dict[str, int]:
    """
    Converts a sequence to a feature vector.
    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    aa_counts = {aa:0 for aa in ORDERED_AMINO_ACIDS}
    for aa in annotation.sequence:
        aa_counts[aa] += 1

    return aa_counts
