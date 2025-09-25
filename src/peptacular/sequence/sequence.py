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

from typing import Any, Callable, Iterable, Mapping, Sequence

from .util import get_annotation_input
from ..constants import ORDERED_AMINO_ACIDS, ModType, ModTypeLiteral
from ..proforma.annotation import (
    ProFormaAnnotation,
)


def sequence_length(sequence: str | ProFormaAnnotation) -> int:
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

    return len(get_annotation_input(sequence, copy=False))


def is_ambiguous(sequence: str | ProFormaAnnotation) -> bool:
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

    return get_annotation_input(sequence, copy=False).contains_sequence_ambiguity()


def is_modified(sequence: str | ProFormaAnnotation) -> bool:
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

    return get_annotation_input(sequence, copy=False).has_mods()


def get_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | ModTypeLiteral | None = None,
) -> dict[str, Any]:
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
        >>> get_mods('PEP[Phospho]T[1]IDE[-3.14]')['internal']
        {2: ('Phospho',), 3: (1,), 6: (-3.14,)}

        >>> get_mods('PEP[Phospho][1.0]TIDE')['internal']
        {2: ('Phospho', 1.0)}

        # N-terminal modifications are mapped to the index 'nterm'
        >>> get_mods('[Acetyl]-PEPTIDE')['nterm']
        ('Acetyl',)

        # C-terminal modifications are mapped to the index 'cterm'
        >>> get_mods('PEPTIDE-[Amide]')['cterm']
        ('Amide',)

        # Isotopic modifications are mapped to the index 'isotope'
        >>> get_mods('<13C>PEPTIDE')['isotope']
        ('13C',)

        # Static modifications are mapped to the index 'static'
        >>> get_mods('<[+1.234]@P>PEPTIDE')['static']
        ('[+1.234]@P',)

        # Labile modifications are mapped to the index 'labile'
        >>> get_mods('{Glycan:Hex}PEPTIDE')['labile']
        ('Glycan:Hex',)

        # Unknown modifications are mapped to the index 'unknown'
        >>> get_mods('[Phospho]^3?PEPTIDE')['unknown']
        ('Phospho', 'Phospho', 'Phospho')

        # Intervals are mapped to the index 'interval'
        >>> get_mods('PEP(TI)[Phospho]DE')['interval']
        (ModInterval(start=3, end=5, ambiguous=False, mods=('Phospho',)),)

        # Charge state is mapped to the index 'charge'
        >>> get_mods('PEPTIDE/+2')['charge']
        2

        # Charge adducts are mapped to the index 'charge_adducts'
        >>> get_mods('PEPTIDE/+2[+2Na+,-H+]')['charge_adducts']
        ('+2Na+,-H+',)

    """

    return get_annotation_input(sequence, copy=True).get_mods(mods)


def add_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any],
    append: bool = True,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
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

        >>> from peptacular import Mod

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
        >>> add_mods('PEPTIDE', {'interval': (3, 5, False, 'Phospho')})
        'PEP(TI)[Phospho]DE'

        # Can also add charge state
        >>> add_mods('PEPTIDE', {'charge': 2})
        'PEPTIDE/2'

        # Can also add charge adducts
        >>> add_mods('PEPTIDE', {'charge': 2, 'charge_adducts': '+2Na+,-H+'})
        'PEPTIDE/2[+2Na+,-H+]'

    """

    annotation = get_annotation_input(sequence, copy=True)
    annotation.add_mods(mods, inplace=True, append=append)
    return annotation.serialize(include_plus=include_plus, precision=precision)


def set_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any] | None,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    return (
        get_annotation_input(sequence, copy=True)
        .set_mods(mods, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def append_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any],
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    return (
        get_annotation_input(sequence, copy=True)
        .append_mods(mods, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def extend_mods(
    sequence: str | ProFormaAnnotation,
    mods: Mapping[ModType | ModTypeLiteral | int, Any],
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    return (
        get_annotation_input(sequence, copy=True)
        .extend_mods(mods, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def condense_static_mods(
    sequence: str | ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
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

    return (
        get_annotation_input(sequence=sequence, copy=True)
        .condense_static_mods(inplace=False)
        .serialize(include_plus=include_plus, precision=precision)
    )


def pop_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
    include_plus: bool = False,
    precision: int | None = None,
) -> tuple[str, dict[str, Any]]:
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
        >>> pop_mods('PEP[phospho]TIDE')[0]
        'PEPTIDE'
        >>> pop_mods('PEP[phospho]TIDE')[1]['internal']
        {2: ('phospho',)}

        # can specify which modifications to pop
        >>> pop_mods('PEP[phospho]TIDE-[+100]', mods=['internal'])
        ('PEPTIDE-[100]', {'internal': {2: ('phospho',)}})

    """
    annotation = get_annotation_input(sequence=sequence, copy=True)
    mod_dict = annotation.pop_mods(mod_types=mods)  # only include keys that have mods
    return (
        annotation.serialize(include_plus=include_plus, precision=precision),
        mod_dict,
    )


def strip_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
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

    annotation = get_annotation_input(sequence=sequence, copy=True)

    return annotation.remove_mods(mods=mods, inplace=True).serialize(
        include_plus=include_plus, precision=precision
    )


def filter_mods(
    sequence: str | ProFormaAnnotation,
    mods: ModType | Iterable[ModType] | None = None,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """
    Keeps only the specified modifications in the sequence, removing all others.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param mods: The modifications to keep. If None, all modifications will be kept.
    :type mods: Optional[Union[str, List[str]]]
    :param include_plus: If True, the modifications will be serialized with a '+' sign for positive values.
    :type include_plus: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence with only the specified modifications kept.
    :rtype: str

    .. code-block:: python

        # Keeps only internal modifications:
        >>> filter_mods('PEP[phospho]TIDE', mods='internal')
        'PEP[phospho]TIDE'

        # Keeps only N and C terminal modifications:
        >>> filter_mods('[Acetyl]-PEPTIDE[1.234]-[Amide]', mods=['nterm', 'cterm'])
        '[Acetyl]-PEPTIDE-[Amide]'

        # Keeps only labile modifications:
        >>> filter_mods('{1.0}[Acetyl]-PEPTIDE[1.234]-[Amide]', mods='labile')
        '{1.0}PEPTIDE'

        # Keeps only isotope notations:
        >>> filter_mods('<C13>[Acetyl]-PEPTIDE[1.234]-[Amide]', mods='isotope')
        '<C13>PEPTIDE'

        # Using a sequence without modifications will return the same sequence:
        >>> filter_mods('PEPTIDE')
        'PEPTIDE'

    """
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .filter_mods(mods=mods, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def reverse(
    sequence: str | ProFormaAnnotation,
    swap_terms: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
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
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .reverse(inplace=True, swap_terms=swap_terms)
        .serialize(include_plus=include_plus, precision=precision)
    )


def shuffle(
    sequence: str | ProFormaAnnotation,
    seed: int | None = None,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
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
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .shuffle(seed=seed, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def shift(
    sequence: str | ProFormaAnnotation, n: int, include_plus: bool = False
) -> str:
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
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .shift(n=n, inplace=True)
        .serialize(include_plus=include_plus)
    )


def span_to_sequence(
    sequence: str | ProFormaAnnotation,
    span: tuple[int, int, int],
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
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

    """

    return (
        get_annotation_input(sequence=sequence, copy=True)
        .slice(span[0], span[1], inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


def split(
    sequence: str | ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> list[str]:
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
    return [
        a.serialize(include_plus=include_plus, precision=precision)
        for a in get_annotation_input(sequence=sequence, copy=True).split()
    ]


def count_residues(sequence: str | ProFormaAnnotation) -> dict[str, int]:
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
        {'P': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1}

        >>> count_residues('[Acetyl]-P[phospho]EP[phospho]TIDE-[Amide]')
        {'[Acetyl]-P[phospho]': 1, 'E': 1, 'P[phospho]': 1, 'T': 1, 'I': 1, 'D': 1, 'E-[Amide]': 1}

        >>> count_residues('<13C>PE[3.14]PTIDE')
        {'<13C>P': 2, '<13C>E[3.14]': 1, '<13C>T': 1, '<13C>I': 1, '<13C>D': 1, '<13C>E': 1}

        >>> count_residues('<[3.14]@C>PEPTIDE')
        {'P': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1}

        >>> count_residues('<[3.14]@E>PEPTIDE')
        {'P': 2, 'E[3.14]': 2, 'T': 1, 'I': 1, 'D': 1}

    """
    return (
        get_annotation_input(sequence, copy=False)
        .condense_static_mods(inplace=True)
        .count_residues()
    )


def percent_residues(
    sequence: str | ProFormaAnnotation, precision: int | None = None
) -> dict[str, float]:
    """
    Calculates the percentage of each amino acid in the input sequence.
    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid
    :return: A dictionary containing the percentage of each amino acid in the input sequence.
    :rtype: Dict[str, float]
    .. code-block:: python

        >>> percent_residues('PEPTIDE', precision=2)
        {'P': 28.57, 'E': 28.57, 'T': 14.29, 'I': 14.29, 'D': 14.29}

    """
    return (
        get_annotation_input(sequence, copy=False)
        .condense_static_mods(inplace=True)
        .percent_residues(precision=precision)
    )


def is_subsequence(
    subsequence: str | ProFormaAnnotation,
    sequence: str | ProFormaAnnotation,
    order: bool = True,
) -> bool:
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
    return all(
        subsequence_counts[aa] <= sequence_counts[aa] for aa in subsequence_counts
    )


def sort(
    sequence: str | ProFormaAnnotation,
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
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
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .sort(inplace=True, key=key, reverse=reverse)
        .serialize(include_plus=include_plus, precision=precision)
    )


def find_subsequence_indices(
    sequence: str | ProFormaAnnotation,
    subsequence: str | ProFormaAnnotation,
    ignore_mods: bool = False,
) -> list[int]:
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
    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequence_annot = get_annotation_input(subsequence, copy=False)
    return subsequence_annot.find_indices(other=sequence_annot, ignore_mods=ignore_mods)


def coverage(
    sequence: str | ProFormaAnnotation,
    subsequences: Iterable[str | ProFormaAnnotation],
    accumulate: bool = False,
    ignore_mods: bool = False,
    ignore_ambiguity: bool = False,
) -> list[int]:
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
    :param ignore_ambiguity: Whether to ignore ambiguity codes when calculating the coverage.
    :type ignore_ambiguity: bool

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

        # By default ambiguity does not add to coverage
        >>> coverage("PEPTIDE", ["P(?EP)"])
        [1, 0, 0, 0, 0, 0, 0]

    """

    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequences_annot = [
        get_annotation_input(subseq, copy=False) for subseq in subsequences
    ]

    return sequence_annot.coverage(
        subsequences_annot,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )


def percent_coverage(
    sequence: str | ProFormaAnnotation,
    subsequences: list[str | ProFormaAnnotation],
    ignore_mods: bool = False,
    accumulate: bool = False,
    ignore_ambiguity: bool = False,
) -> float:
    """
    Calculates the coverage given a list of subsequences.

    :param sequence: The sequence or ProFormaAnnotation object to be covered.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequences: The sequence's or ProFormaAnnotation object's subsequences to be used for coverage.
    :type subsequences: List[Union[str, ProFormaAnnotation]]
    :param ignore_mods: Whether to ignore modifications when calcualting the coverage.
    :type ignore_mods: bool
    :param accumulate: If True, the coverage array will be accumulated
    :param accumulate: bool
    :param ignore_ambiguity: Whether to ignore ambiguity codes when calculating the coverage.
    :type ignore_ambiguity: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The percent coverage.
    :rtype: float

    .. code-block:: python

        >>> round(percent_coverage("PEPTIDE", ["PEP"]), 3)
        0.429

        >>> round(percent_coverage("PEPTIDE", ["PEP", "EPT"]), 3)
        0.571

        >>> round(percent_coverage("PEPTIDE", ["PEPTIDE", "PEPTIDE"], accumulate=True), 3)
        2.0

    """

    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequences_annot = [
        get_annotation_input(subseq, copy=False) for subseq in subsequences
    ]

    return sequence_annot.percent_coverage(
        subsequences_annot,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )


def modification_coverage(
    sequence: str | ProFormaAnnotation,
    subsequences: list[str | ProFormaAnnotation],
    accumulate: bool = False,
) -> dict[int, int]:
    """
    Calculate the modification coverage given a list of subsequences.

    This function identifies which modifications in the main sequence are covered by
    subsequences. It returns a dictionary where each key is a modification position/type
    (matching the format from get_mods) and each value is the number of subsequences
    that cover that modification.

    :param sequence: The sequence or ProFormaAnnotation object representing the protein.
    :type sequence: Union[str, ProFormaAnnotation]
    :param subsequences: The sequence's or ProFormaAnnotation object's subsequences (peptides) to be used for coverage.
    :type subsequences: List[Union[str, ProFormaAnnotation]]
    :param accumulate: If True, the count will accumulate for each subsequence covering the modification.
                       If False, the count will be binary (1 if covered, 0 if not).
    :type accumulate: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: Dictionary mapping modification positions/types to coverage counts.
    :rtype: Dict[Union[int, str], int]

    .. code-block:: python

        >>> modification_coverage("PEPTIDE[Phospho]", ["TIDE"])
        {6: 0}

        >>> modification_coverage("PEPTIDE[Phospho]", ["TIDE[Phospho]"])
        {6: 1}

        >>> modification_coverage("PEP[Phospho]TIDE[Methyl]", ["PEP[Phospho]", "TIDE[Methyl]"])
        {2: 1, 6: 1}

        >>> modification_coverage("PEP[Phospho]TIDE[Methyl]", ["PEP[Phospho]", "TIDE[Methyl]", "PEP[Phospho]"], accumulate=True)
        {2: 2, 6: 1}

    """
    sequence_annot = get_annotation_input(sequence, copy=False)
    subsequence_annots = [
        get_annotation_input(subseq, copy=False) for subseq in subsequences
    ]

    return sequence_annot.modification_coverage(
        annotations=subsequence_annots, accumulate=accumulate
    )


def count_aa(sequence: str | ProFormaAnnotation) -> dict[str, int]:
    """
    Converts a sequence to a feature vector.
    """

    annotation = get_annotation_input(sequence, copy=False)

    aa_counts = {aa: 0 for aa in ORDERED_AMINO_ACIDS}
    for aa in annotation.sequence:
        aa_counts[aa] += 1

    return aa_counts


def annotate_ambiguity(
    sequence: str | ProFormaAnnotation,
    forward_coverage: list[int],
    reverse_coverage: list[int],
    mass_shift: Any | None = None,
    add_mods_to_intervals: bool = False,
    sort_mods: bool = True,
    include_plus: bool = False,
    condense_to_xnotation: bool = False,
    precision: int | None = None,
) -> str:
    """
    Given a peptide sequence and coverage information for forward and reverse fragment ions,
    annotate the sequence with ambiguity intervals and optional mass shift.

    This function identifies regions in the sequence where there is insufficient fragment ion
    coverage and marks them as ambiguous using ProForma notation with parentheses.
    If a mass shift is provided , it will be added to the appropriate location.

    :param sequence: The peptide sequence or ProFormaAnnotation object to annotate.
    :type sequence: Union[str, ProFormaAnnotation]
    :param forward_coverage: Binary list indicating which positions have forward ion coverage (1) or not (0).
    :type forward_coverage: List[int]
    :param reverse_coverage: Binary list indicating which positions have reverse ion coverage (1) or not (0).
    :type reverse_coverage: List[int]
    :param mass_shift: An optional mass shift to be added to the sequence at the appropriate position.
    :type mass_shift: Optional[Any]

    :raises ValueError: If the annotation already contains intervals or if coverage lengths don't match sequence length.

    :return: The annotated sequence with ambiguity regions marked using ProForma notation.
    :rtype: str

    .. code-block:: python

        # Add ambiguity intervals based on fragment ion coverage
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,1,0,0,0], [0,0,0,0,0,1,0])
        '(?PE)PTI(?DE)'

        # With a phosphorylation mass shift
        >>> annotate_ambiguity('PEPTIDE', [1,1,1,0,0,0,0], [0,0,0,0,1,1,1], 79.966)
        'PEPT[79.966]IDE'

        # Handling existing modifications
        >>> annotate_ambiguity('P[10]EPTIDE', [1,1,1,0,0,0,0], [0,0,0,0,0,1,1])
        'P[10]EP(?TI)DE'

        # When mass shift can't be localized to a specific residue
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,0,0,0,0], [0,0,0,0,0,1,0], 120)
        '(?PE)P(?TI)[120](?DE)'

        # When mass shift is completely unlocalized, it becomes a labile modification
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,1,1,0,0], [0,0,1,1,1,1,0], 120)
        '{120}(?PE)PTI(?DE)'

        (?SSGS)IA(?SS)(?YVQ)()[37.959]W(?YQQRPGSA)(?PT)TVIYEDDER(?PS)(?GV)(?PDR)
        Seq: SSGSIASSYVQWYQQRPGSAPTTVIYEDDERPSGVPDR
        For: 00011101001000000000000000000000000000
        Rev: 00000000000110000000101111111111010100
        >>> for_ions = list(map(int, '00011101001000000000000000000000000000'))
        >>> rev_iosn = list(map(int, '00000000000110000000101111111111010100'))
        >>> annotate_ambiguity('SSGSIASSYVQWYQQRPGSAPTTVIYEDDERPSGVPDR', for_ions, rev_iosn, 120)
        '(?SSGS)IA(?SS)(?YVQ)W[120](?YQQRPGSA)(?PT)TVIYEDDER(?PS)(?GV)(?PDR)'
    """
    annot = get_annotation_input(sequence=sequence, copy=True).annotate_ambiguity(
        forward_coverage=forward_coverage,
        reverse_coverage=reverse_coverage,
        mass_shift=mass_shift,
        add_mods_to_intervals=add_mods_to_intervals,
        sort_mods=sort_mods,
        inplace=True,
    )

    if condense_to_xnotation:
        annot.condense_ambiguity_to_xnotation(inplace=True)

    return annot.serialize(include_plus=include_plus, precision=precision)


def join(
    annotations: Sequence[ProFormaAnnotation | str],
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """
    Join a list of ProFormaAnnotation objects into a single annotation.

    :param annotations: The list of ProFormaAnnotation objects to join.
    :type annotations: List[ProFormaAnnotation]

    :raises ValueError: If the list is empty.

    :return: The joined ProFormaAnnotation object.
    :rtype: ProFormaAnnotation
    """
    if not annotations:
        raise ValueError("Cannot join an empty list of annotations.")

    annotations = [get_annotation_input(a, copy=False) for a in annotations]
    return ProFormaAnnotation.join(annotations).serialize(
        include_plus=include_plus, precision=precision
    )
