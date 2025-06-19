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

from typing import Counter as CounterType, Optional, Any
from typing import Dict, List, Tuple, Callable, Union

from ..errors import ProFormaFormatError
from ..constants import ORDERED_AMINO_ACIDS
from ..proforma.proforma_parser import (
    parse,
    ProFormaAnnotation,
    serialize,
    MultiProFormaAnnotation,
)
from ..spans import Span
from ..proforma.input_convert import (
    ModDict,
    fix_list_of_mods,
    fix_intervals_input,
)


def get_annotation_input(
    sequence: Union[str, ProFormaAnnotation], copy: bool = True
) -> ProFormaAnnotation:
    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    elif isinstance(sequence, ProFormaAnnotation):
        if copy:
            annotation = sequence.copy()
        else:
            annotation = sequence
    else:
        raise TypeError(f"Expected str or ProFormaAnnotation, got {type(sequence)}")

    return annotation


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

    return len(get_annotation_input(sequence, copy=False))


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

    return get_annotation_input(sequence, copy=False).contains_sequence_ambiguity()


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

    return get_annotation_input(sequence, copy=False).has_mods()


def get_mods(
    sequence: Union[str, ProFormaAnnotation],
    mods: Optional[Union[str, List[str]]] = None,
) -> ModDict:
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
        {'internal': {2: [Mod('Phospho', 1)], 3: [Mod(1, 1)], 6: [Mod(-3.14, 1)]}}

        >>> get_mods('PEP[Phospho][1.0]TIDE')
        {'internal': {2: [Mod('Phospho', 1), Mod(1.0, 1)]}}

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

    return get_annotation_input(sequence, copy=True).mod_dict(mods)


def add_mods(
    sequence: Union[str, ProFormaAnnotation],
    mods: Dict,
    append: bool = True,
    include_plus: bool = False,
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

    for k in mods:
        if k == "charge":
            continue
        if k == "intervals":
            mods[k] = fix_intervals_input(mods[k])
        else:
            mods[k] = fix_list_of_mods(mods[k])

    return (
        get_annotation_input(sequence, copy=True)
        .add_mod_dict(mods, append=append)
        .serialize(include_plus=include_plus)
    )


def condense_static_mods(
    sequence: Union[str, ProFormaAnnotation], include_plus: bool = False
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
        get_annotation_input(sequence, copy=True)
        .condense_static_mods(inplace=False)
        .serialize(include_plus)
    )


def pop_mods(
    sequence: Union[str, ProFormaAnnotation],
    mods: Optional[Union[str, List[str]]] = None,
    include_plus: bool = False,
) -> Tuple[str, ModDict]:
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
        ('PEPTIDE', {'internal': {2: [Mod('phospho', 1)]}})

        # can specify which modifications to pop
        >>> pop_mods('PEP[phospho]TIDE-[+100]', mods=['internal'])
        ('PEPTIDE-[100]', {'internal': {2: [Mod('phospho', 1)]}})

    """
    annotation = get_annotation_input(sequence, copy=True)
    mod_dict = annotation.pop_mods(mods)
    return annotation.serialize(include_plus=include_plus), mod_dict


def strip_mods(
    sequence: Union[str, ProFormaAnnotation],
    mods: Optional[Union[str, List[str]]] = None,
    include_plus: bool = False,
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

    return pop_mods(sequence, mods, include_plus)[0]


def filter_mods(
    sequence: Union[str, ProFormaAnnotation],
    mods: Optional[Union[str, List[str]]] = None,
    include_plus: bool = False,
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
        get_annotation_input(sequence, copy=True)
        .filter_mods(mods, inplace=True)
        .serialize(include_plus=include_plus)
    )


def remove_mods(
    sequence: Union[str, ProFormaAnnotation],
    mods: Optional[Union[str, List[str]]] = None,
    include_plus: bool = False,
) -> str:
    """
    Removes the specified modifications from the sequence, returning the modified sequence.

    :param sequence: The sequence or ProFormaAnnotation object.
    :type sequence: Union[str, ProFormaAnnotation]
    :param mods: The modifications to remove. If None, all modifications will be removed.
    :type mods: Optional[Union[str, List[str]]]
    :param include_plus: If True, the modifications will be serialized with a '+' sign for positive values.
    :type include_plus: bool

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence with the specified modifications removed.
    :rtype: str

    .. code-block:: python

        # Removes internal modifications:
        >>> remove_mods('PEP[phospho]TIDE', mods='internal')
        'PEPTIDE'

        # Removes N and C terminal modifications:
        >>> remove_mods('[Acetyl]-PEPTIDE[1.234]-[Amide]', mods=['nterm', 'cterm'])
        'PEPTIDE[1.234]'

        # Removes labile modifications:
        >>> remove_mods('{1.0}[Acetyl]-PEPTIDE[1.234]-[Amide]', mods='labile')
        '[Acetyl]-PEPTIDE[1.234]-[Amide]'

        # Removes isotope notations:
        >>> remove_mods('{1.0}<C13>[Acetyl]-PEPTIDE[1.234]-[Amide]', mods='isotope')
        '{1.0}[Acetyl]-PEPTIDE[1.234]-[Amide]'

    """
    return (
        get_annotation_input(sequence, copy=True)
        .remove_mods(mods, inplace=True)
        .serialize(include_plus=include_plus)
    )


def reverse(
    sequence: Union[str, ProFormaAnnotation],
    swap_terms: bool = False,
    include_plus: bool = False,
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
        get_annotation_input(sequence, copy=True)
        .reverse(swap_terms=swap_terms)
        .serialize(include_plus)
    )


def shuffle(
    sequence: Union[str, ProFormaAnnotation],
    seed: int = None,
    include_plus: bool = False,
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
        get_annotation_input(sequence, copy=True).shuffle(seed).serialize(include_plus)
    )


def shift(
    sequence: Union[str, ProFormaAnnotation], n: int, include_plus: bool = False
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
    return get_annotation_input(sequence, copy=True).shift(n).serialize(include_plus)


def span_to_sequence(
    sequence: Union[str, ProFormaAnnotation], span: Span, include_plus: bool = False
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

        >>> span_to_sequence('(PEPT)IDE', (1, 6, 0))
        '(EPT)ID'

    """
    return (
        get_annotation_input(sequence, copy=True)
        .slice(span[0], span[1], inplace=False)
        .serialize(include_plus)
    )


def split(
    sequence: Union[str, ProFormaAnnotation], include_plus: bool = False
) -> List[str]:
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
        a.serialize(include_plus)
        for a in get_annotation_input(sequence, copy=True).split()
    ]


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
    return (
        get_annotation_input(sequence, copy=False)
        .condense_static_mods(inplace=True)
        .count_residues()
    )

def percent_residues(
    sequence: Union[str, ProFormaAnnotation], precision: Optional[int] = None) -> Dict[str, float]:
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
    counts = count_residues(sequence)
    total = sum(counts.values())
    if total == 0:
        return {}
    d = {aa: (count / total) * 100 for aa, count in counts.items()}

    if precision is not None:
        d = {aa: round(value, precision) for aa, value in d.items()}

    return d


def is_subsequence(
    subsequence: Union[str, ProFormaAnnotation],
    sequence: Union[str, ProFormaAnnotation],
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


def _sort_mods(
    mods: ModDict, sort_function: Optional[Callable[[str], str]] = None
) -> None:
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
    return (
        get_annotation_input(sequence, copy=True)
        .sort(inplace=True)
        .serialize(include_plus)
    )


def find_subsequence_indices(
    sequence: Union[str, ProFormaAnnotation],
    subsequence: Union[str, ProFormaAnnotation],
    ignore_mods: bool = False,
) -> List[int]:
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
    sequence: Union[str, ProFormaAnnotation],
    subsequences: List[Union[str, ProFormaAnnotation]],
    accumulate: bool = False,
    ignore_mods: bool = False,
    ignore_ambiguity: bool = False,
) -> List[int]:
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

    cov_arr = [0] * len(sequence_annot)

    for subsequence in subsequences:
        subsequencue_annot = get_annotation_input(subsequence, copy=False)

        ambiguity_intervals = subsequencue_annot.pop_intervals()
        ambiguity_intervals = [i for i in ambiguity_intervals if i.ambiguous]

        subsequence_indicies = subsequencue_annot.find_indices(
            other=sequence_annot, ignore_mods=ignore_mods
        )

        peptide_cov = [1] * len(subsequencue_annot)
        if ignore_ambiguity is False:
            for interval in ambiguity_intervals:
                for i in range(interval.start, interval.end):
                    peptide_cov[i] = 0

        for subsequence_index in subsequence_indicies:
            start = subsequence_index

            for i, cov in enumerate(peptide_cov):
                if accumulate:
                    cov_arr[start + i] += cov
                else:
                    cov_arr[start + i] = cov

    return cov_arr


def percent_coverage(
    sequence: Union[str, ProFormaAnnotation],
    subsequences: List[Union[str, ProFormaAnnotation]],
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

    cov_arr = coverage(
        sequence,
        subsequences,
        accumulate=accumulate,
        ignore_mods=ignore_mods,
        ignore_ambiguity=ignore_ambiguity,
    )

    if len(cov_arr) == 0:
        return 0

    return sum(cov_arr) / len(cov_arr)


def modification_coverage(
    sequence: Union[str, ProFormaAnnotation],
    subsequences: List[Union[str, ProFormaAnnotation]],
    accumulate: bool = False,
) -> Dict[Union[int, str], int]:
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

    # Get all modifications from the main sequence
    sequence_mods = sequence_annot.mod_dict(mods=["internal", "nterm", "cterm"])

    # Initialize coverage dictionary with zeroes for all modification sites
    coverage_dict = {pos: 0 for pos in sequence_mods.get("internal", {})}
    coverage_dict.update(
        {mod_type: 0 for mod_type in ["nterm", "cterm"] if mod_type in sequence_mods}
    )

    # Get the unmodified sequence
    unmodified_sequence = sequence_annot.strip(inplace=False)

    # Process each subsequence
    for subsequence in subsequences:

        subsequence_anot = get_annotation_input(subsequence, copy=False)

        # Get the unmodified subsequence
        unmodified_subsequence = subsequence_anot.strip(inplace=False)

        start_indices = subsequence_anot.find_indices(
            other=sequence_annot, ignore_mods=True
        )

        # Get the modifications of the subsequence
        subsequence_mods = subsequence_anot.mod_dict(
            mods=["internal", "nterm", "cterm"]
        )

        # For each occurrence, check if the subsequence contains modifications at the same positions
        for start_idx in start_indices:
            for mod_pos, mod_values in sequence_mods.get("internal", {}).items():
                # Calculate the relative position in the subsequence
                relative_pos = mod_pos - start_idx

                # Check if this position is within the subsequence
                if 0 <= relative_pos < len(unmodified_subsequence):
                    # Check if the subsequence has a modification at this relative position
                    if relative_pos in subsequence_mods.get("internal", {}):

                        # For each modification in the main sequence at this position
                        for mod_value in mod_values:
                            # Check if any modification in the subsequence matches
                            if any(
                                str(mod_value) == str(subseq_mod)
                                for subseq_mod in subsequence_mods.get("internal", {})[
                                    relative_pos
                                ]
                            ):
                                if accumulate:
                                    coverage_dict[mod_pos] += 1
                                else:
                                    coverage_dict[mod_pos] = 1
                                break

            # Handle terminal modifications
            for mod_type in ["nterm", "cterm"]:
                if mod_type in sequence_mods and mod_type in subsequence_mods:
                    # For terminal mods, check if start or end of subsequence aligns with terminal
                    terminal_match = (mod_type == "nterm" and start_idx == 0) or (
                        mod_type == "cterm"
                        and start_idx + len(unmodified_subsequence)
                        == len(unmodified_sequence)
                    )

                    if terminal_match:
                        # Check if any terminal mod in the subsequence matches
                        for mod_value in sequence_mods[mod_type]:
                            if any(
                                str(mod_value) == str(subseq_mod)
                                for subseq_mod in subsequence_mods[mod_type]
                            ):
                                if accumulate:
                                    coverage_dict[mod_type] += 1
                                else:
                                    coverage_dict[mod_type] = 1
                                break

    return coverage_dict


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
        except ValueError:
            return False
        except ProFormaFormatError:
            return False
    return True


def count_aa(sequence: Union[str, ProFormaAnnotation]) -> Dict[str, int]:
    """
    Converts a sequence to a feature vector.
    """

    annotation = get_annotation_input(sequence, copy=False)

    aa_counts = {aa: 0 for aa in ORDERED_AMINO_ACIDS}
    for aa in annotation.sequence:
        aa_counts[aa] += 1

    return aa_counts


def annotate_ambiguity(
    sequence: Union[str, ProFormaAnnotation],
    forward_coverage: List[int],
    reverse_coverage: List[int],
    mass_shift: Optional[Any] = None,
    include_plus: bool = False,
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
    return (
        get_annotation_input(sequence, copy=True)
        .annotate_ambiguity(
            forward_coverage=forward_coverage,
            reverse_coverage=reverse_coverage,
            mass_shift=mass_shift,
            inplace=True,
        )
        .serialize(include_plus)
    )
