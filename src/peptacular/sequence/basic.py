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

from typing import Any, Sequence, overload

from ..constants import ParrallelMethod, ParrallelMethodLiteral
from ..annotation import (
    ProFormaAnnotation,
)
from .parrallel import parallel_apply_internal
from .util import get_annotation_input
from ..amino_acids import AA_LOOKUP


ORDERED_AMINO_ACIDS = [aa.id for aa in AA_LOOKUP.ordered_amino_acids]


def _parse_single(s: str, validate: bool = False) -> ProFormaAnnotation:
    """Parse a ProForma string into a ProFormaAnnotation object."""
    return ProFormaAnnotation.parse(s, validate=validate)


@overload
def parse(
    s: str,
    validate: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> ProFormaAnnotation: ...


@overload
def parse(
    s: Sequence[str],
    validate: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> list[ProFormaAnnotation]: ...


def parse(
    s: str | Sequence[str],
    validate: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
    reuse_pool: bool = True,
) -> ProFormaAnnotation | list[ProFormaAnnotation]:
    if isinstance(s, Sequence) and not isinstance(s, str):
        return parallel_apply_internal(
            _parse_single,
            s,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            validate=validate,
            reuse_pool=reuse_pool,
        )
    else:
        return _parse_single(s, validate=validate)


def _serialize_single(
    sequence: str | ProFormaAnnotation,
) -> str:
    """Internal function for serializing a single sequence."""
    return get_annotation_input(sequence, copy=False).serialize()


@overload
def serialize(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def serialize(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def serialize(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Serialize a peptide sequence or list of sequences to ProForma string format.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The serialized ProForma string or list of strings.
    :rtype: str | list[str]

    .. code-block:: python

        # Single sequence
        >>> serialize('PEPTIDE')
        'PEPTIDE'

        # Sequence with named modifications
        >>> serialize('PEP[Phospho]TIDE')
        'PEP[Phospho]TIDE'

        # Sequence with float modifications
        >>> serialize('PEP[+79.966]TIDE')
        'PEP[+79.966]TIDE'

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _serialize_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _serialize_single(sequence)


def _sequence_length_single(sequence: str | ProFormaAnnotation) -> int:
    """Internal function for computing sequence length for a single sequence."""
    return len(get_annotation_input(sequence, copy=False))


@overload
def sequence_length(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> int: ...


@overload
def sequence_length(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[int]: ...


def sequence_length(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> int | list[int]:
    """
    Compute the length of the peptide sequence based on the unmodified sequence.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sequence length, or list of lengths for list input.
    :rtype: int | list[int]

    .. code-block:: python

        # The length of the unmodified sequence is the same as the length of the string
        >>> sequence_length("PEPTIDE")
        7

        # Modifications are not counted in the sequence length
        >>> sequence_length("[Acetyl]-PEP[+1.2345]TID[+3.14]E-[Amide]")
        7
        >>> sequence_length("<C13>PEPTIDE[+1.2345]")
        7
        >>> sequence_length("(?PE)PTIDE[+1.2345]")
        7


    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _sequence_length_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _sequence_length_single(sequence)


def _is_ambiguous_single(sequence: str | ProFormaAnnotation) -> bool:
    """Internal function for checking if a single sequence is ambiguous."""
    return get_annotation_input(sequence, copy=False).contains_sequence_ambiguity()


@overload
def is_ambiguous(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> bool: ...


@overload
def is_ambiguous(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[bool]: ...


def is_ambiguous(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> bool | list[bool]:
    """
    Check if the sequence contains ambiguous amino acids.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: True if the sequence contains ambiguous amino acids, False otherwise, or list of bools.
    :rtype: bool | list[bool]

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
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _is_ambiguous_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _is_ambiguous_single(sequence)


def _is_modified_single(sequence: str | ProFormaAnnotation) -> bool:
    """Internal function for checking if a single sequence is modified."""
    return get_annotation_input(sequence, copy=False).has_mods()


@overload
def is_modified(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> bool: ...


@overload
def is_modified(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[bool]: ...


def is_modified(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> bool | list[bool]:
    """
    Check if the sequence contains any modifications.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: True if the sequence contains any modifications, False otherwise, or list of bools.
    :rtype: bool | list[bool]

    .. code-block:: python

        # Any modifications will return True
        >>> is_modified("PEP[Phospho]TIDE")
        True

        # Unmodified sequences will return False
        >>> is_modified("PEPTIDE")
        False


    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _is_modified_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _is_modified_single(sequence)


def _count_residues_single(sequence: str | ProFormaAnnotation) -> dict[str, int]:
    """Internal function for counting residues in a single sequence."""
    return (
        get_annotation_input(sequence, copy=False)
        .condense_static_mods(inplace=True)
        .count_residues()
    )


@overload
def count_residues(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> dict[str, int]: ...


@overload
def count_residues(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[dict[str, int]]: ...


def count_residues(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> dict[str, int] | list[dict[str, int]]:
    """
    Counts the occurrences of each amino acid in the input sequence.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A Counter object or list of Counter objects containing the occurrences of each amino acid.
    :rtype: dict[str, int] | list[dict[str, int]]

    .. code-block:: python

        # Single sequence
        >>> count_residues('PEPTIDE')
        {'P': 2, 'E': 2, 'T': 1, 'I': 1, 'D': 1}


    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _count_residues_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _count_residues_single(sequence)


def _percent_residues_single(
    sequence: str | ProFormaAnnotation, precision: int | None = None
) -> dict[str, float]:
    """Internal function for calculating percent residues in a single sequence."""
    return (
        get_annotation_input(sequence, copy=False)
        .condense_static_mods(inplace=True)
        .percent_residues(precision=precision)
    )


@overload
def percent_residues(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> dict[str, float]: ...


@overload
def percent_residues(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[dict[str, float]]: ...


def percent_residues(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> dict[str, float] | list[dict[str, float]]:
    """
    Calculates the percentage of each amino acid in the input sequence.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param precision: The precision of the percentage. Default is None.
    :type precision: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A dictionary or list of dictionaries containing the percentage of each amino acid.
    :rtype: dict[str, float] | list[dict[str, float]]

    .. code-block:: python

        # Single sequence
        >>> percent_residues('PEPTIDE', precision=2)
        {'P': 28.57, 'E': 28.57, 'T': 14.29, 'I': 14.29, 'D': 14.29}

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _percent_residues_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            precision=precision,
        )
    else:
        return _percent_residues_single(sequence, precision=precision)


def _count_aa_single(sequence: str | ProFormaAnnotation) -> dict[str, int]:
    """Internal function for counting amino acids in a single sequence."""
    annotation = get_annotation_input(sequence, copy=False)

    aa_counts = {aa: 0 for aa in ORDERED_AMINO_ACIDS}
    for aa in annotation.sequence:
        aa_counts[aa] += 1

    return aa_counts


@overload
def count_aa(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> dict[str, int]: ...


@overload
def count_aa(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[dict[str, int]]: ...


def count_aa(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> dict[str, int] | list[dict[str, int]]:
    """
    Converts a sequence to a feature vector.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: A dictionary or list of dictionaries with amino acid counts.
    :rtype: dict[str, int] | list[dict[str, int]]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _count_aa_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    else:
        return _count_aa_single(sequence)


def annotate_ambiguity(
    sequence: str | ProFormaAnnotation,
    forward_coverage: list[int],
    reverse_coverage: list[int],
    mass_shift: Any | None = None,
    add_mods_to_intervals: bool = False,
    sort_mods: bool = True,
    condense_to_xnotation: bool = False,
) -> str:
    """
    Given a peptide sequence and coverage information for forward and reverse fragment ions,
    annotate the sequence with ambiguity intervals and optional mass shift.

    This function identifies regions in the sequence where there is insufficient fragment ion
    coverage and marks them as ambiguous using ProForma notation with parentheses.
    If a mass shift is provided, it will be added to the appropriate location.

    Note: Mass shifts will be formatted with a '+' or '-' sign prefix.

    :param sequence: The peptide sequence or ProFormaAnnotation object to annotate.
    :type sequence: Union[str, ProFormaAnnotation]
    :param forward_coverage: Binary list indicating which positions have forward ion coverage (1) or not (0).
    :type forward_coverage: List[int]
    :param reverse_coverage: Binary list indicating which positions have reverse ion coverage (1) or not (0).
    :type reverse_coverage: List[int]
    :param mass_shift: An optional mass shift to be added to the sequence at the appropriate position.
    :type mass_shift: Optional[Any]
    :param add_mods_to_intervals: Whether to add modifications to interval annotations.
    :type add_mods_to_intervals: bool
    :param sort_mods: Whether to sort modifications.
    :type sort_mods: bool
    :param condense_to_xnotation: Whether to condense ambiguity to X notation.
    :type condense_to_xnotation: bool

    :raises ValueError: If the annotation already contains intervals or if coverage lengths don't match sequence length.

    :return: The annotated sequence with ambiguity regions marked using ProForma notation.
    :rtype: str

    .. code-block:: python

        # Add ambiguity intervals based on fragment ion coverage
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,1,0,0,0], [0,0,0,0,0,1,0])
        '(?PE)PTI(?DE)'

        # With a phosphorylation mass shift (note the '+' sign)
        >>> annotate_ambiguity('PEPTIDE', [1,1,1,0,0,0,0], [0,0,0,0,1,1,1], 79.966)
        'PEPT[+79.966]IDE'

        # Handling existing modifications
        >>> annotate_ambiguity('P[+10]EPTIDE', [1,1,1,0,0,0,0], [0,0,0,0,0,1,1])
        'P[+10]EP(?TI)DE'

        # When mass shift can't be localized to a specific residue
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,0,0,0,0], [0,0,0,0,0,1,0], 120)
        '(?PE)P(?TI)[+120](?DE)'

        # When mass shift is completely unlocalized, it becomes a labile modification
        >>> annotate_ambiguity('PEPTIDE', [0,1,1,1,1,0,0], [0,0,1,1,1,1,0], 120)
        '{+120}(?PE)PTI(?DE)'

        # Complex example with multiple intervals
        >>> for_ions = list(map(int, '00011101001000000000000000000000000000'))
        >>> rev_ions = list(map(int, '00000000000110000000101111111111010100'))
        >>> annotate_ambiguity('SSGSIASSYVQWYQQRPGSAPTTVIYEDDERPSGVPDR', for_ions, rev_ions, 120)
        '(?SSGS)IA(?SS)(?YVQ)W[+120](?YQQRPGSA)(?PT)TVIYEDDER(?PS)(?GV)(?PDR)'
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

    return annot.serialize()
