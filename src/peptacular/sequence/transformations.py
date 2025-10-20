from typing import Any, Callable, Sequence, overload

from .util import get_annotation_input
from ..proforma.annotation import (
    ProFormaAnnotation,
)
from .parrallel import parallel_apply_internal
from ..constants import ParrallelMethodLiteral, ParrallelMethod


def _reverse_single(
    sequence: str | ProFormaAnnotation,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """Internal function for reversing a single sequence."""
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .reverse(inplace=True, keep_nterm=keep_nterm, keep_cterm=keep_cterm)
        .serialize(include_plus=include_plus, precision=precision)
    )


@overload
def reverse(
    sequence: str | ProFormaAnnotation,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def reverse(
    sequence: Sequence[str | ProFormaAnnotation],
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def reverse(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Reverses the sequence, while preserving the position of any modifications.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param swap_terms: If True, the N- and C-terminal modifications will be swapped.
    :type swap_terms: bool
    :param keep_nterm: Number of N-terminal residues to keep unchanged. Default is 0.
    :type keep_nterm: int
    :param keep_cterm: Number of C-terminal residues to keep unchanged. Default is 0.
    :type keep_cterm: int
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The reversed sequence(s) with modifications preserved.
    :rtype: str | list[str]

    .. code-block:: python

        # Single sequence
        >>> reverse('PEPTIDE')
        'EDITPEP'

        # Keep first 2 residues unchanged
        >>> reverse('PEPTIDE', keep_nterm=2)
        'PEDITPE'

        # Keep first and last residue unchanged
        >>> reverse('PEPTIDE', keep_nterm=1, keep_cterm=1)
        'PDITPEE'

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> reversed_seqs = reverse(sequences, keep_nterm=1)
        >>> len(reversed_seqs)
        3

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _reverse_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
            include_plus=include_plus,
            precision=precision,
        )
    else:
        return _reverse_single(
            sequence, keep_nterm, keep_cterm, include_plus, precision
        )


def _shuffle_single(
    sequence: str | ProFormaAnnotation,
    seed: int | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """Internal function for shuffling a single sequence."""
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .shuffle(seed=seed, keep_nterm=keep_nterm, keep_cterm=keep_cterm, inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


@overload
def shuffle(
    sequence: str | ProFormaAnnotation,
    seed: int | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def shuffle(
    sequence: Sequence[str | ProFormaAnnotation],
    seed: int | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def shuffle(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    seed: int | None = None,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Shuffles the sequence, while preserving the position of any modifications.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param seed: Seed for the random number generator.
    :type seed: int | None
    :param keep_nterm: Number of N-terminal residues to keep unchanged. Default is 0.
    :type keep_nterm: int
    :param keep_cterm: Number of C-terminal residues to keep unchanged. Default is 0.
    :type keep_cterm: int
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The shuffled sequence(s) with modifications preserved.
    :rtype: str | list[str]

    .. code-block:: python

        # Single sequence
        >>> shuffle('PEPTIDE', seed=0)
        'IPEPDTE'

        # Keep first 2 residues unchanged
        >>> shuffle('PEPTIDE', seed=0, keep_nterm=2)
        'PEIPDTE'

        # Keep first and last residue unchanged
        >>> shuffle('PEPTIDE', seed=0, keep_nterm=1, keep_cterm=1)
        'PIPEDTE'

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> shuffled_seqs = shuffle(sequences, seed=0, keep_nterm=1)
        >>> len(shuffled_seqs)
        3

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _shuffle_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            seed=seed,
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
            include_plus=include_plus,
            precision=precision,
        )
    else:
        return _shuffle_single(
            sequence, seed, keep_nterm, keep_cterm, include_plus, precision
        )


def _shift_single(
    sequence: str | ProFormaAnnotation,
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
) -> str:
    """Internal function for shifting a single sequence."""
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .shift(n=n, keep_nterm=keep_nterm, keep_cterm=keep_cterm, inplace=True)
        .serialize(include_plus=include_plus)
    )


@overload
def shift(
    sequence: str | ProFormaAnnotation,
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def shift(
    sequence: Sequence[str | ProFormaAnnotation],
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def shift(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n: int,
    keep_nterm: int = 0,
    keep_cterm: int = 0,
    include_plus: bool = False,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Shifts the sequence to the left by a given number of positions, while preserving the position of any modifications.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param n: The number of positions to shift the sequence to the left.
    :type n: int
    :param keep_nterm: Number of N-terminal residues to keep unchanged. Default is 0.
    :type keep_nterm: int
    :param keep_cterm: Number of C-terminal residues to keep unchanged. Default is 0.
    :type keep_cterm: int
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The shifted sequence(s) with modifications preserved.
    :rtype: str | list[str]

    .. code-block:: python

        # Single sequence
        >>> shift('PEPTIDE', 2)
        'PTIDEPE'

        # Keep first 2 residues unchanged
        >>> shift('PEPTIDE', 2, keep_nterm=2)
        'PEPTIDE'

        # Keep first and last residue unchanged - shift only the middle
        >>> shift('PEPTIDE', 2, keep_nterm=1, keep_cterm=1)
        'PTIDEPE'

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> shifted_seqs = shift(sequences, 2, keep_nterm=1)
        >>> len(shifted_seqs)
        3

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _shift_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            n=n,
            keep_nterm=keep_nterm,
            keep_cterm=keep_cterm,
            include_plus=include_plus,
        )
    else:
        return _shift_single(sequence, n, keep_nterm, keep_cterm, include_plus)


def _span_to_sequence_single(
    sequence: str | ProFormaAnnotation,
    span: tuple[int, int, int],
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """Internal function for extracting span from a single sequence."""
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .slice(span[0], span[1], inplace=True)
        .serialize(include_plus=include_plus, precision=precision)
    )


@overload
def span_to_sequence(
    sequence: str | ProFormaAnnotation,
    span: tuple[int, int, int],
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def span_to_sequence(
    sequence: Sequence[str | ProFormaAnnotation],
    span: tuple[int, int, int],
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def span_to_sequence(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    span: tuple[int, int, int],
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Extracts a subsequence from the input sequence based on the provided span.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param span: A tuple representing the span of the subsequence to be extracted.
    :type span: tuple[int, int, int]
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The subsequence(s) of the input sequence defined by the span.
    :rtype: str | list[str]

    .. code-block:: python

        # Single sequence
        >>> span_to_sequence('PEPTIDE', (0, 4, 0))
        'PEPT'

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> spans = span_to_sequence(sequences, (1, 4, 0))
        >>> len(spans)
        3

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _span_to_sequence_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            span=span,
            include_plus=include_plus,
            precision=precision,
        )
    else:
        return _span_to_sequence_single(sequence, span, include_plus, precision)


def _split_single(
    sequence: str | ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
) -> list[str]:
    """Internal function for splitting a single sequence."""
    return [
        a.serialize(include_plus=include_plus, precision=precision)
        for a in get_annotation_input(sequence=sequence, copy=True).split()
    ]


@overload
def split(
    sequence: str | ProFormaAnnotation,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def split(
    sequence: Sequence[str | ProFormaAnnotation],
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def split(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Splits sequence into a list of amino acids, preserving modifications.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A list of amino acids or list of lists for multiple sequences.
    :rtype: list[str] | list[list[str]]

    .. code-block:: python

        # Single sequence
        >>> split('PEPTIDE')
        ['P', 'E', 'P', 'T', 'I', 'D', 'E']

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN']
        >>> split_seqs = split(sequences)
        >>> len(split_seqs)
        2
        >>> len(split_seqs[0])
        7

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _split_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            include_plus=include_plus,
            precision=precision,
        )
    else:
        return _split_single(sequence, include_plus, precision)


def _sort_single(
    sequence: str | ProFormaAnnotation,
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
) -> str:
    """Internal function for sorting a single sequence."""
    return (
        get_annotation_input(sequence=sequence, copy=True)
        .sort(inplace=True, key=key, reverse=reverse)
        .serialize(include_plus=include_plus, precision=precision)
    )


@overload
def sort(
    sequence: str | ProFormaAnnotation,
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str: ...


@overload
def sort(
    sequence: Sequence[str | ProFormaAnnotation],
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


def sort(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    key: Callable[[str], Any] | None = None,
    reverse: bool = False,
    include_plus: bool = False,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> str | list[str]:
    """
    Sorts the input sequence using the provided sort function. Terminal sequences are kept in place.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: The sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param key: Optional key function for sorting.
    :type key: Callable[[str], Any] | None
    :param reverse: If True, sort in reverse order.
    :type reverse: bool
    :param include_plus: Whether to include the plus sign for positive mass modifications.
    :type include_plus: bool
    :param precision: The precision of the mass. Default is None.
    :type precision: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: The sorted sequence(s).
    :rtype: str | list[str]

    .. code-block:: python

        # Single sequence
        >>> sort('PEPTIDE')
        'DEEIPPT'

        # Multiple sequences (automatic parallel processing)
        >>> sequences = ['PEPTIDE', 'PROTEIN', 'SEQUENCE']
        >>> sorted_seqs = sort(sequences)
        >>> len(sorted_seqs)
        3

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _sort_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            key=key,
            reverse=reverse,
            include_plus=include_plus,
            precision=precision,
        )
    else:
        return _sort_single(sequence, key, reverse, include_plus, precision)


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
