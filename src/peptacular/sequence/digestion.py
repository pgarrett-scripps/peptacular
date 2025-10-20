"""
Digest.py contains functions for generating enzymatic and non-enzymatic peptides from a given sequence. All functions
in this module are designed to work with both standard sequences and ProFormaAnnotations.

Since creating ProFormaAnnotation objects takes more resources than working with strings, ensure that return_type is
set to 'str' when working with large datasets to avoid unnecessary overhead from creating ProFormaAnnotation objects.

Valid DigestReturnType's:
    - 'str': Returns the digested sequences as strings.
    - 'annotation': Returns the digested sequences as ProFormaAnnotations.
    - 'span': Returns the digested sequences as Span objects.
    - 'str-span': Returns the digested sequences as a tuple of strings and Span objects.
    - 'annotation-span': Returns the digested sequences as a tuple of ProFormaAnnotations and Span objects.

"""

from typing import Literal, overload
from collections.abc import Sequence

from ..digestion import EnzymeConfig, DigestReturnType

from ..proforma.annotation import ProFormaAnnotation
from . import get_annotation_input
from .parrallel import parallel_apply_internal
from ..constants import ParrallelMethodLiteral, ParrallelMethod


def _get_left_semi_enzymatic_sequences_single(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[str]:
    """Get left semi-enzymatic sequences for a single sequence"""
    return list(
        get_annotation_input(sequence, copy=False).get_left_semi_enzymatic_sequences(
            min_len=min_len,
            max_len=max_len,
            return_type=DigestReturnType.STR,
        )
    )


@overload
def get_left_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def get_left_semi_enzymatic_sequences(
    sequence: Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def get_left_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Builds all left-hand semi-enzymatic subsequences derived from the input `sequence`.

    Automatically uses parallel processing when a list of sequences is provided.
    When method=None (default), automatically detects if GIL is disabled and uses
    threading for better performance, otherwise uses multiprocessing.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param min_len: Minimum length for the subsequences (inclusive). If None, defaults to 1.
    :type min_len: int | None
    :param max_len: Maximum length for the subsequences (inclusive). If None, goes up to len(sequence) - 1.
    :type max_len: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: The left-hand semi-enzymatic subsequences as strings, or list of lists.
    :rtype: list[str] | list[list[str]]

    .. code-block:: python

        >>> get_left_semi_enzymatic_sequences('PEPTIDE')
        ['PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _get_left_semi_enzymatic_sequences_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _get_left_semi_enzymatic_sequences_single(
            sequence=sequence,
            min_len=min_len,
            max_len=max_len,
        )


def _get_right_semi_enzymatic_sequences_single(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[str]:
    """Get right semi-enzymatic sequences for a single sequence"""
    return list(
        get_annotation_input(sequence, copy=False).get_right_semi_enzymatic_sequences(
            min_len=min_len,
            max_len=max_len,
            return_type=DigestReturnType.STR,
        )
    )


@overload
def get_right_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def get_right_semi_enzymatic_sequences(
    sequence: Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def get_right_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Builds all right-hand semi-enzymatic subsequences derived from the input `sequence`.

    Automatically uses parallel processing when a list of sequences is provided.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param min_len: Minimum length for the subsequences (inclusive). If None, defaults to 1.
    :type min_len: int | None
    :param max_len: Maximum length for the subsequences (inclusive). If None, goes up to len(sequence) - 1.
    :type max_len: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: The right-hand semi-enzymatic subsequences as strings, or list of lists.
    :rtype: list[str] | list[list[str]]

    .. code-block:: python

        >>> get_right_semi_enzymatic_sequences('PEPTIDE')
        ['EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _get_right_semi_enzymatic_sequences_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _get_right_semi_enzymatic_sequences_single(
            sequence=sequence,
            min_len=min_len,
            max_len=max_len,
        )


def _get_semi_enzymatic_sequences_single(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[str]:
    """Get semi-enzymatic sequences for a single sequence"""
    return list(
        get_annotation_input(sequence, copy=False).get_semi_enzymatic_sequences(
            min_len=min_len,
            max_len=max_len,
            return_type=DigestReturnType.STR,
        )
    )


@overload
def get_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def get_semi_enzymatic_sequences(
    sequence: Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def get_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Builds all semi-enzymatic sequences from the given input `sequence`.
    Equivalent to combining left and right semi-enzymatic sequences.

    Automatically uses parallel processing when a list of sequences is provided.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param min_len: Minimum length for the subsequences (inclusive). If None, defaults to 1.
    :type min_len: int | None
    :param max_len: Maximum length for the subsequences (inclusive). If None, goes up to len(sequence) - 1.
    :type max_len: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: Semi-enzymatic subsequences as strings, or list of lists.
    :rtype: list[str] | list[list[str]]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _get_semi_enzymatic_sequences_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _get_semi_enzymatic_sequences_single(
            sequence=sequence,
            min_len=min_len,
            max_len=max_len,
        )


def _get_non_enzymatic_sequences_single(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[str]:
    """Get non-enzymatic sequences for a single sequence"""
    return list(
        get_annotation_input(sequence, copy=False).get_non_enzymatic_sequences(
            min_len=min_len,
            max_len=max_len,
            return_type=DigestReturnType.STR,
        )
    )


@overload
def get_non_enzymatic_sequences(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def get_non_enzymatic_sequences(
    sequence: Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def get_non_enzymatic_sequences(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Builds all non-enzymatic sequences from the given input `sequence`.

    Automatically uses parallel processing when a list of sequences is provided.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param min_len: Minimum length for the subsequences (inclusive). If None, defaults to 1.
    :type min_len: int | None
    :param max_len: Maximum length for the subsequences (inclusive). If None, goes up to len(sequence) - 1.
    :type max_len: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: Non-enzymatic subsequences as strings, or list of lists.
    :rtype: list[str] | list[list[str]]

    .. code-block:: python

        >>> get_non_enzymatic_sequences('PEP')
        ['P', 'PE', 'E', 'EP', 'P']
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _get_non_enzymatic_sequences_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _get_non_enzymatic_sequences_single(
            sequence=sequence,
            min_len=min_len,
            max_len=max_len,
        )


def _get_cleavage_sites_single(
    sequence: str | ProFormaAnnotation, enzyme_regex: str
) -> list[int]:
    """Get cleavage sites for a single sequence"""
    return list(
        get_annotation_input(sequence, copy=False).get_cleavage_sites(
            enzyme_regex=enzyme_regex
        )
    )


@overload
def get_cleavage_sites(
    sequence: str | ProFormaAnnotation,
    enzyme_regex: str,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[int]: ...


@overload
def get_cleavage_sites(
    sequence: Sequence[str | ProFormaAnnotation],
    enzyme_regex: str,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[int]]: ...


def get_cleavage_sites(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    enzyme_regex: str,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[int] | list[list[int]]:
    """
    Return positions where cleavage occurs in input `sequence` based on the provided enzyme regex.

    Automatically uses parallel processing when a list of sequences is provided.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param enzyme_regex: Regular expression defining enzyme's cleavage rule.
    :type enzyme_regex: str
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: List of positions where cleavage occurs, or list of lists.
    :rtype: list[int] | list[list[int]]

    .. code-block:: python

        >>> get_cleavage_sites('TIDERTIDEKTIDE', 'trypsin')
        [5, 10]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _get_cleavage_sites_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            enzyme_regex=enzyme_regex,
        )
    else:
        return _get_cleavage_sites_single(
            sequence=sequence,
            enzyme_regex=enzyme_regex,
        )


def _digest_by_regex_single(
    sequence: str | ProFormaAnnotation,
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    complete_digestion: bool = True,
    sort_output: bool = True,
) -> list[str]:
    """Digest by regex for a single sequence"""
    return list(
        get_annotation_input(sequence, copy=False).regex_digest(
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
            return_type=DigestReturnType.STR,
        )
    )


@overload
def digest_by_regex(
    sequence: str | ProFormaAnnotation,
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def digest_by_regex(
    sequence: Sequence[str | ProFormaAnnotation],
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def digest_by_regex(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    enzyme_regex: str,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Returns digested sequences using a regular expression to define cleavage sites.

    Automatically uses parallel processing when a list of sequences is provided.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param enzyme_regex: Regular expression defining enzyme's cleavage rules.
    :type enzyme_regex: str
    :param missed_cleavages: Maximum number of missed cleavages, defaults to 0.
    :type missed_cleavages: int
    :param semi: Whether to include semi-enzymatic peptides, defaults to False.
    :type semi: bool
    :param min_len: Minimum length for the subsequences. If None, defaults to 1.
    :type min_len: int | None
    :param max_len: Maximum length for the subsequences. If None, no upper limit.
    :type max_len: int | None
    :param complete_digestion: If True, excludes original sequence. If False, includes it.
    :type complete_digestion: bool
    :param sort_output: Whether to sort the output sequences, defaults to True.
    :type sort_output: bool
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: List of digested peptides as strings, or list of lists.
    :rtype: list[str] | list[list[str]]

    .. code-block:: python

        >>> digest_by_regex('TIDERTIDEKTIDE', '([KR])', missed_cleavages=1)
        ['TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        >>> digest_by_regex('PEPCTIDE', '(?=C)', missed_cleavages=0)
        ['PEP', 'CTIDE']
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _digest_by_regex_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
        )
    else:
        return _digest_by_regex_single(
            sequence=sequence,
            enzyme_regex=enzyme_regex,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
        )


def _digest_single(
    sequence: str | ProFormaAnnotation,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    complete_digestion: bool = True,
    sort_output: bool = True,
) -> list[str]:
    """Digest for a single sequence"""
    return list(
        get_annotation_input(sequence, copy=False).digest(
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
            return_type=DigestReturnType.STR,
        )
    )


@overload
def digest(
    sequence: str | ProFormaAnnotation,
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def digest(
    sequence: Sequence[str | ProFormaAnnotation],
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def digest(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    cleave_on: str,
    restrict_before: str = "",
    restrict_after: str = "",
    cterminal: bool = True,
    missed_cleavages: int = 0,
    semi: bool = False,
    min_len: int | None = None,
    max_len: int | None = None,
    *,
    complete_digestion: bool = True,
    sort_output: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Returns digested sequences using amino acid specifications with optional restrictions.

    Automatically uses parallel processing when a list of sequences is provided.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences.
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param cleave_on: Amino acids to cleave on (e.g., 'KR' for trypsin-like).
    :type cleave_on: str
    :param restrict_before: Amino acids that prevent cleavage when preceding cleavage site.
    :type restrict_before: str
    :param restrict_after: Amino acids that prevent cleavage when following cleavage site.
    :type restrict_after: str
    :param cterminal: If True, cleave after the residue. If False, cleave before.
    :type cterminal: bool
    :param missed_cleavages: Maximum number of missed cleavages, defaults to 0.
    :type missed_cleavages: int
    :param semi: Whether to include semi-enzymatic peptides, defaults to False.
    :type semi: bool
    :param min_len: Minimum length for the subsequences. If None, defaults to 1.
    :type min_len: int | None
    :param max_len: Maximum length for the subsequences. If None, no upper limit.
    :type max_len: int | None
    :param complete_digestion: If True, excludes original sequence. If False, includes it.
    :type complete_digestion: bool
    :param sort_output: Whether to sort the output sequences, defaults to True.
    :type sort_output: bool
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: List of digested peptides as strings, or list of lists.
    :rtype: list[str] | list[list[str]]

    .. code-block:: python

        >>> digest('TIDERTIDEKTIDE', 'KR', cterminal=True)
        ['TIDER', 'TIDEK', 'TIDE']

        >>> digest('KPKRKP', 'K', restrict_after='P', cterminal=True)
        ['KPK', 'RKP']
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _digest_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
        )
    else:
        return _digest_single(
            sequence=sequence,
            cleave_on=cleave_on,
            restrict_before=restrict_before,
            restrict_after=restrict_after,
            cterminal=cterminal,
            missed_cleavages=missed_cleavages,
            semi=semi,
            min_len=min_len,
            max_len=max_len,
            complete_digestion=complete_digestion,
            sort_output=sort_output,
        )


def _digest_from_config_single(
    sequence: str | ProFormaAnnotation,
    config: EnzymeConfig,
    min_len: int | None = None,
    max_len: int | None = None,
    sort_output: bool = True,
) -> list[str]:
    """Digest from config for a single sequence"""
    return _digest_by_regex_single(
        sequence=sequence,
        enzyme_regex=config.enzyme_regex,
        missed_cleavages=config.missed_cleavages,
        semi=config.semi_enzymatic,
        min_len=min_len,
        max_len=max_len,
        complete_digestion=config.complete_digestion,
        sort_output=sort_output,
    )


@overload
def digest_from_config(
    sequence: str | ProFormaAnnotation,
    config: EnzymeConfig,
    min_len: int | None = None,
    max_len: int | None = None,
    sort_output: bool = True,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def digest_from_config(
    sequence: Sequence[str | ProFormaAnnotation],
    config: EnzymeConfig,
    min_len: int | None = None,
    max_len: int | None = None,
    sort_output: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def digest_from_config(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    config: EnzymeConfig,
    min_len: int | None = None,
    max_len: int | None = None,
    sort_output: bool = True,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Same as digest() but with a simplified configuration object for a single enzyme.

    Automatically uses parallel processing when a list of sequences is provided.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _digest_from_config_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            config=config,
            min_len=min_len,
            max_len=max_len,
            sort_output=sort_output,
        )
    else:
        return _digest_from_config_single(
            sequence=sequence,
            config=config,
            min_len=min_len,
            max_len=max_len,
            sort_output=sort_output,
        )


def _sequential_digest_single(
    sequence: str | ProFormaAnnotation,
    enzyme_configs: list[EnzymeConfig],
    min_len: int | None = None,
    max_len: int | None = None,
) -> list[str]:
    """Sequential digest for a single sequence"""
    return list(
        get_annotation_input(sequence, copy=False).sequential_digest(
            enzyme_configs=enzyme_configs,
            min_len=min_len,
            max_len=max_len,
            return_type=DigestReturnType.STR,
        )
    )


@overload
def sequential_digest(
    sequence: str | ProFormaAnnotation,
    enzyme_configs: list[EnzymeConfig],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str]: ...


@overload
def sequential_digest(
    sequence: Sequence[str | ProFormaAnnotation],
    enzyme_configs: list[EnzymeConfig],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[list[str]]: ...


def sequential_digest(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    enzyme_configs: list[EnzymeConfig],
    min_len: int | None = None,
    max_len: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[str] | list[list[str]]:
    """
    Returns digested sequences using sequential digestion with multiple enzymes.

    Automatically uses parallel processing when a list of sequences is provided.

    :param sequence: A sequence, ProFormaAnnotation, or list of sequences
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param enzyme_configs: List of EnzymeConfig objects to apply sequentially
    :type enzyme_configs: list[EnzymeConfig]
    :param min_len: Minimum length for peptides. If None, defaults to 1.
    :type min_len: int | None
    :param max_len: Maximum length for peptides. If None, no upper limit.
    :type max_len: int | None
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None

    :return: List of digested peptides as strings, or list of lists.
    :rtype: list[str] | list[list[str]]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _sequential_digest_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            enzyme_configs=enzyme_configs,
            min_len=min_len,
            max_len=max_len,
        )
    else:
        return _sequential_digest_single(
            sequence=sequence,
            enzyme_configs=enzyme_configs,
            min_len=min_len,
            max_len=max_len,
        )
