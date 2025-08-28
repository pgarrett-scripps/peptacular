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

from typing import Generator

from ..digestion import EnzymeConfig, DigestReturnType

from ..proforma.annotation import ProFormaAnnotation
from . import get_annotation_input


def get_left_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[str, None, None]:
    """
    Builds all left-hand semi-enzymatic subsequences derived from the input `sequence`.

    :param sequence: A sequence or ProFormaAnnotation.
    :param min_len: Minimum length for the subsequences (inclusive). If None, defaults to 1.
    :param max_len: Maximum length for the subsequences (inclusive). If None, goes up to len(sequence) - 1.

    :return: The left-hand semi-enzymatic subsequences as strings.

    .. code-block:: python

        >>> list(get_left_semi_enzymatic_sequences('PEPTIDE'))
        ['PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']
    """
    return get_annotation_input(sequence, copy=False).get_left_semi_enzymatic_sequences(
        min_len=min_len,
        max_len=max_len,
        return_type=DigestReturnType.STR,
    )


def get_right_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[str, None, None]:
    """
    Builds all right-hand semi-enzymatic subsequences derived from the input `sequence`.

    :param sequence: A sequence or ProFormaAnnotation.
    :param min_len: Minimum length for the subsequences (inclusive). If None, defaults to 1.
    :param max_len: Maximum length for the subsequences (inclusive). If None, goes up to len(sequence) - 1.

    :return: The right-hand semi-enzymatic subsequences as strings.

    .. code-block:: python

        >>> list(get_right_semi_enzymatic_sequences('PEPTIDE'))
        ['EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']
    """

    return get_annotation_input(
        sequence, copy=False
    ).get_right_semi_enzymatic_sequences(
        min_len=min_len,
        max_len=max_len,
        return_type=DigestReturnType.STR,
    )


def get_semi_enzymatic_sequences(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[str, None, None]:
    """
    Builds all semi-enzymatic sequences from the given input `sequence`.
    Equivalent to combining left and right semi-enzymatic sequences.

    :param sequence: A sequence or ProFormaAnnotation.
    :param min_len: Minimum length for the subsequences (inclusive). If None, defaults to 1.
    :param max_len: Maximum length for the subsequences (inclusive). If None, goes up to len(sequence) - 1.

    :return: Semi-enzymatic subsequences as strings.
    """

    return get_annotation_input(sequence, copy=False).get_semi_enzymatic_sequences(
        min_len=min_len,
        max_len=max_len,
        return_type=DigestReturnType.STR,
    )


def get_non_enzymatic_sequences(
    sequence: str | ProFormaAnnotation,
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[str, None, None]:
    """
    Builds all non-enzymatic sequences from the given input `sequence`.

    :param sequence: A sequence or ProFormaAnnotation.
    :param min_len: Minimum length for the subsequences (inclusive). If None, defaults to 1.
    :param max_len: Maximum length for the subsequences (inclusive). If None, goes up to len(sequence) - 1.

    :return: Non-enzymatic subsequences as strings.

    .. code-block:: python

        >>> list(get_non_enzymatic_sequences('PEP'))
        ['P', 'PE', 'E', 'EP', 'P']
    """
    return get_annotation_input(sequence, copy=False).get_non_enzymatic_sequences(
        min_len=min_len,
        max_len=max_len,
        return_type=DigestReturnType.STR,
    )


def get_cleavage_sites(
    sequence: str | ProFormaAnnotation, enzyme_regex: str
) -> Generator[int, None, None]:
    """
    Return positions where cleavage occurs in input `sequence` based on the provided enzyme regex.

    :param sequence: A sequence or ProFormaAnnotation.
    :param enzyme_regex: Regular expression defining enzyme's cleavage rule.

    :return: List of positions where cleavage occurs.

    .. code-block:: python

        >>> list(get_cleavage_sites('TIDERTIDEKTIDE', 'trypsin'))
        [5, 10]
    """
    return get_annotation_input(sequence, copy=False).get_cleavage_sites(
        enzyme_regex=enzyme_regex
    )


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
) -> Generator[str, None, None]:
    """
    Returns digested sequences using a regular expression to define cleavage sites.

    :param sequence: A sequence or ProFormaAnnotation.
    :param enzyme_regex: Regular expression defining enzyme's cleavage rules.
    :param missed_cleavages: Maximum number of missed cleavages, defaults to 0.
    :param semi: Whether to include semi-enzymatic peptides, defaults to False.
    :param min_len: Minimum length for the subsequences. If None, defaults to 1.
    :param max_len: Maximum length for the subsequences. If None, no upper limit.
    :param complete_digestion: If True, excludes original sequence. If False, includes it.
    :param sort_output: Whether to sort the output sequences, defaults to True.

    :return: List of digested peptides as strings.

    .. code-block:: python

        >>> list(digest_by_regex('TIDERTIDEKTIDE', '([KR])', missed_cleavages=1))
        ['TIDER', 'TIDERTIDEK', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        >>> list(digest_by_regex('PEPCTIDE', '(?=C)', missed_cleavages=0))
        ['PEP', 'CTIDE']
    """
    return get_annotation_input(sequence, copy=False).regex_digest(
        enzyme_regex=enzyme_regex,
        missed_cleavages=missed_cleavages,
        semi=semi,
        min_len=min_len,
        max_len=max_len,
        complete_digestion=complete_digestion,
        sort_output=sort_output,
        return_type=DigestReturnType.STR,
    )


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
) -> Generator[str, None, None]:
    """
    Returns digested sequences using amino acid specifications with optional restrictions.

    :param sequence: A sequence or ProFormaAnnotation.
    :param cleave_on: Amino acids to cleave on (e.g., 'KR' for trypsin-like).
    :param restrict_before: Amino acids that prevent cleavage when preceding cleavage site.
    :param restrict_after: Amino acids that prevent cleavage when following cleavage site.
    :param cterminal: If True, cleave after the residue. If False, cleave before.
    :param missed_cleavages: Maximum number of missed cleavages, defaults to 0.
    :param semi: Whether to include semi-enzymatic peptides, defaults to False.
    :param min_len: Minimum length for the subsequences. If None, defaults to 1.
    :param max_len: Maximum length for the subsequences. If None, no upper limit.
    :param complete_digestion: If True, excludes original sequence. If False, includes it.
    :param sort_output: Whether to sort the output sequences, defaults to True.

    :return: List of digested peptides as strings.

    .. code-block:: python

        >>> list(digest('TIDERTIDEKTIDE', 'KR', cterminal=True))
        ['TIDER', 'TIDEK', 'TIDE']

        >>> list(digest('KPKRKP', 'K', restrict_after='P', cterminal=True))
        ['KPK', 'RKP']
    """
    return get_annotation_input(sequence, copy=False).digest(
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


def digest_from_config(
    sequence: str | ProFormaAnnotation,
    config: EnzymeConfig,
    min_len: int | None = None,
    max_len: int | None = None,
    sort_output: bool = True,
) -> Generator[str, None, None]:
    """
    Same as digest() but with a simplified configuration object for a single enzyme.
    """

    return digest_by_regex(
        sequence=sequence,
        enzyme_regex=config.enzyme_regex,
        missed_cleavages=config.missed_cleavages,
        semi=config.semi_enzymatic,
        min_len=min_len,
        max_len=max_len,
        complete_digestion=config.complete_digestion,
        sort_output=sort_output,
    )


def sequential_digest(
    sequence: str | ProFormaAnnotation,
    enzyme_configs: list[EnzymeConfig],
    min_len: int | None = None,
    max_len: int | None = None,
) -> Generator[str, None, None]:
    """
    Returns digested sequences using sequential digestion with multiple enzymes.

    :param sequence: A sequence or ProFormaAnnotation
    :param enzyme_configs: List of EnzymeConfig objects to apply sequentially
    :param min_len: Minimum length for peptides. If None, defaults to 1.
    :param max_len: Maximum length for peptides. If None, no upper limit.

    :return: List of digested peptides as strings.
    """
    return get_annotation_input(sequence, copy=False).sequential_digest(
        enzyme_configs=enzyme_configs,
        min_len=min_len,
        max_len=max_len,
        return_type=DigestReturnType.STR,
    )
