"""
Digest.py contains functions to handle all your digestion needs and of course it works with peptacular's standard
modification notation. From generating left/right semi enzymatic sequences to calculating protease cleavage sites
there is likely a function for it here.

Proteases can be a name of a protease in peptacular.constants.PROTEASES, or they can be a custom regex string. When
specifying more than one protease, all cleavage sites will be combined (as if both proteases were present at the
same time)
"""

from __future__ import annotations

from typing import Union, List

from peptacular.spans import Span
from peptacular.constants import PROTEASES_COMPILED
from peptacular.proforma.proforma import ProFormaAnnotation, create_annotation
from peptacular.spans import build_left_semi_spans, build_right_semi_spans, build_non_enzymatic_spans, build_spans
from peptacular.sequence.sequence import sequence_to_annotation
from peptacular.util import get_regex_match_indices


def _return_digest(annotation: ProFormaAnnotation, spans: List[Span],
                   return_str: bool) -> List[str] | List[ProFormaAnnotation]:
    if not annotation.has_mods():
        sequences = [annotation.sequence[span[0]: span[1]] for span in spans]

        if return_str:
            return sequences

        return [create_annotation(sequence) for sequence in sequences]

    annotations = [annotation.slice(span[0], span[1]) for span in spans]
    if return_str:
        return [annotation.serialize() for annotation in annotations]

    return annotations


def get_left_semi_enzymatic_sequences(sequence: str | ProFormaAnnotation,
                                      min_len: int | None = None,
                                      max_len: int | None = None,
                                      return_str: bool = True) -> List[str] | List[ProFormaAnnotation]:
    """
    Builds all left-hand semi-enzymatic subsequences derived from the input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param min_len: Minimum length for the subsequences (inclusive), defaults to [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]
    :param return_str: Whether to return the digested sequences as strings or as ProFormaAnnotations, defaults to [True].
    :type return_str: bool

    :return: Left-hand semi-enzymatic subsequences.
    :rtype: List[str]

    .. code-block:: python

        # Generates all left-hand semi enzymatic sequences (Returned values does not include input sequence)
        >>> get_left_semi_enzymatic_sequences('PEPTIDE')
        ['PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        # For single-letter or empty sequences, the function returns an empty list:
        >>> get_left_semi_enzymatic_sequences('P')
        []
        >>> get_left_semi_enzymatic_sequences('')
        []

        # Modifications are preserved:
        >>> get_left_semi_enzymatic_sequences('[1]-P[2]EPTIDE')
        ['[1]-P[2]EPTID', '[1]-P[2]EPTI', '[1]-P[2]EPT', '[1]-P[2]EP', '[1]-P[2]E', '[1]-P[2]']

        >>> get_left_semi_enzymatic_sequences('<13C>TIDE')
        ['<13C>TID', '<13C>TI', '<13C>T']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
        input_type = str
    elif isinstance(sequence, ProFormaAnnotation):
        annotation = sequence
        input_type = ProFormaAnnotation
    else:
        raise ValueError(f"Unsupported input type: {type(sequence)}")

    s = (0, len(annotation), 0)
    spans = build_left_semi_spans(span=s, min_len=min_len, max_len=max_len)

    return _return_digest(annotation, spans, return_str)


def get_right_semi_enzymatic_sequences(sequence: str | ProFormaAnnotation,
                                       min_len: int | None = None,
                                       max_len: int | None = None,
                                       return_str: bool = True) -> List[str]:
    """
    Builds all right-hand semi-enzymatic subsequences derived from the input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]
    :param return_str: Whether to return the digested sequences as strings or as ProFormaAnnotations, defaults to [True].
    :type return_str: bool

    :return: Right-hand semi-enzymatic subsequences
    :rtype: List[str]

    .. code-block:: python

        # Generates all right-hand semi enzymatic sequences (Returned values does not include input sequence)
        >>> get_right_semi_enzymatic_sequences('PEPTIDE')
        ['EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        # For single-letter or empty sequences, the function returns an empty list:
        >>> get_right_semi_enzymatic_sequences('P')
        []
        >>> get_right_semi_enzymatic_sequences('')
        []

        # Modifications are preserved:
        >>> get_right_semi_enzymatic_sequences('PEPTIDE[1]-[2]')
        ['EPTIDE[1]-[2]', 'PTIDE[1]-[2]', 'TIDE[1]-[2]', 'IDE[1]-[2]', 'DE[1]-[2]', 'E[1]-[2]']

        >>> get_right_semi_enzymatic_sequences('<13C>TIDE')
        ['<13C>IDE', '<13C>DE', '<13C>E']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
        input_type = str
    elif isinstance(sequence, ProFormaAnnotation):
        annotation = sequence
        input_type = ProFormaAnnotation
    else:
        raise ValueError(f"Unsupported input type: {type(sequence)}")

    s = (0, len(annotation), 0)
    spans = build_right_semi_spans(span=s, min_len=min_len, max_len=max_len)

    return _return_digest(annotation, spans, return_str)


def get_semi_enzymatic_sequences(sequence: str | ProFormaAnnotation,
                                 min_len: int | None = None,
                                 max_len: int | None = None,
                                 return_str: bool = True) -> List[str] | List[ProFormaAnnotation]:
    """
    Builds allsemi-enzymatic sequences from the given input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up  to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]

    :return: Semi-enzymatic subsequences.
    :rtype: List[str]

    .. code-block:: python

        # Equivalent to build_left_semi_sequences + build_right_semi_sequences
        >>> res = get_left_semi_enzymatic_sequences('PEPTIDE') + get_right_semi_enzymatic_sequences('PEPTIDE')
        >>> get_semi_enzymatic_sequences('PEPTIDE') == res
        True

    """

    return get_left_semi_enzymatic_sequences(sequence=sequence, min_len=min_len, max_len=max_len) + \
        get_right_semi_enzymatic_sequences(sequence=sequence, min_len=min_len, max_len=max_len)


def get_non_enzymatic_sequences(sequence: str | ProFormaAnnotation,
                                min_len: int | None = None,
                                max_len: int | None = None,
                                return_str: bool = True) -> List[str] | List[ProFormaAnnotation]:
    """
    Builds all non-enzymatic sequences from the given input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up  to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]
    :param return_str: Whether to return the digested sequences as strings or as ProFormaAnnotations, defaults to [True].
    :type return_str: bool

    :return: Non-enzymatic subsequences
    :rtype: List[str]

    .. code-block:: python

        # Generates non-enzymatic sequences (Returned values does not include input sequence):
        >>> get_non_enzymatic_sequences('PEP')
        ['P', 'PE', 'E', 'EP', 'P']

        # For single-letter or empty sequences, the function returns an empty list:
        >>> get_non_enzymatic_sequences('P')
        []
        >>> get_non_enzymatic_sequences('')
        []

         # Sequences with modifications are processed preserving those modifications:
        >>> get_non_enzymatic_sequences('[Acetyl]-P[1.0]EP[1.0]-[Amide]')
        ['[Acetyl]-P[1.0]', '[Acetyl]-P[1.0]E', 'E', 'EP[1.0]-[Amide]', 'P[1.0]-[Amide]']

        >>> get_non_enzymatic_sequences('<13C>PEP')
        ['<13C>P', '<13C>PE', '<13C>E', '<13C>EP', '<13C>P']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    elif isinstance(sequence, ProFormaAnnotation):
        annotation = sequence
    else:
        raise ValueError(f"Unsupported input type: {type(sequence)}")

    s = (0, len(annotation), 0)
    spans = build_non_enzymatic_spans(span=s, min_len=min_len, max_len=max_len)
    return _return_digest(annotation, spans, return_str)


def get_enzymatic_sequences(sequence: str | ProFormaAnnotation,
                            enzyme_regex: str,
                            missed_cleavages: int = 0,
                            semi: bool = False,
                            min_len: int | None = None,
                            max_len: int | None = None,
                            return_str: bool = True) -> List[str] | List[ProFormaAnnotation]:
    """
    Builds all enzymatic sequences from the given input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param enzyme_regex: Regular expression or key in PROTEASES dictionary defining the enzyme's cleavage rule.
    :type enzyme_regex: str
    :param missed_cleavages: Maximum number of missed cleavages.
    :type missed_cleavages: int
    :param semi: Whether to include semi-enzymatic peptides, defaults to [False].
    :type semi: bool
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up  to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]
    :param return_str: Whether to return the digested sequences as strings or as ProFormaAnnotations, defaults to [True].
    :type return_str: bool

    :return: Enzymatic subsequences.
    :rtype: List[str]

    .. code-block:: python

        # Can specify the name of a protease in the PROTEASES dictionary:
        >>> get_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        # Or specify a regular expression:
        >>> get_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        # Filter sequences by max length:
        >>> get_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, max_len=5)
        ['TIDER', 'TIDEK', 'TIDE']

        # Filter sequences by min length:
        >>> get_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, min_len=6)
        ['TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEKTIDE']

        # Modifications are preserved:
        >>> get_enzymatic_sequences(sequence='[1]-TIDERT[1.0]IDEKTI-[2]', enzyme_regex='([KR])', missed_cleavages=2)
        ['[1]-TIDER', '[1]-TIDERT[1.0]IDEK', '[1]-TIDERT[1.0]IDEKTI-[2]', 'T[1.0]IDEK', 'T[1.0]IDEKTI-[2]', 'TI-[2]']

        >>> get_enzymatic_sequences(sequence='<13C>TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2)
        ['<13C>TIDER', '<13C>TIDERTIDEK', '<13C>TIDERTIDEKTIDE', '<13C>TIDEK', '<13C>TIDEKTIDE', '<13C>TIDE']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
        input_type = str
    elif isinstance(sequence, ProFormaAnnotation):
        annotation = sequence
        input_type = ProFormaAnnotation
    else:
        raise ValueError(f"Unsupported input type: {type(sequence)}")

    cleavage_sites = [0] + get_cleavage_sites(sequence=annotation, enzyme_regex=enzyme_regex) + [len(annotation)]
    spans = build_spans(max_index=len(annotation), enzyme_sites=cleavage_sites, missed_cleavages=missed_cleavages,
                        min_len=min_len, max_len=max_len, semi=semi)

    return _return_digest(annotation, spans, return_str)


def get_cleavage_sites(sequence: str | ProFormaAnnotation, enzyme_regex: str) -> List[int]:
    """
    Return a list of positions where cleavage occurs in input `sequence` based on the provided enzyme regex.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param enzyme_regex: Regular expression or key in PROTEASES dictionary defining enzyme's cleavage rule.
    :type enzyme_regex: str

    :return: List of positions where cleavage occurs in the sequence.
    :rtype: List[int]

    .. code-block:: python

        # Can use a key in PROTEASES to specify the enzyme_regex:
        >>> get_cleavage_sites(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P')
        [5, 10]

        # Or specify a regular expression:
        >>> get_cleavage_sites(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])')
        [5, 10]

        # No cleavage sites are identified if the enzyme_regex does not match the sequence:
        >>> get_cleavage_sites(sequence='TIDEPTIDEPTIDE', enzyme_regex='trypsin/P')
        []

        # If the protease cleaves at the N-terminus, the first position is included:
        >>> get_cleavage_sites(sequence='KPEPTIDEK', enzyme_regex='lys-n')
        [8]

        # Similarly, if the protease cleaves at the C-terminus, the last position is included:
        >>> get_cleavage_sites(sequence='KPEPTIDEK', enzyme_regex='lys-c')
        [1]

        # Will also work with modified sequences
        >>> get_cleavage_sites(sequence='[Acetyl]-TIDERT[1.0]IDEKTIDE-[Amide]', enzyme_regex='trypsin/P')
        [5, 10]

        # Non-specific cleavage sites are also identified
        >>> get_cleavage_sites(sequence='PEPTIDE', enzyme_regex='non-specific')
        [1, 2, 3, 4, 5, 6]

        >>> get_cleavage_sites(sequence='PPPP', enzyme_regex='PP')
        [1, 2, 3]

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    if enzyme_regex in PROTEASES_COMPILED:
        enzyme_regex = PROTEASES_COMPILED[enzyme_regex]

    last_index = len(annotation)
    return [i + 1 for i in get_regex_match_indices(input_str=annotation.sequence, regex_str=enzyme_regex)
            if i < last_index - 1]


def digest(sequence: str | ProFormaAnnotation,
           enzyme_regex: List[str] | str,
           missed_cleavages: int = 0,
           semi: bool = False,
           min_len: int | None = None,
           max_len: int | None = None,
           return_str: bool = True) -> List[str] | List[ProFormaAnnotation]:
    """
    Returns a list of digested sequences derived from the input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param enzyme_regex: Regular expression or list of regular expressions representing enzyme's cleavage rules.
    :type enzyme_regex: Union[List[str], str]
    :param missed_cleavages: Maximum number of missed cleavages, defaults to [0].
    :type missed_cleavages: int
    :param semi: Whether to include semi-enzymatic peptides, defaults to [False].
    :type semi: bool
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up  to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]
    :param return_str: Whether to return the digested sequences as strings or as ProFormaAnnotations, defaults to [True].
    :type return_str: bool

    :return: List of digested peptides.
    :rtype: List[str]

    .. code-block:: python

        # Can use a key in PROTEASES to specify the enzyme_regex:
        >>> digest(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        # Or specify a regular expression:
        >>> digest(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        # Filter sequecnes by max length:
        >>> digest(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, max_len=5)
        ['TIDER', 'TIDEK', 'TIDE']

        # Filter sequecnes by min length:
        >>> digest(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, min_len=6)
        ['TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEKTIDE']

        # Generate semi-enzymatic sequences:
        >>> digest(sequence='TIDERTIDEK[1]TIDE-[2]', enzyme_regex='([KR])', missed_cleavages=1, min_len=9, semi=True)
        ['TIDERTIDEK[1]', 'TIDEK[1]TIDE-[2]', 'TIDERTIDE', 'IDERTIDEK[1]']

        # Generate sequences with a non-specific enzyme
        >>> digest(sequence='PEPT', enzyme_regex='non-specific')
        ['P', 'PE', 'PEP', 'E', 'EP', 'EPT', 'P', 'PT', 'T']

        # Modifications are preserved:
        >>> digest(sequence='<13C>T[1][2]IDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, min_len=6)
        ['<13C>T[1][2]IDERTIDEK', '<13C>T[1][2]IDERTIDEKTIDE', '<13C>TIDEKTIDE']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    elif isinstance(sequence, ProFormaAnnotation):
        annotation = sequence
    else:
        raise ValueError(f"Unsupported input type: {type(sequence)}")

    seq_len = len(annotation)

    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = seq_len

    if isinstance(enzyme_regex, str):
        enzyme_regex = [enzyme_regex]

    cleavage_sites = []
    for regex in enzyme_regex:
        cleavage_sites.extend(get_cleavage_sites(sequence=annotation, enzyme_regex=regex))

    spans = build_spans(max_index=seq_len, enzyme_sites=cleavage_sites, missed_cleavages=missed_cleavages,
                        min_len=min_len, max_len=max_len, semi=semi)

    return _return_digest(annotation, spans, return_str)
