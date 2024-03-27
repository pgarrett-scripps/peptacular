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

from typing import Union, List, Optional, TypeAlias, Literal, Tuple

from peptacular.spans import Span
from peptacular.constants import PROTEASES_COMPILED
from peptacular.proforma.proforma import ProFormaAnnotation, create_annotation
from peptacular.spans import build_left_semi_spans, build_right_semi_spans, build_non_enzymatic_spans, build_spans
from peptacular.sequence.sequence import sequence_to_annotation
from peptacular.util import get_regex_match_indices

DigestReturnType: TypeAlias = Literal["str", "annotation", "span", "str-span", "annotation-span"]

DIGEST_RETURN_TYPING = Union[List[str], List[ProFormaAnnotation], List[Span], List[Tuple[str, Span]],
                             List[Tuple[ProFormaAnnotation, Span]]]


def _return_digested_sequences(annotation: ProFormaAnnotation,
                               spans: List[Span],
                               return_type: DigestReturnType) -> DIGEST_RETURN_TYPING:
    """
    Helper function to return the digested sequences as strings or ProFormaAnnotations.
    """

    if return_type == 'span':
        return spans

    elif return_type == 'annotation':
        if not annotation.has_mods():
            return [create_annotation(annotation.sequence[span[0]: span[1]]) for span in spans]

        return [annotation.slice(span[0], span[1]) for span in spans]

    elif return_type == 'str':
        if not annotation.has_mods():
            return [annotation.sequence[span[0]: span[1]] for span in spans]

        return [annotation.slice(span[0], span[1]).serialize() for span in spans]

    elif return_type == 'str-span':
        if not annotation.has_mods():
            return [(annotation.sequence[span[0]: span[1]], span) for span in spans]

        return [(annotation.slice(span[0], span[1]).serialize(), span) for span in spans]

    elif return_type == 'annotation-span':
        if not annotation.has_mods():
            return [(create_annotation(annotation.sequence[span[0]: span[1]]), span) for span in spans]

        return [(annotation.slice(span[0], span[1]), span) for span in spans]

    else:
        raise ValueError(f"Unsupported return type: {return_type}")


def get_left_semi_enzymatic_sequences(sequence: Union[str, ProFormaAnnotation],
                                      min_len: Optional[int] = None,
                                      max_len: Optional[int] = None,
                                      return_type: DigestReturnType = 'str') -> DIGEST_RETURN_TYPING:
    """
    Builds all left-hand semi-enzymatic subsequences derived from the input `sequence`.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: Union[str, ProFormaAnnotation]
    :param min_len: Minimum length for the subsequences (inclusive), defaults to [None]. If None, min_len will be
                    equal to 1.
    :type min_len: Optional[int]
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Optional[int]
    :param return_type: The type of the returned values.
    :type return_type: DigestReturnType

    :return: The left-hand semi-enzymatic subsequences. Return type is determined by the `return_type` parameter.
    :rtype: List[str] | List[ProFormaAnnotation] | List[Span] |List[(str, Span)] | List[(ProFormaAnnotation, Span)]

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
    elif isinstance(sequence, ProFormaAnnotation):
        annotation = sequence
    else:
        raise ValueError(f"Unsupported input type: {type(sequence)}")

    s = (0, len(annotation), 0)
    spans = build_left_semi_spans(span=s, min_len=min_len, max_len=max_len)

    return _return_digested_sequences(annotation, spans, return_type)


def get_right_semi_enzymatic_sequences(sequence: Union[str, ProFormaAnnotation],
                                       min_len: Optional[int] = None,
                                       max_len: Optional[int] = None,
                                       return_type: DigestReturnType = 'str') -> DIGEST_RETURN_TYPING:
    """
    Builds all right-hand semi-enzymatic subsequences derived from the input `sequence`.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: Union[str, ProFormaAnnotation]
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: Optional[int]
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Optional[int]
    :param return_type: The type of the returned values.
    :type return_type: DigestReturnType

    :return: The right-hand semi-enzymatic subsequences. Return type is determined by the `return_type` parameter.
    :rtype: List[str] | List[ProFormaAnnotation] | List[Span] |List[(str, Span)] | List[(ProFormaAnnotation, Span)]

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
    elif isinstance(sequence, ProFormaAnnotation):
        annotation = sequence
    else:
        raise ValueError(f"Unsupported input type: {type(sequence)}")

    s = (0, len(annotation), 0)
    spans = build_right_semi_spans(span=s, min_len=min_len, max_len=max_len)

    return _return_digested_sequences(annotation, spans, return_type)


def get_semi_enzymatic_sequences(sequence: Union[str, ProFormaAnnotation],
                                 min_len: Optional[int] = None,
                                 max_len: Optional[int] = None,
                                 return_type: DigestReturnType = 'str') -> DIGEST_RETURN_TYPING:
    """
    Builds allsemi-enzymatic sequences from the given input `sequence`.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: Union[str, ProFormaAnnotation]
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: Optional[int]
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up  to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Optional[int]
    :param return_type: The type of the returned values.
    :type return_type: DigestReturnType

    :return: Semi-enzymatic subsequences. Return type is determined by the `return_type` parameter.
    :rtype: List[str] | List[ProFormaAnnotation] | List[Span] |List[(str, Span)] | List[(ProFormaAnnotation, Span)]

    .. code-block:: python

        # Equivalent to build_left_semi_sequences + build_right_semi_sequences
        >>> res = get_left_semi_enzymatic_sequences('PEPTIDE') + get_right_semi_enzymatic_sequences('PEPTIDE')
        >>> get_semi_enzymatic_sequences('PEPTIDE') == res
        True

    """

    return get_left_semi_enzymatic_sequences(sequence=sequence, min_len=min_len, max_len=max_len,
                                             return_type=return_type) + \
        get_right_semi_enzymatic_sequences(sequence=sequence, min_len=min_len, max_len=max_len,
                                           return_type=return_type)


def get_non_enzymatic_sequences(sequence: Union[str, ProFormaAnnotation],
                                min_len: Optional[int] = None,
                                max_len: Optional[int] = None,
                                return_type: DigestReturnType = 'str') -> DIGEST_RETURN_TYPING:
    """
    Builds all non-enzymatic sequences from the given input `sequence`.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: Union[str, ProFormaAnnotation]
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: Optional[int]
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up  to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Optional[int]
    :param return_type: The type of the returned values.
    :type return_type: DigestReturnType

    :return: Non-enzymatic subsequences. Return type is determined by the `return_type` parameter.
    :rtype: List[str] | List[ProFormaAnnotation] | List[Span] |List[(str, Span)] | List[(ProFormaAnnotation, Span)]

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
    return _return_digested_sequences(annotation, spans, return_type)


def get_cleavage_sites(sequence: Union[str, ProFormaAnnotation], enzyme_regex: str) -> List[int]:
    """
    Return a list of positions where cleavage occurs in input `sequence` based on the provided enzyme regex.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: Union[str, ProFormaAnnotation]
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


def digest(sequence: Union[str, ProFormaAnnotation],
           enzyme_regex: Union[List[str], str],
           missed_cleavages: int = 0,
           semi: bool = False,
           min_len: Optional[int] = None,
           max_len: Optional[int] = None,
           return_type: DigestReturnType = 'str') -> DIGEST_RETURN_TYPING:
    """
    Returns a list of digested sequences derived from the input `sequence`.

    :param sequence: A sequence or ProFormaAnnotation.
    :type sequence: Union[str, ProFormaAnnotation]
    :param enzyme_regex: Regular expression or list of regular expressions representing enzyme's cleavage rules.
    :type enzyme_regex: Union[List[str], str]
    :param missed_cleavages: Maximum number of missed cleavages, defaults to [0].
    :type missed_cleavages: int
    :param semi: Whether to include semi-enzymatic peptides, defaults to [False].
    :type semi: bool
    :param min_len: Minimum length for the subsequences (inclusive), default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: Optional[int]
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up  to
                    1 - length of the `sequence`, defaults to [None].
    :type max_len: Optional[int]
    :param return_type: Whether to return the digested sequences as strings, ProFormaAnnotations, or Spans,
    defaults to ['str'].
    :type return_type: DigestReturnType

    :return: List of digested peptides. Return type is determined by the `return_type` parameter.
    :rtype: List[str] | List[ProFormaAnnotation] | List[Span] |List[(str, Span)] | List[(ProFormaAnnotation, Span)]

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

    if isinstance(enzyme_regex, str):
        enzyme_regex = [enzyme_regex]

    cleavage_sites = []
    for regex in enzyme_regex:
        cleavage_sites.extend(get_cleavage_sites(sequence=annotation, enzyme_regex=regex))

    spans = build_spans(max_index=len(annotation), enzyme_sites=cleavage_sites, missed_cleavages=missed_cleavages,
                        min_len=min_len, max_len=max_len, semi=semi)

    return _return_digested_sequences(annotation, spans, return_type)
