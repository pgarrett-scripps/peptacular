from typing import Union, List
import regex as reg

from peptacular._spans import build_left_semi_spans, build_right_semi_spans, span_to_sequence, \
    build_non_enzymatic_spans, build_spans
from peptacular.constants import PROTEASES


def _build_left_semi_sequences(sequence: str,
                               min_len: Union[int, None] = None, max_len: Union[int, None] = None) -> List[str]:
    """
    Returns left subsequences of the amino acid sequence `sequence` with lengths between `min_len` and `max_len`.

    :param sequence: Amino acid sequence to generate subsequences from.
    :type sequence: str
    :param min_len: Minimum length of the subsequence. If None, defaults to 1.
    :type min_len: Union[int, None]
    :param max_len: Maximum length of the subsequence. If None, defaults to the length of the sequence.
    :type max_len: Union[int, None]

    :return: List of left subsequences of the input sequence.
    :rtype: List[str]

    :Example:

    .. code-block:: python

        >>> _build_left_semi_sequences('PEPTIDE')
        ['PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        >>> _build_left_semi_sequences('P')
        []

        >>> _build_left_semi_sequences('')
        []

    """
    spans = build_left_semi_spans((0, len(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def _build_right_semi_sequences(sequence: str,
                                min_len: Union[int, None] = None, max_len: Union[int, None] = None) -> List[str]:
    """
    Returns right subsequences of the amino acid sequence `sequence` with lengths between `min_len` and `max_len`.

    :param sequence: Amino acid sequence to generate subsequences from.
    :type sequence: str
    :param min_len: Minimum length of the subsequence. If None, defaults to 1.
    :type min_len: Union[int, None]
    :param max_len: Maximum length of the subsequence. If None, defaults to the length of the sequence.
    :type max_len: Union[int, None]

    :return: List of right subsequences of the input sequence.
    :rtype: List[str]

    :Example:

    .. code-block:: python

        >>> _build_right_semi_sequences('PEPTIDE')
        ['EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> _build_right_semi_sequences('P')
        []

        >>> _build_right_semi_sequences('')
        []

    """

    spans = build_right_semi_spans((0, len(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def build_semi_sequences(sequence: str,
                         min_len: Union[int, None] = None, max_len: Union[int, None] = None) -> List[str]:
    """
    Generate a list of all semi-enzymatic sequences from the given amino acid sequence.

    :param sequence: Amino acid sequence to be processed.
    :type sequence: str
    :param min_len: Minimum length of the semi-enzymatic sequence. Default is None.
    :type min_len: Union[int, None]
    :param max_len: Maximum length of the semi-enzymatic sequence. Default is None.
    :type max_len: Union[int, None]

    :return: A list of semi-enzymatic sequences within the defined length bounds.
    :rtype: List[str]

    :Example:

    .. code-block:: python

        >>> build_semi_sequences('PEPTIDE')
        ['PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P', 'EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> build_semi_sequences('P')
        []

        >>> build_semi_sequences('')
        []

    """

    return _build_left_semi_sequences(sequence, min_len, max_len) + _build_right_semi_sequences(sequence, min_len, max_len)


def build_non_enzymatic_sequences(sequence: str,
                                  min_len: Union[int, None] = None, max_len: Union[int, None] = None) -> List[str]:
    """
    Generate a list of all non-enzymatic amino acid sequences from a given string sequence.

    :param sequence: Amino acid sequence to be processed.
    :type sequence: str
    :param min_len: Minimum length of the non-enzymatic sequence. Default is None.
    :type min_len: Union[int, None]
    :param max_len: Maximum length of the non-enzymatic sequence. Default is None.
    :type max_len: Union[int, None]

    :return: A list of non-enzymatic sequences within the defined length bounds.
    :rtype: List[str]

    :Example:

    .. code-block:: python

        >>> build_non_enzymatic_sequences('PEP')
        ['P', 'PE', 'PEP', 'E', 'EP', 'P']

        >>> build_non_enzymatic_sequences('P')
        ['P']

        >>> build_non_enzymatic_sequences('')
        []

    """
    spans = build_non_enzymatic_spans((0, len(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def build_enzymatic_sequences(sequence: str, enzyme_regex: str, missed_cleavages: int, semi: bool = False,
                              min_len: Union[int, None] = None, max_len: Union[int, None] = None) -> List[str]:
    """
    Generate a list of all enzymatic amino acid sequences based on the provided enzyme_regex.

    :param sequence: Amino acid sequence to be digested.
    :type sequence: str
    :param enzyme_regex: Regular expression or key in PROTEASES dictionary defining the enzyme's cleavage rule.
    :type enzyme_regex: str
    :param missed_cleavages: Maximum number of missed cleavages.
    :type missed_cleavages: int
    :param semi: Whether to include semi-enzymatic peptides.
    :type semi: bool
    :param min_len: Minimum length of the enzymatic sequence. Default is None.
    :type min_len: int, optional
    :param max_len: Maximum length of the enzymatic sequence. Default is None.
    :type max_len: int, optional

    :return: List of enzymatic sequences satisfying enzyme cleavage rules and within the defined length bounds.
    :rtype: List[str]

    :Example:

    .. code-block:: python

        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, max_len=5)
        ['TIDER', 'TIDEK', 'TIDE']

        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, min_len=6)
        ['TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEKTIDE']
    """

    cleavage_sites = identify_cleavage_sites(sequence, enzyme_regex)
    spans = build_spans(len(sequence), cleavage_sites, missed_cleavages, min_len, max_len, semi)
    return [span_to_sequence(sequence, span) for span in spans]


def digest_sequence(sequence: str, enzyme_regex: Union[List[str], str], missed_cleavages: int,
                    min_len: Union[int, None] = None, max_len: Union[int, None] = None, semi: bool = False) -> List[str]:
    """
    Digests a given amino acid sequence using specified enzyme rules, returning a list of resulting peptides.

    :param sequence: Amino acid sequence to be digested.
    :type sequence: str
    :param enzyme_regex: Regular expression or list of regular expressions representing enzyme's cleavage rules.
    :type enzyme_regex: Union[List[str], str]
    :param missed_cleavages: Number of allowed missed cleavages.
    :type missed_cleavages: int
    :param min_len: Minimum length of the returned peptides. Default is None.
    :type min_len: int, optional
    :param max_len: Maximum length of the returned peptides. Default is None.
    :type max_len: int, optional
    :param semi: Whether to include semi-digested peptides.
    :type semi: bool

    :return: List of digested peptides.
    :rtype: List[str]

    :Example:

    .. code-block:: python

        >>> digest_sequence(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        >>> digest_sequence(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        >>> digest_sequence(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, max_len=5)
        ['TIDER', 'TIDEK', 'TIDE']

        >>> digest_sequence(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, min_len=6)
        ['TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEKTIDE']

        >>> digest_sequence(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=1, min_len=8, semi=True)
        ['TIDERTIDEK', 'TIDEKTIDE', 'TIDERTIDE', 'TIDERTID', 'TIDEKTID', 'IDERTIDEK', 'DERTIDEK', 'IDEKTIDE']

    """

    if isinstance(enzyme_regex, str):
        enzyme_regex = [enzyme_regex]

    cleavage_sites = []
    for regex in enzyme_regex:
        cleavage_sites.extend(identify_cleavage_sites(sequence, regex))
    cleavage_sites.sort()

    spans = build_spans(len(sequence), cleavage_sites, missed_cleavages, min_len, max_len, semi)
    sequences = [span_to_sequence(sequence, span) for span in spans]

    return sequences


def identify_cleavage_sites(sequence: str, enzyme_regex: str) -> List[int]:
    """
    Identifies cleavage sites in an amino acid sequence based on the provided enzyme_regex pattern.
    If enzyme_regex is a key in PROTEASES, the corresponding regex pattern will be used.
    Cleavage sites are identified as positions after matching residues.

    :param sequence: Amino acid sequence to be processed.
    :type sequence: str
    :param enzyme_regex: Regular expression or key in PROTEASES dictionary defining enzyme's cleavage rule.
    :type enzyme_regex: str

    :return: List of positions where cleavage occurs in the sequence.
    :rtype: List[int]

    :Example:

    .. code-block:: python

        >>> identify_cleavage_sites(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P')
        [5, 10]

        >>> identify_cleavage_sites(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])')
        [5, 10]

        >>> identify_cleavage_sites(sequence='TIDEPTIDEPTIDE', enzyme_regex='trypsin/P')
        []

        >>> identify_cleavage_sites(sequence='KPEPTIDEK', enzyme_regex='lys-n')
        [0, 8]

        >>> identify_cleavage_sites(sequence='KPEPTIDEK', enzyme_regex='lys-c')
        [1, 9]

    """

    if enzyme_regex in PROTEASES:
        enzyme_regex = PROTEASES[enzyme_regex]

    enzyme_sites = []
    for site in reg.finditer(enzyme_regex, 'x' + sequence + 'x', overlapped=True):
        enzyme_sites.append(site.span(0))
    return [site[0] for site in enzyme_sites]
