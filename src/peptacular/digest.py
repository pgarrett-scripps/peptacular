from typing import Union, List

from peptacular._spans import build_left_semi_spans, build_right_semi_spans, build_non_enzymatic_spans, build_spans
from peptacular.constants import PROTEASES
from peptacular.sequence import strip_modifications, calculate_sequence_length, span_to_sequence
from peptacular.util import identify_regex_indexes


def build_left_semi_sequences(sequence: str, min_len: int = None, max_len: int = None) -> List[str]:
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

    .. code-block:: python

        # Works with unmodified sequences
        >>> build_left_semi_sequences('PEPTIDE')
        ['PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        >>> build_left_semi_sequences('P')
        []

        >>> build_left_semi_sequences('')
        []

        # but will also preserve modifications, including terminal modifications
        >>> build_left_semi_sequences('[1]P(2)EPTIDE')
        ['[1]P(2)EPTID', '[1]P(2)EPTI', '[1]P(2)EPT', '[1]P(2)EP', '[1]P(2)E', '[1]P(2)']

        >>> build_left_semi_sequences('PEPTIDE(1.0)[Amide]')
        ['PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

    """

    spans = build_left_semi_spans((0, calculate_sequence_length(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def build_right_semi_sequences(sequence: str, min_len: int = None, max_len: int = None) -> List[str]:
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

    .. code-block:: python

        # Works with unmodified sequences
        >>> build_right_semi_sequences('PEPTIDE')
        ['EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> build_right_semi_sequences('P')
        []

        >>> build_right_semi_sequences('')
        []

        # but will also preserve modifications, including terminal modifications
        >>> build_right_semi_sequences('[Acetyl]P(1.0)EPTIDE')
        ['EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        >>> build_right_semi_sequences('PEPTIDE(1)[2]')
        ['EPTIDE(1)[2]', 'PTIDE(1)[2]', 'TIDE(1)[2]', 'IDE(1)[2]', 'DE(1)[2]', 'E(1)[2]']

    """

    spans = build_right_semi_spans((0, calculate_sequence_length(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def build_semi_sequences(sequence: str,  min_len: int = None, max_len: int = None) -> List[str]:
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

    .. code-block:: python

        # Equivalent to build_left_semi_sequences + build_right_semi_sequences
        >>> res = build_left_semi_sequences('PEPTIDE') + build_right_semi_sequences('PEPTIDE')
        >>> build_semi_sequences('PEPTIDE') == res
        True

    """

    return build_left_semi_sequences(sequence, min_len, max_len) + build_right_semi_sequences(sequence, min_len,
                                                                                              max_len)


def build_non_enzymatic_sequences(sequence: str, min_len: int = None, max_len: int = None) -> List[str]:
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

    .. code-block:: python

        # Works with unmodified sequences
        >>> build_non_enzymatic_sequences('PEP')
        ['P', 'PE', 'PEP', 'E', 'EP', 'P']

        >>> build_non_enzymatic_sequences('P')
        ['P']

        >>> build_non_enzymatic_sequences('')
        []

        # but will also preserve modifications, including terminal modifications
        >>> build_non_enzymatic_sequences('[Acetyl]P(1.0)EP(1.0)[Amide]')
        ['[Acetyl]P(1.0)', '[Acetyl]P(1.0)E', '[Acetyl]P(1.0)EP(1.0)[Amide]', 'E', 'EP(1.0)[Amide]', 'P(1.0)[Amide]']

    """

    spans = build_non_enzymatic_spans((0, calculate_sequence_length(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def build_enzymatic_sequences(sequence: str, enzyme_regex: str, missed_cleavages: int, semi: bool = False,
                              min_len: int = None, max_len: int = None) -> List[str]:
    """
    Generate a list of all enzymatic amino acid sequences based on the provided enzyme regex. If enzyme_regex is a
    key in PROTEASES, the corresponding regex pattern will be used.


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

    .. code-block:: python

        Can specify the name of a protease in the PROTEASES dictionary:
        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        Or specify a regular expression:
        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        Can specify a maximum length:
        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, max_len=5)
        ['TIDER', 'TIDEK', 'TIDE']

        Can specify a minimum length:
        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, min_len=6)
        ['TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEKTIDE']

        # Will also work with modified sequences
        >>> build_enzymatic_sequences(sequence='[1]TIDERT(1.0)IDEKTIDE[2]', enzyme_regex='([KR])', missed_cleavages=2)
        ['[1]TIDER', '[1]TIDERT(1.0)IDEK', '[1]TIDERT(1.0)IDEKTIDE[2]', 'T(1.0)IDEK', 'T(1.0)IDEKTIDE[2]', 'TIDE[2]']

    """

    cleavage_sites = identify_cleavage_sites(sequence, enzyme_regex)
    spans = build_spans(calculate_sequence_length(sequence), cleavage_sites, missed_cleavages, min_len, max_len, semi)
    return [span_to_sequence(sequence, span) for span in spans]


def identify_cleavage_sites(sequence: str, enzyme_regex: str) -> List[int]:
    """
    Identifies cleavage sites in an amino acid sequence based on the provided enzyme regex pattern.
    If enzyme_regex is a key in PROTEASES, the corresponding regex pattern will be used.
    Cleavage sites are identified as positions after matching residues.

    Note: 'x' is appended to the start and end of sequence to ensure that the first and last residues are matched.
    Ensure that the regex pattern does not match 'x' if this is not desired.

    :param sequence: Amino acid sequence to be processed.
    :type sequence: str
    :param enzyme_regex: Regular expression or key in PROTEASES dictionary defining enzyme's cleavage rule.
    :type enzyme_regex: str

    :return: List of positions where cleavage occurs in the sequence.
    :rtype: List[int]

    .. code-block:: python

        Can use a key in PROTEASES to specify the enzyme_regex:
        >>> identify_cleavage_sites(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P')
        [5, 10]

        Or specify a regular expression:
        >>> identify_cleavage_sites(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])')
        [5, 10]

        No cleavage sites are identified if the enzyme_regex does not match the sequence:
        >>> identify_cleavage_sites(sequence='TIDEPTIDEPTIDE', enzyme_regex='trypsin/P')
        []

        If the protease cleaves at the N-terminus, the first position is included:
        >>> identify_cleavage_sites(sequence='KPEPTIDEK', enzyme_regex='lys-n')
        [0, 8]

        Similarly, if the protease cleaves at the C-terminus, the last position is included:
        >>> identify_cleavage_sites(sequence='KPEPTIDEK', enzyme_regex='lys-c')
        [1, 9]

        # Will also work with modified sequences
        >>> identify_cleavage_sites(sequence='[Acetyl]TIDERT(1.0)IDEKTIDE[Amide]', enzyme_regex='trypsin/P')
        [5, 10]

    """

    stripped_sequence = strip_modifications(sequence)
    if stripped_sequence != sequence:
        sequence = stripped_sequence

    if enzyme_regex in PROTEASES:
        enzyme_regex = PROTEASES[enzyme_regex]

    return identify_regex_indexes('x' + sequence + 'x', enzyme_regex)


def digest_sequence(sequence: str, enzyme_regex: Union[List[str], str], missed_cleavages: int,
                    min_len: int = None, max_len: int = None, semi: bool = False) -> List[str]:
    """
    Digests a given amino acid sequence using specified enzyme regex pattern. If enzyme_regex is a key in
    PROTEASES, the corresponding regex pattern will be used.


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