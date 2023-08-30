from typing import Union, List

from peptacular.spans import build_left_semi_spans, build_right_semi_spans, build_non_enzymatic_spans, build_spans
from peptacular.constants import PROTEASES
from peptacular.sequence import strip_modifications, calculate_sequence_length, span_to_sequence
from peptacular.util import identify_regex_indexes


def build_left_semi_sequences(sequence: str, min_len: int = 1, max_len: int = None) -> List[str]:
    """
    Builds all left-hand semi-enzymatic subsequences derived from the input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param min_len: Minimum length for the subsequences (inclusive), defaults to [1].
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to the
                    length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]

    :return: Left-hand semi-enzymatic subsequences.
    :rtype: List[str]

    .. code-block:: python

        # Generates all left-hand semi enzymatic sequences (Returned values does not include input sequence)
        >>> build_left_semi_sequences('PEPTIDE')
        ['PEPTID', 'PEPTI', 'PEPT', 'PEP', 'PE', 'P']

        # For single-letter or empty sequences, the function returns an empty list:
        >>> build_left_semi_sequences('P')
        []
        >>> build_left_semi_sequences('')
        []

        # Modifications are preserved:
        >>> build_left_semi_sequences('[1]P(2)EPTIDE')
        ['[1]P(2)EPTID', '[1]P(2)EPTI', '[1]P(2)EPT', '[1]P(2)EP', '[1]P(2)E', '[1]P(2)']

    """

    spans = build_left_semi_spans((0, calculate_sequence_length(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def build_right_semi_sequences(sequence: str, min_len: int = 1, max_len: int = None) -> List[str]:
    """
    Builds all right-hand semi-enzymatic subsequences derived from the input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param min_len: Minimum length for the subsequences (inclusive), defaults to [1].
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to the
                    length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]

    :return: Right-hand semi-enzymatic subsequences
    :rtype: List[str]

    .. code-block:: python

        # Generates all right-hand semi enzymatic sequences (Returned values does not include input sequence)
        >>> build_right_semi_sequences('PEPTIDE')
        ['EPTIDE', 'PTIDE', 'TIDE', 'IDE', 'DE', 'E']

        # For single-letter or empty sequences, the function returns an empty list:
        >>> build_right_semi_sequences('P')
        []
        >>> build_right_semi_sequences('')
        []

        # Modifications are preserved:
        >>> build_right_semi_sequences('PEPTIDE(1)[2]')
        ['EPTIDE(1)[2]', 'PTIDE(1)[2]', 'TIDE(1)[2]', 'IDE(1)[2]', 'DE(1)[2]', 'E(1)[2]']

    """

    spans = build_right_semi_spans((0, calculate_sequence_length(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def build_semi_sequences(sequence: str,  min_len: int = 1, max_len: int = None) -> List[str]:
    """
    Builds allsemi-enzymatic sequences from the given input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param min_len: Minimum length for the subsequences (inclusive), defaults to [1].
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to the
                    length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]

    :return: Semi-enzymatic subsequences.
    :rtype: List[str]

    .. code-block:: python

        # Equivalent to build_left_semi_sequences + build_right_semi_sequences
        >>> res = build_left_semi_sequences('PEPTIDE') + build_right_semi_sequences('PEPTIDE')
        >>> build_semi_sequences('PEPTIDE') == res
        True

    """

    return build_left_semi_sequences(sequence, min_len, max_len) + build_right_semi_sequences(sequence, min_len,
                                                                                              max_len)


def build_non_enzymatic_sequences(sequence: str, min_len: int = 1, max_len: int = None) -> List[str]:
    """
    Builds all non-enzymatic sequences from the given input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param min_len: Minimum length for the subsequences (inclusive), defaults to [1].
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to the
                    length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]

    :return: Non-enzymatic subsequences
    :rtype: List[str]

    .. code-block:: python

        # Generates non-enzymatic sequences (Returned values does not include input sequence):
        >>> build_non_enzymatic_sequences('PEP')
        ['P', 'PE', 'E', 'EP', 'P']

        # For single-letter or empty sequences, the function returns an empty list:
        >>> build_non_enzymatic_sequences('P')
        []
        >>> build_non_enzymatic_sequences('')
        []

         # Sequences with modifications are processed preserving those modifications:
        >>> build_non_enzymatic_sequences('[Acetyl]P(1.0)EP(1.0)[Amide]')
        ['[Acetyl]P(1.0)', '[Acetyl]P(1.0)E', 'E', 'EP(1.0)[Amide]', 'P(1.0)[Amide]']

    """

    spans = build_non_enzymatic_spans((0, calculate_sequence_length(sequence), 0), min_len, max_len)
    return [span_to_sequence(sequence, span) for span in spans]


def build_enzymatic_sequences(sequence: str, enzyme_regex: str, missed_cleavages: int, semi: bool = False,
                              min_len: int = 1, max_len: int = None) -> List[str]:
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
    :param min_len: Minimum length for the subsequences (inclusive), defaults to [1].
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to the
                    length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]

    :return: Enzymatic subsequences.
    :rtype: List[str]

    .. code-block:: python

        # Can specify the name of a protease in the PROTEASES dictionary:
        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        # Or specify a regular expression:
        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2)
        ['TIDER', 'TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEK', 'TIDEKTIDE', 'TIDE']

        # Filter sequences by max length:
        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, max_len=5)
        ['TIDER', 'TIDEK', 'TIDE']

        # Filter sequences by min length:
        >>> build_enzymatic_sequences(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])', missed_cleavages=2, min_len=6)
        ['TIDERTIDEK', 'TIDERTIDEKTIDE', 'TIDEKTIDE']

        # Modifications are preserved:
        >>> build_enzymatic_sequences(sequence='[1]TIDERT(1.0)IDEKTIDE[2]', enzyme_regex='([KR])', missed_cleavages=2)
        ['[1]TIDER', '[1]TIDERT(1.0)IDEK', '[1]TIDERT(1.0)IDEKTIDE[2]', 'T(1.0)IDEK', 'T(1.0)IDEKTIDE[2]', 'TIDE[2]']

    """

    cleavage_sites = identify_cleavage_sites(sequence, enzyme_regex)
    spans = build_spans(calculate_sequence_length(sequence), cleavage_sites, missed_cleavages, min_len, max_len, semi)
    return [span_to_sequence(sequence, span) for span in spans]


def identify_cleavage_sites(sequence: str, enzyme_regex: str) -> List[int]:
    """
    Return a list of positions where cleavage occurs in input `sequence` based on the provided enzyme regex.

    Note: 'x' is appended to the start and end of sequence to ensure that the first and last residues are matched.
    Ensure that the regex pattern does not match 'x' if this is not desired.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param enzyme_regex: Regular expression or key in PROTEASES dictionary defining enzyme's cleavage rule.
    :type enzyme_regex: str

    :return: List of positions where cleavage occurs in the sequence.
    :rtype: List[int]

    .. code-block:: python

        # Can use a key in PROTEASES to specify the enzyme_regex:
        >>> identify_cleavage_sites(sequence='TIDERTIDEKTIDE', enzyme_regex='trypsin/P')
        [5, 10]

        # Or specify a regular expression:
        >>> identify_cleavage_sites(sequence='TIDERTIDEKTIDE', enzyme_regex='([KR])')
        [5, 10]

        # No cleavage sites are identified if the enzyme_regex does not match the sequence:
        >>> identify_cleavage_sites(sequence='TIDEPTIDEPTIDE', enzyme_regex='trypsin/P')
        []

        # If the protease cleaves at the N-terminus, the first position is included:
        >>> identify_cleavage_sites(sequence='KPEPTIDEK', enzyme_regex='lys-n')
        [0, 8]

        # Similarly, if the protease cleaves at the C-terminus, the last position is included:
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


def digest(sequence: str, enzyme_regex: Union[List[str], str], missed_cleavages: int,
           semi: bool = False, min_len: int = 1, max_len: int = None) -> List[str]:
    """
    Returns a list of digested sequences derived from the input `sequence`.

    :param sequence: The amino acid sequence, which can include modifications.
    :type sequence: str
    :param enzyme_regex: Regular expression or list of regular expressions representing enzyme's cleavage rules.
    :type enzyme_regex: Union[List[str], str]
    :param missed_cleavages: Maximum number of missed cleavages.
    :type missed_cleavages: int
    :param semi: Whether to include semi-enzymatic peptides, defaults to [False].
    :type semi: bool
    :param min_len: Minimum length for the subsequences (inclusive), defaults to [1].
    :type min_len: int
    :param max_len: Maximum length for the subsequences (inclusive). If None, the subsequences will go up to the
                    length of the `sequence`, defaults to [None].
    :type max_len: Union[int, None]

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
        >>> digest(sequence='TIDERTIDEK(1)TIDE[2]', enzyme_regex='([KR])', missed_cleavages=1, min_len=9, semi=True)
        ['TIDERTIDEK(1)', 'TIDEK(1)TIDE[2]', 'TIDERTIDE', 'IDERTIDEK(1)']
    """

    if isinstance(enzyme_regex, str):
        enzyme_regex = [enzyme_regex]

    cleavage_sites = []
    for regex in enzyme_regex:
        cleavage_sites.extend(identify_cleavage_sites(sequence, regex))
    cleavage_sites.sort()

    spans = build_spans(calculate_sequence_length(sequence), cleavage_sites, missed_cleavages, min_len, max_len, semi)
    sequences = [span_to_sequence(sequence, span) for span in spans]

    return sequences
