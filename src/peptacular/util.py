"""
Utils.py
"""

import warnings
from typing import Union, List, Tuple, Dict, Generator
import regex


def convert_type(val: str) -> Union[str, int, float]:
    """
    Convert the given value to the appropriate type.

    :param val: The value to convert.
    :type val: str

    :return: The converted value.
    :rtype: Union[str, int, float]

    .. code-block:: python

        >>> convert_type("1.234")
        1.234

        >>> convert_type("1")
        1

        >>> convert_type("+1")
        1

        >>> convert_type("-1")
        -1

        >>> convert_type("1.0")
        1.0

        >>> convert_type("abc")
        'abc'

    """

    if isinstance(val, (int, float)):
        return val

    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            return val


def merge_dicts(d1: Dict, d2: Dict) -> Dict:
    """
    Merge two dictionaries. And remove any keys with value 0.
    """
    d = {}
    for k, v in d1.items():
        d[k] = v
    for k, v in d2.items():
        d[k] = d.get(k, 0) + v

    # remove any keys with value 0
    d = {k: v for k, v in d.items() if v != 0}

    return d


def get_regex_match_indices(input_str: str, regex_str: Union[str, regex.Pattern], offset: int = 0) \
        -> Generator[int, None, None]:
    """
    Identify the starting indexes of occurrences of a given regex pattern within a string.

    :param input_str: The sequence in which to search.
    :type input_str: str
    :param regex_str: The regex pattern to search for.
    :type regex_str: str
    :param offset: An optional offset to add to each identified index. Default is 0.
    :type offset: int
    :return: A list of starting indexes where the regex pattern is found in the sequence.
    :rtype: List[int]

    .. code-block:: python

        >>> list(get_regex_match_indices("PEPTIDE", "P"))
        [1, 3]

        >>> list(get_regex_match_indices("PEPTIDE", regex.compile("E")))
        [2, 7]

        >>> list(get_regex_match_indices("PEPTIDE", 'E'))
        [2, 7]

        # More complex regex
        >>> list(get_regex_match_indices("PEPTIDE", "P[ST]"))
        [3]

        >>> list(get_regex_match_indices("PPPP", "PP"))
        [1, 2, 3]

        >>> list(get_regex_match_indices("PEPCTIDE", "(?=C)"))
        [3]

        >>> list(get_regex_match_indices("PEPCTIDE", "C"))
        [4]

        >>> list(get_regex_match_indices("PEPCTIDE", "(?<=C)"))
        [4]

    """

    if not isinstance(regex_str, regex.Pattern):
        regex_pattern = regex.compile(regex_str)
    else:
        regex_pattern = regex_str

    for match in regex_pattern.finditer(input_str, overlapped=True):
        if match.start() != match.end():
            warnings.warn("The regex pattern has a match with a none zero length. Using start index + 1 for the match.")
            yield match.start() + offset + 1
        else:
            yield match.start() + offset


def get_regex_match_range(input_str: str, regex_str: Union[str, regex.Pattern], offset: int = 0) \
        -> List[Tuple[int, int]]:
    """
    Identify the starting indexes of occurrences of a given regex pattern within a string.

    :param input_str: The sequence in which to search.
    :type input_str: str
    :param regex_str: The regex pattern to search for.
    :type regex_str: str
    :param offset: An optional offset to add to each identified index. Default is 0.
    :type offset: int
    :return: A list of starting indexes where the regex pattern is found in the sequence.
    :rtype: List[int]

    .. code-block:: python

        >>> get_regex_match_range("PEPTIDE", "P")
        [(0, 1), (2, 3)]

        >>> get_regex_match_range("PEPTIDE", regex.compile("E"))
        [(1, 2), (6, 7)]

        # More complex regex
        >>> get_regex_match_range("PEPTIDE", "P[ST]")
        [(2, 4)]

        >>> get_regex_match_range("PPPP", "PP")
        [(0, 2), (1, 3), (2, 4)]

    """

    # if regex is compiled
    if isinstance(regex_str, regex.Pattern):
        return [(i.start() + offset, i.end() + offset) for i in regex_str.finditer(input_str)]

    return [(match.start() + offset, match.end() + offset) for match in regex.finditer(regex_str, input_str,
                                                                                       overlapped=True)]


def _validate_span(span: Tuple[int, int, int]) -> None:
    """
    Validates if a given span is valid.

    :param span: A tuple representing the span.
    :type span: Tuple[int, int, int]
    :raises ValueError: If the span is not valid.

    .. code-block:: python

        >>> _validate_span((0, 5, 0))  # No error raised

        >>> _validate_span((0, 0, 0))  # No error raised

        >>> _validate_span((5, 0, 0))  # Raises ValueError
        Traceback (most recent call last):
        ValueError: Start of span: 5, should be less than or equal to end of span: 0.

        >>> _validate_span((-1, 0, 0))  # Raises ValueError
        Traceback (most recent call last):
        ValueError: Start of span should be non-negative, got -1.

        >>> _validate_span((0, -1, 0))  # Raises ValueError
        Traceback (most recent call last):
        ValueError: End of span should be non-negative, got -1.

    """

    start, end, _ = span
    if start < 0:
        raise ValueError(f'Start of span should be non-negative, got {start}.')
    if end < 0:
        raise ValueError(f'End of span should be non-negative, got {end}.')
    if start > end:
        raise ValueError(f'Start of span: {start}, should be less than or equal to end of span: {end}.')
