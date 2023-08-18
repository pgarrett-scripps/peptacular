"""
Used in digest.py, and exclusively work with unmodified sequences!
"""

from functools import wraps
from typing import Tuple, List, Optional, Callable
from itertools import groupby


# TODO: Remove wrapper function? Its confusing and hurts readability

def span_processing(func: Callable):
    """
    A decorator to enforce constraints on the span and min_len and max_len arguments of decorated functions.

    :param func: The function to be decorated.
    :type func: Callable
    :return: A wrapped function with constraints enforced on its arguments.
    :rtype: Callable
    """

    @wraps(func)
    def wrapper(span: Tuple[int, int, int], min_len: Optional[int] = None, max_len: Optional[int] = None):
        assert isinstance(span, Tuple) and len(span) == 3, 'span should be a tuple of length 3.'
        assert all(isinstance(arg, (int, type(None))) for arg in
                   (min_len, max_len)), 'min_len and max_len should be None or an integer.'
        _validate_span(span)

        if min_len is None:
            min_len = 1
        if max_len is None:
            max_len = span[1] - span[0]

        return func(span, min_len, max_len)

    return wrapper


def _validate_span(span: Tuple[int, int, int]) -> None:
    """
    Validates if a given span is valid.

    :param span: A tuple representing the span.
    :type span: Tuple[int, int, int]
    :raises ValueError: If the span is not valid.

    :Example:

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


@span_processing
def build_non_enzymatic_spans(span: Tuple[int, int, int], min_len: int, max_len: int) -> List[Tuple[int, int, int]]:
    """
    Generates and returns all possible sub-spans of the given span.

    :param span: A tuple representing the original span.
    :type span: Tuple[int, int, int]
    :param min_len: The minimum length of sub-spans to be generated.
    :type min_len: int
    :param max_len: The maximum length of sub-spans to be generated.
    :type max_len: int
    :return: A list of all possible sub-spans as tuples.
    :rtype: List[Tuple[int, int, int]]

    :Example:

    .. code-block:: python

        >>> build_non_enzymatic_spans((0, 3, 0), 1, 2)
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

    """

    start, end, _ = span
    return [(i, j, 0) for i in range(start, end + 1) for j in range(i + min_len, min(end + 1, i + max_len + 1))]


@span_processing
def build_left_semi_spans(span: Tuple[int, int, int], min_len: int, max_len: int) -> List[Tuple[int, int, int]]:
    """
    Generates and returns all possible sub-spans of the given span starting from the left.

    :param span: A tuple representing the original span.
    :type span: Tuple[int, int, int]
    :param min_len: The minimum length of sub-spans to be generated.
    :type min_len: int
    :param max_len: The maximum length of sub-spans to be generated.
    :type max_len: int
    :return: A list of all possible left sub-spans as tuples.
    :rtype: List[Tuple[int, int, int]]
    """

    if min_len > max_len:
        return []

    start, end, value = span
    new_end = min(start + max_len, end - 1)
    return [(start, i, value) for i in range(new_end, start - 1, -1) if i - start >= min_len]


@span_processing
def build_right_semi_spans(span: Tuple[int, int, int], min_len: int, max_len: int) -> List[Tuple[int, int, int]]:
    """
    Generates and returns all possible sub-spans of the given span starting from the right.

    :param span: A tuple representing the original span.
    :type span: Tuple[int, int, int]
    :param min_len: The minimum length of sub-spans to be generated.
    :type min_len: int
    :param max_len: The maximum length of sub-spans to be generated.
    :type max_len: int
    :return: A list of all possible right sub-spans as tuples.
    :rtype: List[Tuple[int, int, int]]
    """

    if min_len > max_len:
        return []

    start, end, value = span
    new_start = max(start + 1, end - max_len)
    return [(i, end, value) for i in range(new_start, end + 1) if end - i >= min_len]


def span_to_sequence(sequence: str, span: Tuple[int, int, int]) -> str:
    """
    Extracts a subsequence from the input sequence based on the provided span.

    :param sequence: The original sequence.
    :type sequence: str
    :param span: A tuple representing the span of the subsequence to be extracted.
    :type span: Tuple[int, int, int]
    :return: The subsequence of the input sequence defined by the span.
    :rtype: str
    """

    _validate_span(span)
    return sequence[span[0]:span[1]]


def get_enzymatic_spans(max_index: int, enzyme_sites: List[int], missed_cleavages: int,
                        min_len: int = None, max_len: int = None) -> List[Tuple[int, int, int]]:
    """
    Computes enzymatic spans for the given enzyme sites and missed cleavages.

    :param max_index: The max index of the span.
    :type max_index: int
    :param enzyme_sites: The list of indices representing enzyme sites.
    :type enzyme_sites: List[int]
    :param missed_cleavages: The number of allowed missed cleavages.
    :type missed_cleavages: int
    :param min_len: The minimum length of an enzymatic span.
    :type min_len: int, optional
    :param max_len: The maximum length of an enzymatic span.
    :type max_len: int, optional
    :return: A list of tuples representing enzymatic spans.
    :rtype: List[Tuple[int, int, int]]
    """

    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = max_index

    enzyme_sites = set(enzyme_sites)
    if 0 not in enzyme_sites:
        enzyme_sites.add(0)

    if max_index not in enzyme_sites:
        enzyme_sites.add(max_index)

    enzyme_sites = sorted(list(enzyme_sites))

    spans = []
    for i, start_site in enumerate(enzyme_sites):
        spans.extend(
            [(start_site, end_site, j) for j, end_site in enumerate(enzyme_sites[i + 1:i + missed_cleavages + 2])])

    # Filter spans based on length
    spans = [span for span in spans if min_len <= span[1] - span[0] <= max_len]

    return spans


def _get_all_left_semi_spans(spans: List[Tuple[int, int, int]], min_len: int, max_len: int) -> \
        List[Tuple[int, int, int]]:
    """
    Get all left semi-spans from the given list of spans that have a length within the specified range.

    :param spans: A list of tuples, where each tuple represents a span with three integers (start, end, value).
    :param min_len: The minimum length of the left semi-spans.
    :param max_len: The maximum length of the left semi-spans.
    :return: A list of tuples representing all left semi-spans that are within the specified length range.

    This function groups the spans by the start position, for each group, it checks every span to see if its length
    is at least the minimum length. If it is, it calculates the new maximum length, which is the smallest of the max_len
    and the span's length. Then it adds the left semi-spans to the list of semi-spans.
    """
    semi_spans = []
    spans = sorted(spans, key=lambda x: (x[0], -x[2]))
    for _, group in groupby(spans, key=lambda x: x[0]):
        group = list(group)
        for i, span in enumerate(group):
            span_len = span[1] - span[0]

            if span_len <= min_len:
                break

            new_max_len = min(max_len, span_len - 1)
            if i == len(group) - 1:
                semi_spans.extend(build_left_semi_spans(span, min_len, new_max_len))
            else:
                next_span = group[i + 1]
                next_span_len = next_span[1] - next_span[0]
                new_min = max(min_len, next_span_len + 1)
                semi_spans.extend(build_left_semi_spans(span, new_min, new_max_len))

    return semi_spans


def _get_all_right_semi_spans(spans: List[Tuple[int, int, int]], min_len: int, max_len: int) -> \
        List[Tuple[int, int, int]]:
    """
    Get all right semi-spans from the given list of spans that have a length within the specified range.

    :param spans: A list of tuples, where each tuple represents a span with three integers (start, end, value).
    :param min_len: The minimum length of the right semi-spans.
    :param max_len: The maximum length of the right semi-spans.
    :return: A list of tuples representing all right semi-spans that are within the specified length range.

    This function groups the spans by the start position, for each group, it checks every span to see if its length
    is at least the minimum length. If it is, it calculates the new maximum length, which is the smallest of the max_len
    and the span's length. Then it adds the right semi-spans to the list of semi-spans.
    """

    semi_spans = []
    spans = sorted(spans, key=lambda x: (x[1], -x[2]))
    for _, group in groupby(spans, key=lambda x: x[1]):
        group = list(group)
        for i, span in enumerate(group):
            span_len = span[1] - span[0]

            if span_len < min_len:
                break

            new_max_len = min(max_len, span_len - 1)
            if i == len(group) - 1:
                semi_spans.extend(build_right_semi_spans(span, min_len, new_max_len))
            else:
                next_span = group[i + 1]
                next_span_len = next_span[1] - next_span[0]
                new_min = max(min_len, next_span_len + 1)
                semi_spans.extend(build_right_semi_spans(span, new_min, new_max_len))

    return semi_spans


def get_semi_spans(spans: List[Tuple[int, int, int]], min_len: int, max_len: int) -> \
        List[Tuple[int, int, int]]:
    """
    Computes semi spans for a given list of spans based on length criteria.

    :param spans: The list of spans.
    :type spans: List[Tuple[int, int, int]]
    :param min_len: The minimum length of a semi span.
    :type min_len: int
    :param max_len: The maximum length of a semi span.
    :type max_len: int
    :return: A list of tuples representing semi spans.
    :rtype: List[Tuple[int, int, int]]
    """

    semi_spans = []

    spans = sorted(spans, key=lambda x: (x[0], -x[2]))
    for _, group in groupby(spans, key=lambda x: x[0]):
        semi_spans.extend(_get_all_left_semi_spans(list(group), min_len, max_len))

    spans = sorted(spans, key=lambda x: (x[1], -x[2]))
    for _, group in groupby(spans, key=lambda x: x[1]):
        semi_spans.extend(_get_all_right_semi_spans(list(group), min_len, max_len))

    return semi_spans


def build_spans(max_index: int, enzyme_sites: List[int], missed_cleavages: int, min_len: int,
                max_len: int, semi: bool) -> List[Tuple[int, int, int]]:
    """
    Builds spans for a given sequence based on enzyme sites, missed cleavages, and length criteria.

    :param max_index: The max index of the span.
    :type max_index: int
    :param enzyme_sites: The list of indices representing enzyme sites.
    :type enzyme_sites: List[int]
    :param missed_cleavages: The number of allowed missed cleavages.
    :type missed_cleavages: int
    :param min_len: The minimum length of a span.
    :type min_len: int
    :param max_len: The maximum length of a span.
    :type max_len: int
    :param semi: Whether to compute and include semi spans.
    :type semi: bool
    :return: A list of tuples representing spans.
    :rtype: List[Tuple[int, int, int]]
    """

    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = max_index

    if semi:
        spans = get_enzymatic_spans(max_index, enzyme_sites, missed_cleavages, None, None)
        semi_spans = get_semi_spans(spans, min_len, max_len)

        # filter spans by len
        spans = [span for span in spans if max_len >= span[1] - span[0] >= min_len] + semi_spans
    else:
        spans = get_enzymatic_spans(max_index, enzyme_sites, missed_cleavages, min_len, max_len)
    return spans
