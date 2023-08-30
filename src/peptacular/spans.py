"""
Used in digest.py, and exclusively work with unmodified sequences!
"""

from typing import Tuple, List
from itertools import groupby


def build_non_enzymatic_spans(span: Tuple[int, int, int], min_len: int = 1, max_len: int = None) \
        -> List[Tuple[int, int, int]]:
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

    .. code-block:: python

        >>> build_non_enzymatic_spans((0, 3, 0), 1, 2)
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

        >>> build_non_enzymatic_spans((0, 3, 0))
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

        >>> build_non_enzymatic_spans((0, 3, 0), 1, 10)
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

    """

    if min_len is None:
        min_len = 1

    max_span = span[1] - span[0] - 1
    if max_len is None:
        max_len = max_span
    max_len = min(max_len, max_span)

    start, end, _ = span
    return [(i, j, 0) for i in range(start, end + 1) for j in range(i + min_len, min(end + 1, i + max_len + 1))]


def build_left_semi_spans(span: Tuple[int, int, int], min_len: int = 1, max_len: int = None) \
        -> List[Tuple[int, int, int]]:
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

    .. code-block:: python

        >>> build_left_semi_spans((0, 3, 0))
        [(0, 2, 0), (0, 1, 0)]

        # Keeps the value of the original span
        >>> build_left_semi_spans((0, 3, 2))
        [(0, 2, 2), (0, 1, 2)]

        # Set min and max length
        >>> build_left_semi_spans((0, 3, 0), 1, 1)
        [(0, 1, 0)]
        >>> build_left_semi_spans((0, 3, 0), 1, 10)
        [(0, 2, 0), (0, 1, 0)]

    """
    if min_len is None:
        min_len = 1

    if max_len is None:
        max_len = span[1] - span[0]

    start, end, value = span
    new_end = min(start + max_len, end - 1)
    return [(start, i, value) for i in range(new_end, start - 1, -1) if i - start >= min_len]


def build_right_semi_spans(span: Tuple[int, int, int], min_len: int = 1, max_len: int = None) \
        -> List[Tuple[int, int, int]]:
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

    .. code-block:: python

        >>> build_right_semi_spans((0, 3, 0))
        [(1, 3, 0), (2, 3, 0)]

        # Keeps the value of the original span
        >>> build_right_semi_spans((0, 3, 2))
        [(1, 3, 2), (2, 3, 2)]

        # Set min and max length
        >>> build_right_semi_spans((0, 3, 0), 1, 1)
        [(2, 3, 0)]
        >>> build_right_semi_spans((0, 3, 0), 1, 10)
        [(1, 3, 0), (2, 3, 0)]

    """

    if min_len is None:
        min_len = 1

    if max_len is None:
        max_len = span[1] - span[0]

    start, end, value = span
    new_start = max(start + 1, end - max_len)
    return [(i, end, value) for i in range(new_start, end + 1) if end - i >= min_len]


def build_enzymatic_spans(max_index: int, enzyme_sites: List[int], missed_cleavages: int,
                          min_len: int = 1, max_len: int = None) -> List[Tuple[int, int, int]]:
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

    .. code-block:: python

        >>> build_enzymatic_spans(5, [0, 3, 5], 1)
        [(0, 3, 0), (0, 5, 1), (3, 5, 0)]

        # Set min length
        >>> build_enzymatic_spans(5, [0, 3, 5], 1, min_len=5)
        [(0, 5, 1)]

        # Set max length
        >>> build_enzymatic_spans(5, [0, 3, 5], 1, max_len=3)
        [(0, 3, 0), (3, 5, 0)]

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


def _grouped_left_semi_span_builder(spans: List[Tuple[int, int, int]], min_len: int, max_len: int) -> \
        List[Tuple[int, int, int]]:
    """
    Get all left semi-spans from the given list of spans that have a length within the specified range.

    This function groups the spans by the start position, for each group, it checks every span to see if its length
    is at least the minimum length. If it is, it calculates the new maximum length, which is the smallest of the max_len
    and the span's length. Then it adds the left semi-spans to the list of semi-spans.

    :param spans: A list of tuples, where each tuple represents a span with three integers (start, end, value).
    :type spans: List[Tuple[int, int, int]]
    :param min_len: The minimum length of the left semi-spans.
    :type min_len: int
    :param max_len: The maximum length of the left semi-spans.
    :type max_len: int

    :return: A list of tuples representing all left semi-spans that are within the specified length range.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        >>> _grouped_left_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], 1, 5)
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0)]

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


def _grouped_right_semi_span_builder(spans: List[Tuple[int, int, int]], min_len: int, max_len: int) -> \
        List[Tuple[int, int, int]]:
    """
    Get all right semi-spans from the given list of spans that have a length within the specified range.

    This function groups the spans by the start position, for each group, it checks every span to see if its length
    is at least the minimum length. If it is, it calculates the new maximum length, which is the smallest of the max_len
    and the span's length. Then it adds the right semi-spans to the list of semi-spans.

    :param spans: A list of tuples, where each tuple represents a span with three integers (start, end, value).
    :type spans: List[Tuple[int, int, int]]
    :param min_len: The minimum length of the right semi-spans.
    :type min_len: int
    :param max_len: The maximum length of the right semi-spans.
    :type max_len: int

    :return: A list of tuples representing all right semi-spans that are within the specified length range.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        >>> _grouped_right_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], 1, 5)
        [(1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

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


def build_semi_spans(spans: List[Tuple[int, int, int]], min_len: int = 1, max_len: int = None) -> \
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

    .. code-block:: python

        >>> build_semi_spans([(0, 3, 0), (0, 5, 1), (3, 5, 0)], 1, 5)
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0), (1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

    """

    if min_len is None:
        min_len = 1

    if max_len is None:
        max_len = max(span[1] - span[0] for span in spans)

    return _grouped_left_semi_span_builder(spans, min_len, max_len) + _grouped_right_semi_span_builder(spans, min_len, max_len)


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

    if len(enzyme_sites) == max_index:  # non-enzymatic
        return build_non_enzymatic_spans((0, max_index, 0), min_len, max_len)

    if semi:
        spans = build_enzymatic_spans(max_index, enzyme_sites, missed_cleavages, min_len, None)
        semi_spans = build_semi_spans(spans, min_len, max_len)

        # filter spans by len
        spans = [span for span in spans if max_len >= span[1] - span[0] >= min_len] + semi_spans
    else:
        spans = build_enzymatic_spans(max_index, enzyme_sites, missed_cleavages, min_len, max_len)
    return spans


def calculate_span_coverage(spans: List[Tuple[int, int, int]], max_index: int, accumulate: bool = False) -> list[int]:
    """
    Calculates the coverage array for a given list of spans.

    :param spans: The list of spans.
    :type spans: List[Tuple[int, int, int]]
    :param max_index: The max index of the span.
    :type max_index: int
    :param accumulate: Whether to accumulate the coverage array.
    :type accumulate: bool
    :return: The coverage array.
    :rtype: list[int]

    .. code-block:: python

        >>> calculate_span_coverage([(0, 3, 0), (3, 6, 0), (6, 9, 0)], 9)
        [1, 1, 1, 1, 1, 1, 1, 1, 1]

        >>> calculate_span_coverage([(0, 3, 0), (0, 3, 0), (6, 9, 0)], 9)
        [1, 1, 1, 0, 0, 0, 1, 1, 1]

        >>> calculate_span_coverage([(0, 3, 0), (0, 3, 0), (6, 9, 0)], 9, accumulate=True)
        [2, 2, 2, 0, 0, 0, 1, 1, 1]

        >>> calculate_span_coverage([(0, 3, 0), (3, 6, 0), (6, 9, 0)], 12)
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]

        >>> calculate_span_coverage([(0, 3, 0), (3, 6, 0), (6, 9, 0)], 6)
        Traceback (most recent call last):
        ...
        IndexError: list assignment index out of range

    """

    cov_array = [0] * max_index
    for span in spans:
        for i in range(span[0], span[1]):
            if accumulate:
                cov_array[i] += 1
            else:
                cov_array[i] = 1

    return cov_array
