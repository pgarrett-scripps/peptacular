"""
The span module contains multiple functions for generating and working with spans. This module is mostly used within
the digest module, since spans essentially represent peptide sequences. Each span has a start, end and value component.
The start and end values reference the start and end index of the peptide within the protein, while the value component
is used to denote the number of missed cleavages the span contains.

Working with spans can be a more efficient way of processing peptide data, since a peptide can be reference with only
3 integers.
"""

from typing import Tuple, List
from itertools import groupby


def build_non_enzymatic_spans(span: Tuple[int, int, int], min_len: int = None, max_len: int = None) \
        -> List[Tuple[int, int, int]]:
    """
    Generates non-enymatic spans with span lengths <= max_len and >= min_len

    :param span: The input span for which to generate sub-spans.
    :type span: Tuple[int, int, int]
    :param min_len: The minimum length of spans to be generated, default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: The maximum length of spans to be generated, default is [None]. If None or any value greater than
                    span length - 1, max_len will be equal to span length - 1.
    :type max_len: int

    :return: All non-enymatic spans.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        # By default all spans are returned with lengths >= 1 and <= span length - 1
        >>> build_non_enzymatic_spans((0, 3, 0))
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

        # The span value for non-enymatic spans will always be 0
        >>> build_non_enzymatic_spans((0, 3, 2))
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

        # Can also explicitly specify min_len and max_len
        >>> build_non_enzymatic_spans((0, 3, 0), min_len=1, max_len=1)
        [(0, 1, 0), (1, 2, 0), (2, 3, 0)]

        # But it is not possible to generate spans >= span length - 1
        >>> build_non_enzymatic_spans((0, 3, 0), max_len=10)
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


def build_left_semi_spans(span: Tuple[int, int, int], min_len: int = None, max_len: int = None) \
        -> List[Tuple[int, int, int]]:
    """
    Generates left-semi spans with span lengths <= max_len and >= min_len. A left-semi span is any span
    which has the same start position as the parent span.

    :param span: The input span for which to generate sub-spans.
    :type span: Tuple[int, int, int]
    :param min_len: The minimum length of spans to be generated, default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: The maximum length of spans to be generated, default is [None]. If None or any value greater than
                    span length - 1, max_len will be equal to span length - 1.
    :type max_len: int

    :return: A list of all left-semi spans.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        # By default all spans are returned with lengths >= 1 and <= span length - 1
        >>> build_left_semi_spans((0, 3, 0))
        [(0, 2, 0), (0, 1, 0)]

        # Keeps the value of the original span
        >>> build_left_semi_spans((0, 3, 2))
        [(0, 2, 2), (0, 1, 2)]

        # Can also explicitly specify min_len and max_len
        >>> build_left_semi_spans((0, 3, 0), min_len=1, max_len=1)
        [(0, 1, 0)]

        # But it is not possible to generate spans >= span length - 1
        >>> build_left_semi_spans((0, 3, 0), max_len=10)
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
    Generates right-semi spans with span lengths <= max_len and >= min_len. A right-semi span is any span
    which has the same end position as the parent span.

    :param span: The input span for which to generate sub-spans.
    :type span: Tuple[int, int, int]
    :param min_len: The minimum length of spans to be generated, default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: The maximum length of spans to be generated, default is [None]. If None or any value greater than
                    span length - 1, max_len will be equal to span length - 1.
    :type max_len: int

    :return: A list of all right-semi spans.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        # By default all spans are returned with lengths >= 1 and <= span length - 1
        >>> build_right_semi_spans((0, 3, 0))
        [(1, 3, 0), (2, 3, 0)]

        # Keeps the value of the original span
        >>> build_right_semi_spans((0, 3, 2))
        [(1, 3, 2), (2, 3, 2)]

        # Can also explicitly specify min_len and max_len
        >>> build_right_semi_spans((0, 3, 0), min_len=1, max_len=1)
        [(2, 3, 0)]

        # But it is not possible to generate spans >= span length - 1
        >>> build_right_semi_spans((0, 3, 0), min_len=1, max_len=10)
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
                          min_len: int = None, max_len: int = None) -> List[Tuple[int, int, int]]:
    """
    Computes enzymatic spans for the given enzyme sites and missed cleavages.

    :param max_index: The max index of the span.
    :type max_index: int
    :param enzyme_sites: The list of indices representing enzyme sites.
    :type enzyme_sites: List[int]
    :param missed_cleavages: The number of allowed missed cleavages.
    :type missed_cleavages: int
    :param min_len: The minimum length of an enzymatic span, default is [None]. If None min_len will be equal to 1.
    :type min_len: int
    :param max_len: The maximum length of an enzymatic span, default is [None]. If None max_len will be equal to
                    max_index.
    :type max_len: int

    :return: A list of all enzymatic spans.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        >>> build_enzymatic_spans(5, [3], 1)
        [(0, 3, 0), (0, 5, 1), (3, 5, 0)]

        # Set min length
        >>> build_enzymatic_spans(5, [3], 1, min_len=5)
        [(0, 5, 1)]

        # Set max length
        >>> build_enzymatic_spans(5, [3], 1, max_len=3)
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


def _grouped_left_semi_span_builder(spans: List[Tuple[int, int, int]], min_len: int = None, max_len: int = None) -> \
        List[Tuple[int, int, int]]:
    """
    Efficiently generates all left-semi-spans from the given list of spans that have a length within the specified
    range. The input spans must be enzymatic spans where the values of the span represents the number of missed
    cleavages.

    :param spans: The input spans for which to generate sub-spans.
    :type spans: List[Tuple[int, int, int]]
    :param min_len: The minimum length of a span, default is [None]. If None min_len will be equal to 1.
    :type min_len: int
    :param max_len: The maximum length of a span, default is [None]. If None max_len will be equal to the largest span.
    :type max_len: int

    :return: A list of all left-semi spans.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        >>> _grouped_left_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=1, max_len=5)
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0)]

        >>> _grouped_left_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=None, max_len=None)
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0)]

    """

    if min_len is None:
        min_len = 1

    semi_spans = []
    spans = sorted(spans, key=lambda x: (x[0], -x[2]))
    for _, group in groupby(spans, key=lambda x: x[0]):
        group = list(group)
        for i, span in enumerate(group):
            span_len = span[1] - span[0]

            if span_len <= min_len:
                break

            new_max_len = span_len - 1
            if max_len is not None:
                new_max_len = min(max_len, new_max_len)

            if i == len(group) - 1:
                semi_spans.extend(build_left_semi_spans(span, min_len, new_max_len))
            else:
                next_span = group[i + 1]
                next_span_len = next_span[1] - next_span[0]
                new_min = max(min_len, next_span_len + 1)
                semi_spans.extend(build_left_semi_spans(span, new_min, new_max_len))

    return semi_spans


def _grouped_right_semi_span_builder(spans: List[Tuple[int, int, int]], min_len: int = None, max_len: int = None) -> \
        List[Tuple[int, int, int]]:
    """
    Efficiently generates all right-semi-spans from the given list of spans that have a length within the specified
    range. The input spans must be enzymatic spans where the values of the span represents the number of missed
    cleavages.

    :param spans: The input spans for which to generate sub-spans.
    :type spans: List[Tuple[int, int, int]]
    :param min_len: The minimum length of a span, default is [None]. If None min_len will be equal to 1.
    :type min_len: int
    :param max_len: The maximum length of a span, default is [None]. If None max_len will be equal to the largest span.
    :type max_len: int

    :return: A list of all right-semi spans.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        >>> _grouped_right_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=1, max_len=5)
        [(1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

        >>> _grouped_right_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=None, max_len=None)
        [(1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

    """

    if min_len is None:
        min_len = 1

    semi_spans = []
    spans = sorted(spans, key=lambda x: (x[1], -x[2]))
    for _, group in groupby(spans, key=lambda x: x[1]):
        group = list(group)
        for i, span in enumerate(group):
            span_len = span[1] - span[0]

            if span_len < min_len:
                break

            new_max_len = span_len - 1
            if max_len is not None:
                new_max_len = min(max_len, new_max_len)

            if i == len(group) - 1:
                semi_spans.extend(build_right_semi_spans(span, min_len, new_max_len))
            else:
                next_span = group[i + 1]
                next_span_len = next_span[1] - next_span[0]
                new_min = max(min_len, next_span_len + 1)
                semi_spans.extend(build_right_semi_spans(span, new_min, new_max_len))

    return semi_spans


def build_semi_spans(spans: List[Tuple[int, int, int]], min_len: int = None, max_len: int = None) -> \
        List[Tuple[int, int, int]]:
    """
    Efficiently generates all semi-spans from the given list of spans that have a length within the specified
    range. The input spans must be enzymatic spans where the values of the span represents the number of missed
    cleavages.


    :param spans: The list of spans.
    :type spans: List[Tuple[int, int, int]]
    :param min_len: The minimum length of spans to be generated, default is [None]. If None, min_len will be
                    equal to 1.
    :type min_len: int
    :param max_len: The maximum length of spans to be generated, default is [None]. If None or any value greater than
                    the max span length will be equal to max span length.

    :return: A list of all semi spans.
    :rtype: List[Tuple[int, int, int]]

    .. code-block:: python

        >>> build_semi_spans([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=1, max_len=5)
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0), (1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

        >>> build_semi_spans([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=None, max_len=None)
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0), (1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

    """

    return _grouped_left_semi_span_builder(spans, min_len, max_len) + \
        _grouped_right_semi_span_builder(spans, min_len, max_len)


def build_spans(max_index: int, enzyme_sites: List[int], missed_cleavages: int, min_len: int = None,
                max_len: int = None, semi: bool = False) -> List[Tuple[int, int, int]]:
    """
    Builds all spans for the given digestion parameters and enzyme sites

    :param max_index: The max index of the span.
    :type max_index: int
    :param enzyme_sites: The list of indices representing enzyme sites.
    :type enzyme_sites: List[int]
    :param missed_cleavages: The number of allowed missed cleavages.
    :type missed_cleavages: int
    :param min_len: The minimum length of a span, default is [None]. If None min_len will be equal to 1.
    :type min_len: int
    :param max_len: The maximum length of a span, default is [None]. If None max_len will be equal to max_index.
    :type max_len: int
    :param semi: Whether to compute and include semi spans, default is False.
    :type semi: bool

    :return: A list of all generated spans.
    :rtype: List[Tuple[int, int, int]]
    """

    if min_len is None:
        min_len = 1

    if max_len is None:
        max_len = max_index

    if len(enzyme_sites) == max_index-1:  # non-enzymatic
        return build_non_enzymatic_spans((0, max_index, 0), min_len, max_len)

    if semi:
        spans = build_enzymatic_spans(max_index, enzyme_sites, missed_cleavages, min_len, None)
        semi_spans = build_semi_spans(spans, min_len, max_len)

        # filter spans by len
        spans = [span for span in spans if max_len >= span[1] - span[0] >= min_len] + semi_spans
    else:
        spans = build_enzymatic_spans(max_index, enzyme_sites, missed_cleavages, min_len, max_len)

    return spans


def calculate_span_coverage(spans: List[Tuple[int, int, int]], max_index: int, accumulate: bool = False) -> List[int]:
    """
    Calculates the coverage array for a given list of spans.

    :param spans: The list of spans.
    :type spans: List[Tuple[int, int, int]]
    :param max_index: The max index of the span.
    :type max_index: int
    :param accumulate: Whether to accumulate the coverage array.
    :type accumulate: bool

    :return: The coverage array.
    :rtype: List[int]

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
