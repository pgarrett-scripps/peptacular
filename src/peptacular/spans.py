"""
The span module contains multiple functions for generating and working with spans. This module is mostly used within
the digest module, since spans essentially represent peptide sequences. Each span has a start, end and value component.
The start and end values reference the start and end index of the peptide within the protein, while the value component
is used to denote the number of missed cleavages the span contains.

Working with spans can be a more efficient way of processing the data, since a peptide can be reference with only
3 ints.
"""

from typing import List, Optional, Generator, Iterable
from itertools import groupby

from .types import Span


def build_non_enzymatic_spans(
    span: Span, min_len: Optional[int] = None, max_len: Optional[int] = None
) -> Generator[Span, None, None]:
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
        >>> list(build_non_enzymatic_spans((0, 3, 0)))
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

        # The span value for non-enymatic spans will always be 0
        >>> list(build_non_enzymatic_spans((0, 3, 2)))
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

        # Can also explicitly specify min_len and max_len
        >>> list(build_non_enzymatic_spans((0, 3, 0), min_len=1, max_len=1))
        [(0, 1, 0), (1, 2, 0), (2, 3, 0)]

        # But it is not possible to generate spans >= span length - 1
        >>> list(build_non_enzymatic_spans((0, 3, 0), max_len=10))
        [(0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 3, 0), (2, 3, 0)]

    """

    if min_len is None:
        min_len = 1

    max_span = span[1] - span[0] - 1
    if max_len is None:
        max_len = max_span
    max_len = min(max_len, max_span)

    start, end, _ = span
    return (
        (i, j, 0)
        for i in range(start, end)
        for j in range(i + min_len, min(end + 1, i + max_len + 1))
    )


def build_left_semi_spans(
    span: Span, min_len: Optional[int] = None, max_len: Optional[int] = None
) -> Generator[Span, None, None]:
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
        >>> list(build_left_semi_spans((0, 3, 0)))
        [(0, 2, 0), (0, 1, 0)]

        # Keeps the value of the original span
        >>> list(build_left_semi_spans((0, 3, 2)))
        [(0, 2, 2), (0, 1, 2)]

        # Can also explicitly specify min_len and max_len
        >>> list(build_left_semi_spans((0, 3, 0), min_len=1, max_len=1))
        [(0, 1, 0)]

        # But it is not possible to generate spans >= span length - 1
        >>> list(build_left_semi_spans((0, 3, 0), max_len=10))
        [(0, 2, 0), (0, 1, 0)]

    """

    if min_len is None:
        min_len = 1

    if max_len is None:
        max_len = span[1] - span[0]

    start, end, value = span
    new_end = min(start + max_len, end - 1)
    return (
        (start, i, value) for i in range(new_end, start - 1, -1) if i - start >= min_len
    )


def build_right_semi_spans(
    span: Span, min_len: int = 1, max_len: Optional[int] = None
) -> Generator[Span, None, None]:
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
        >>> list(build_right_semi_spans((0, 3, 0)))
        [(1, 3, 0), (2, 3, 0)]

        # Keeps the value of the original span
        >>> list(build_right_semi_spans((0, 3, 2)))
        [(1, 3, 2), (2, 3, 2)]

        # Can also explicitly specify min_len and max_len
        >>> list(build_right_semi_spans((0, 3, 0), min_len=1, max_len=1))
        [(2, 3, 0)]

        # But it is not possible to generate spans >= span length - 1
        >>> list(build_right_semi_spans((0, 3, 0), min_len=1, max_len=10))
        [(1, 3, 0), (2, 3, 0)]

    """

    if min_len is None:
        min_len = 1

    if max_len is None:
        max_len = span[1] - span[0]

    start, end, value = span
    new_start = max(start + 1, end - max_len)
    return ((i, end, value) for i in range(new_start, end + 1) if end - i >= min_len)


def build_enzymatic_spans(
    max_index: int,
    enzyme_sites: List[int],
    missed_cleavages: int,
    min_len: Optional[int] = None,
    max_len: Optional[int] = None,
) -> Generator[Span, None, None]:
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

        >>> list(build_enzymatic_spans(5, [3], 1))
        [(0, 3, 0), (0, 5, 1), (3, 5, 0)]

        # Set min length
        >>> list(build_enzymatic_spans(5, [3], 1, min_len=5))
        [(0, 5, 1)]

        # Set max length
        >>> list(build_enzymatic_spans(5, [3], 1, max_len=3))
        [(0, 3, 0), (3, 5, 0)]

    """

    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = max_index

    enzyme_sites = set(enzyme_sites)
    enzyme_sites.add(0)
    enzyme_sites.add(max_index)
    enzyme_sites = sorted(enzyme_sites)

    for i, start_site in enumerate(enzyme_sites):
        for j, end_site in enumerate(enzyme_sites[i + 1 : i + missed_cleavages + 2]):
            if min_len <= (end_site - start_site) <= max_len:
                yield start_site, end_site, j


def _grouped_left_semi_span_builder(
    spans: List[Span], min_len: Optional[int] = None, max_len: Optional[int] = None
) -> Generator[Span, None, None]:
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

        >>> list(_grouped_left_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=1, max_len=5))
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0)]

        >>> list(_grouped_left_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=None, max_len=None))
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0)]

    """

    if min_len is None:
        min_len = 1

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
                yield from build_left_semi_spans(span, min_len, new_max_len)
            else:
                next_span = group[i + 1]
                next_span_len = next_span[1] - next_span[0]
                new_min = max(min_len, next_span_len + 1)
                yield from build_left_semi_spans(span, new_min, new_max_len)


def _grouped_right_semi_span_builder(
    spans: List[Span], min_len: Optional[int] = None, max_len: Optional[int] = None
) -> Generator[Span, None, None]:
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

        >>> list(_grouped_right_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=1, max_len=5))
        [(1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

        >>> list(_grouped_right_semi_span_builder([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=None, max_len=None))
        [(1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

    """

    if min_len is None:
        min_len = 1

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
                yield from build_right_semi_spans(span, min_len, new_max_len)
            else:
                next_span = group[i + 1]
                next_span_len = next_span[1] - next_span[0]
                new_min = max(min_len, next_span_len + 1)
                yield from build_right_semi_spans(span, new_min, new_max_len)


def build_semi_spans(
    spans: List[Span], min_len: Optional[int] = None, max_len: Optional[int] = None
) -> Generator[Span, None, None]:
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

        >>> list(build_semi_spans([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=1, max_len=5))
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0), (1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

        >>> list(build_semi_spans([(0, 3, 0), (0, 5, 1), (3, 5, 0)], min_len=None, max_len=None))
        [(0, 4, 1), (0, 2, 0), (0, 1, 0), (3, 4, 0), (1, 3, 0), (2, 3, 0), (1, 5, 1), (2, 5, 1), (4, 5, 0)]

    """

    yield from _grouped_left_semi_span_builder(spans, min_len, max_len)
    yield from _grouped_right_semi_span_builder(spans, min_len, max_len)


def build_spans(
    max_index: int,
    enzyme_sites: Iterable[int],
    missed_cleavages: int,
    min_len: Optional[int] = None,
    max_len: Optional[int] = None,
    semi: bool = False,
) -> Generator[Span, None, None]:
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

    enzyme_sites = sorted(set(enzyme_sites))

    if len(enzyme_sites) == max_index + 1:  # non-enzymatic case
        yield from build_non_enzymatic_spans((0, max_index, 0), min_len, max_len)
        return  # Exit early since we only need non-enzymatic spans

    spans = list(
        build_enzymatic_spans(
            max_index,
            enzyme_sites,
            missed_cleavages,
            min_len,
            None if semi else max_len,
        )
    )

    if semi:
        semi_spans = build_semi_spans(spans, min_len, max_len)
        for span in spans:
            if max_len >= span[1] - span[0] >= min_len:
                yield span
        yield from semi_spans
    else:
        yield from spans


def calculate_span_coverage(
    spans: List[Span], max_index: int, accumulate: bool = False
) -> List[int]:
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
