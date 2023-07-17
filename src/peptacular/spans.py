from functools import wraps
from typing import Tuple, List, Any, Optional
from itertools import groupby


def span_processing(func):
    """
    A decorator to enforce constraints on the span and min_len and max_len arguments of decorated functions.
    The span should be a tuple of length 3 and min_len and max_len should be None or an integer.
    """

    @wraps(func)
    def wrapper(span: Tuple[int, int, Any], min_len: Optional[int] = None, max_len: Optional[int] = None):
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


def _validate_span(span: Tuple[int, int, Any]) -> None:
    """
    This function checks if a given span is valid.
    Raises:
        ValueError: If the span is not valid.

    """
    start, end, value = span
    if start < 0:
        raise ValueError(f'Start of span should be non-negative, got {start}.')
    if end < 0:
        raise ValueError(f'End of span should be non-negative, got {end}.')
    if start > end:
        raise ValueError(f'Start of span {start} should be less than or equal to end of span {end}.')


@span_processing
def build_non_enzymatic_spans(span: Tuple[int, int, Any], min_len: int, max_len: int) -> List[Tuple[int, int, Any]]:
    """
    This function generates and returns all possible sub-spans of the given span. These sub-spans have lengths ranging
    from `min_len` to `max_len`. The sub-spans are "non-enzymatic", meaning they are direct subsets of the given span
    and not created through any enzymatic action or process.

    Args:
        span (Tuple[int, int, int]): A tuple representing the original span, structured as (start, end, value).
        min_len (int, optional): The minimum length of sub-spans to be generated. Defaults to 1.
        max_len (int, optional): The maximum length of sub-spans to be generated. If not provided, it defaults to the
                                 length of the input span.

    Returns:
        List[Tuple[int, int, int]]: A list of all possible sub-spans as tuples, each structured as (start, end, value).
    """
    start, end, value = span
    return [(i, j, 0) for i in range(start, end + 1) for j in range(i + min_len, min(end + 1, i + max_len + 1))]


@span_processing
def build_left_semi_spans(span: Tuple[int, int, Any], min_len: int, max_len: int) -> List[Tuple[int, int, Any]]:
    """
    This function generates and returns all possible sub-spans of the given span starting from the left.
    These sub-spans have lengths ranging from `min_len` to `max_len`.

    Args:
        span (Tuple[int, int, int]): A tuple representing the original span, structured as (start, end, value).
        min_len (int, optional): The minimum length of sub-spans to be generated. Defaults to 1.
        max_len (int, optional): The maximum length of sub-spans to be generated. If not provided, it defaults to the
                                 length of the input span.

    Returns:
        List[Tuple[int, int, int]]: A list of all possible left sub-spans as tuples, each structured as
                                    (start, end, value).
    """
    if min_len > max_len:
        return []

    start, end, value = span
    new_end = min(start + max_len, end - 1)
    return [(start, i, value) for i in range(new_end, start - 1, -1) if i - start >= min_len]


@span_processing
def build_right_semi_spans(span: Tuple[int, int, Any], min_len: int, max_len: int) -> List[Tuple[int, int, Any]]:
    """
    This function generates and returns all possible sub-spans of the given span starting from the right.
    These sub-spans have lengths ranging from `min_len` to `max_len`.

    Args:
        span (Tuple[int, int, int]): A tuple representing the original span, structured as (start, end, value).
        min_len (int, optional): The minimum length of sub-spans to be generated. Defaults to 1.
        max_len (int, optional): The maximum length of sub-spans to be generated. If not provided, it defaults to the
                                 length of the input span.

    Returns:
        List[Tuple[int, int, int]]: A list of all possible right sub-spans as tuples, each structured as
                                    (start, end, value).
    """

    if min_len > max_len:
        return []

    start, end, value = span
    new_start = max(start + 1, end - max_len)
    return [(i, end, value) for i in range(new_start, end + 1) if end - i >= min_len]


def span_to_sequence(sequence: str, span: Tuple[int, int, Any]) -> str:
    """
    This function takes a sequence and a span as input, then returns the subsequence of the input sequence
    that corresponds to the provided span.

    Args:
        sequence (str): The original sequence from which a subsequence will be extracted.
        span (Tuple[int, int, int]): A tuple representing the span of the subsequence to be extracted, structured as
                                     (start, end, value).

    Returns:
        str: The subsequence of the input sequence as defined by the start and end of the provided span.
    """

    _validate_span(span)
    return sequence[span[0]:span[1]]


def get_enzymatic_spans(max_index: int, enzyme_sites: List[int], missed_cleavages: int,
                        min_len: int = None, max_len: int = None) -> List[Tuple[int, int, int]]:
    """
    Computes the enzymatic spans for the given enzymatic sites and the number of missed
    cleavages.

    Parameters:
        max_index (int): The max index of the span.
        enzyme_sites (List[int]): The list of indices in the protein sequence that represent enzyme sites.
        missed_cleavages (int): The number of allowed missed cleavages.
        min_len (int): The minimum length of an enzymatic span.
        max_len (int): The maximum length of an enzymatic span.

    Returns:
        List[Tuple[int, int, int]]: A list of tuples, where each tuple represents a span and contains three integers:
                                    the start index, the end index, and the number of missed cleavages.
    """

    if min_len is None:
        min_len = 1
    if max_len is None:
        max_len = max_index

    enzyme_sites = [0] + enzyme_sites + [max_index]

    spans = []
    for i, start_site in enumerate(enzyme_sites):
        spans.extend([(start_site, end_site, j) for j, end_site in enumerate(enzyme_sites[i+1:i+missed_cleavages+2])])

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
    is at least the minimum length. If it is, it calculates the new maximum length, which is the smaller of the max_len
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
    is at least the minimum length. If it is, it calculates the new maximum length, which is the smaller of the max_len
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
    Computes the semi spans for a given list of spans based on the given minimum and maximum length criteria.

    Parameters:
        spans (List[Tuple[int, int, int]]): The list of spans for which the semi spans will be computed.
        min_len (int): The minimum length of a semi span.
        max_len (int): The maximum length of a semi span.

    Returns:
        List[Tuple[int, int, int]]: A list of tuples, where each tuple represents a semi span and contains three
                                    integers: the start index, the end index, and the number of missed cleavages.
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
    Builds spans for a given sequence based on the given enzyme sites, the number of missed cleavages, and the minimum
    and maximum length criteria. If semi is set to True, semi spans are also computed and added to the list of spans.

    Parameters:
        max_index (int): The max index of the span.
        enzyme_sites (List[int]): The list of indices in the sequence that represent enzyme sites.
        missed_cleavages (int): The number of allowed missed cleavages.
        min_len (int): The minimum length of a span.
        max_len (int): The maximum length of a span.
        semi (bool): Whether to compute and include semi spans.

    Returns:
        List[Tuple[int, int, int]]: A list of tuples, where each tuple represents a span and contains three integers:
                                    the start index, the end index, and the number of missed cleavages.
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
