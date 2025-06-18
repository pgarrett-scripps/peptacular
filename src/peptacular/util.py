"""
Utils.py
"""

import warnings
from typing import Union, List, Tuple, Dict, Generator, Optional
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


def get_regex_match_indices(
    input_str: str, regex_str: Union[str, regex.Pattern], offset: int = 0
) -> Generator[int, None, None]:
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
            warnings.warn(
                message="The regex pattern has a match with a none zero length. Using start index + 1 for the match.",
            )
            yield match.start() + offset + 1
        else:
            yield match.start() + offset


def get_regex_match_range(
    input_str: str, regex_str: Union[str, regex.Pattern], offset: int = 0
) -> List[Tuple[int, int]]:
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
        return [
            (i.start() + offset, i.end() + offset)
            for i in regex_str.finditer(input_str)
        ]

    return [
        (match.start() + offset, match.end() + offset)
        for match in regex.finditer(regex_str, input_str, overlapped=True)
    ]


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
        raise ValueError(f"Start of span should be non-negative, got {start}.")
    if end < 0:
        raise ValueError(f"End of span should be non-negative, got {end}.")
    if start > end:
        raise ValueError(
            f"Start of span: {start}, should be less than or equal to end of span: {end}."
        )


def _construct_ambiguity_intervals(
    counts: List[int], reverse: bool
) -> List[Tuple[int, int]]:
    """
    Construct intervals for sequences of zeros in the counts list. When reverse is false, start from the left hand side
    and move to the right. When reverse is true, start from the right hand side and move to the left. Intervals start
    at 0 and end on any positive value. Both are inclusive. Returned intervals should be in forwards format, that is
    have a starting value less than the ending value.

    :param counts: List of integers (typically counts)
    :param reverse: If True, reverse the list before processing
    :return: List of intervals [start, end] indicating runs of zeros

    .. code-block:: python

        # [0, 1, 1, 1, 0, 0, 0]
        # [1, 1, 0, 0, 1, 1, 1] # ambiguity
        >>> _construct_ambiguity_intervals([0, 1, 1, 1, 0, 0, 0], reverse=False)
        [(0, 1), (4, 6)]

        # [0, 0, 1, 1, 1, 1, 0]
        # [1, 1, 0, 0, 0, 1, 1] # ambiguity
        >>> _construct_ambiguity_intervals([0, 0, 1, 1, 1, 0, 0], reverse=True)
        [(0, 1), (4, 6)]

        >>> _construct_ambiguity_intervals([0, 1, 1, 1, 0, 0, 1], reverse=False)
        [(0, 1), (4, 6)]
    """

    if reverse:
        ambiguity_intervals = _construct_ambiguity_intervals(
            counts[::-1], reverse=False
        )
        ambiguity_intervals = [
            (len(counts) - 1 - end, len(counts) - 1 - start)
            for start, end in ambiguity_intervals
        ]
        # sort the intervals
        ambiguity_intervals.sort(key=lambda x: x[0])
        return ambiguity_intervals

    ambiguity_intervals = []
    current_interval = None
    for i, cnt in enumerate(counts):
        if cnt == 0:
            if current_interval is not None:
                current_interval = (current_interval[0], i)
            else:
                current_interval = (i, i)

        else:
            if current_interval is not None:
                current_interval = (current_interval[0], i)
                ambiguity_intervals.append(current_interval)
                current_interval = None
            else:
                continue

    if current_interval is not None:
        current_interval = (current_interval[0], len(counts) - 1)
        ambiguity_intervals.append(current_interval)

    return ambiguity_intervals


def _combine_ambiguity_intervals(
    *interval_lists: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    """
    Merge multiple lists of ambiguity intervals into a single list of common ambiguity intervals.

    This function identifies positions that are ambiguous across all provided interval lists.
    For a position to be considered ambiguous in the result, it must be contained in at least
    one interval from each input list. The function then constructs optimized intervals
    covering these common ambiguous positions.

    Intervals are represented as tuples (start, end) where:
    - start is inclusive
    - end is exclusive

    Intervals with identical start and end values (zero-length intervals) are removed.

    :param interval_lists: Variable number of lists containing ambiguity intervals
    :type interval_lists: List[Tuple[int, int]]

    :return: A list of merged intervals representing positions that are ambiguous across all input lists
    :rtype: List[Tuple[int, int]]

    .. code-block:: python

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0,1)])
        [(0, 1)]

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0,1), (4,5)])
        [(0, 1), (4, 5)]

        >>> _combine_ambiguity_intervals([(0, 1), (4, 6)], [(0, 4), (5, 6)])
        [(0, 1), (5, 6)]

        >>> _combine_ambiguity_intervals([(2, 5)], [(3, 6)])
        [(3, 5)]

        >>> _combine_ambiguity_intervals([(0, 1)], [(4, 6)])
        []
    """

    # First, collect all unique intervals from input lists
    all_intervals = set()
    for interval_list in interval_lists:
        for interval in interval_list:
            all_intervals.add(interval)

    # Remove intervals where start == end (these are not ambiguous)
    filtered_intervals = {(start, end) for start, end in all_intervals if start != end}

    # Find all possible indices that are covered by any interval
    all_indices = set()
    for start, end in filtered_intervals:
        for i in range(start, end):
            all_indices.add(i)

    # For each index, check if it's contained in at least one interval from each input list
    common_indices = set()
    for idx in all_indices:
        is_common = True
        for interval_list in interval_lists:
            if not any(start <= idx < end for start, end in interval_list):
                is_common = False
                break
        if is_common:
            common_indices.add(idx)

    # If no common indices found, return empty list
    if not common_indices:
        return []

    # Construct new intervals from the common indices
    result = []
    if common_indices:
        sorted_indices = sorted(common_indices)
        start = sorted_indices[0]
        for i in range(1, len(sorted_indices)):
            if sorted_indices[i] > sorted_indices[i - 1] + 1:
                # Gap found, close the current interval and start a new one
                result.append((start, sorted_indices[i - 1] + 1))
                start = sorted_indices[i]
        # Add the last interval
        result.append((start, sorted_indices[-1] + 1))

    return result


def _get_mass_shift_interval(
    forward_coverage: List[int], reverse_coverage: List[int]
) -> Optional[Tuple[int, int]]:
    """
        Determine the interval where a mass shift should be placed based on fragment ion coverage.

    This function examines the forward and reverse ion coverage to identify the region where
    a mass shift (such as a modification) should be positioned. It returns the start and end
    indices (inclusive) of this region, or None if no suitable region is found.

    The mass shift interval is determined by:
    1. Finding the highest position with forward ion coverage
    2. Finding the lowest position with reverse ion coverage
    3. The mass shift belongs between these two positions

    :param forward_coverage: Binary list indicating forward ion coverage (1) or no coverage (0)
    :type forward_coverage: List[int]
    :param reverse_coverage: Binary list indicating reverse ion coverage (1) or no coverage (0)
    :type reverse_coverage: List[int]

    :return: A tuple containing the start and end indices (inclusive) for the mass shift,
             or None if there is no valid interval
    :rtype: Optional[Tuple[int, int]]

    .. code-block:: python

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,0,1,1,1])
        (3, 3)

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,1,1,1,1])
        (3, 3)

        >>> _get_mass_shift_interval([1,1,0,0,0,0,0], [0,0,0,0,1,1,1])
        (2, 3)

        >>> _get_mass_shift_interval([0,0,0,0,0,0,0], [0,0,0,0,1,1,1])
        (0, 3)

        >>> _get_mass_shift_interval([1,1,1,0,0,0,0], [0,0,0,0,0,0,0])
        (3, 6)

        >>> _get_mass_shift_interval([1,1,1,1,1,0,0], [0,0,0,0,1,1,1]) # None

    """

    highest_forward_fragment = [i for i, cnt in enumerate(forward_coverage) if cnt > 0]
    if len(highest_forward_fragment) == 0:
        highest_forward_fragment = -1
    else:
        highest_forward_fragment = max(highest_forward_fragment)

    highest_reverse_fragment = [i for i, cnt in enumerate(reverse_coverage) if cnt > 0]
    if len(highest_reverse_fragment) == 0:
        highest_reverse_fragment = len(reverse_coverage)
    else:
        highest_reverse_fragment = min(highest_reverse_fragment)

    if highest_forward_fragment >= highest_reverse_fragment:
        return None

    if highest_forward_fragment == highest_reverse_fragment - 1:
        return (highest_forward_fragment + 1, highest_forward_fragment + 1)

    # if there is a gap between the two, return the gap
    return (highest_forward_fragment + 1, highest_reverse_fragment - 1)
