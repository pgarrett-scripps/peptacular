from typing import Union, List, Tuple, Dict
import regex as re


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

    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            return val


def get_regex_match_indices(input_str: str, regex_str: str, offset: int = 0) -> List[int]:
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

        >>> get_regex_match_indices("PEPTIDE", "P")
        [0, 2]

        >>> get_regex_match_indices("PEPTIDE", "E")
        [1, 6]

        # More complex regex
        >>> get_regex_match_indices("PEPTIDE", "P[ST]")
        [2]

        >>> get_regex_match_indices("PPPP", "PP")
        [0, 1, 2]

    """

    return [i[0] for i in get_regex_match_range(input_str, regex_str, offset)]


def get_regex_match_range(input_str: str, regex_str: str, offset: int = 0) -> List[Tuple[int, int]]:
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

        >>> get_regex_match_range("PEPTIDE", "E")
        [(1, 2), (6, 7)]

        # More complex regex
        >>> get_regex_match_range("PEPTIDE", "P[ST]")
        [(2, 4)]

        >>> get_regex_match_range("PPPP", "PP")
        [(0, 2), (1, 3), (2, 4)]

    """

    return [(match.start() + offset, match.end() + offset) for match in re.finditer(regex_str, input_str,
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


def _map_brackets_to_preceding_char(text: str) -> List[Tuple[int, int]]:
    """
    Find matching pairs of brackets in a string.
    
    :param text: The input string.
    :type text: str
    
    :return: A list of tuples containing the start and end indices of matching pairs of brackets.
    :rtype: List[Tuple[int, int]]
    
    .. code-block:: python
    
        >>> _map_brackets_to_preceding_char("a[b]c")
        [(1, 3)]

        >>> _map_brackets_to_preceding_char("a[b][x]c")
        [(1, 3), (4, 6)]

        >>> _map_brackets_to_preceding_char("a[b[c]d]e")
        [(1, 7)]

        >>> _map_brackets_to_preceding_char("a[b[c]d]e]]")
        Traceback (most recent call last):
        ValueError: Unmatched brackets in input text.

        >>> _map_brackets_to_preceding_char("a][bc")
        Traceback (most recent call last):
        ValueError: Unmatched brackets in input text.

        >>> _map_brackets_to_preceding_char("a]bc")
        Traceback (most recent call last):
        ValueError: Unmatched brackets in input text.

        >>> _map_brackets_to_preceding_char("a[bc")
        Traceback (most recent call last):
        ValueError: Unmatched brackets in input text.

    """

    if text.count('[') != text.count(']'):
        raise ValueError('Unmatched brackets in input text.')

    stack = []  # Stack to keep track of opening bracket indices for nesting
    pairs = []  # List to store pairs of start and end indices for outermost brackets

    for i, char in enumerate(text):
        if char == '[':  # Opening bracket
            if not stack:  # If stack is empty, this is an outermost bracket
                stack.append((i, 'start'))
            else:
                stack.append((i, 'nested'))
        elif char == ']':  # Closing bracket
            if stack:  # Ensure stack is not empty
                start_index, bracket_type = stack.pop()  # Get last opening bracket
                if not stack:  # If stack is empty now, we've closed an outermost bracket
                    pairs.append((start_index, i))
                # If stack is not empty, we're still inside a nested structure
            else:
                raise ValueError('Unmatched brackets in input text.')

    # No need to handle unmatched opening brackets as they are ignored for outermost pairs

    return pairs


def _strip_bracket_spans(text: str, pairs: List[Tuple[int, int]]) -> str:
    """
    Remove spans from a string.

    :param text: The input string.
    :type text: str
    :param pairs: A list of tuples containing the start and end indices of spans to remove.
    :type pairs: List[Tuple[int, int]]

    :return: The input string with the spans removed.
    :rtype: str

    .. code-block:: python

        >>> _strip_bracket_spans("a[b]c", [(1, 3)])
        'ac'

        >>> _strip_bracket_spans("a[b][x]c", [(1, 3), (4, 6)])
        'ac'

        >>> _strip_bracket_spans("[x]a[b]c[x]", [(0, 2), (4, 6), (8, 10)])
        'ac'

        >>> _strip_bracket_spans("a[b[c]d]e", [(1, 7)])
        'ae'

    """

    result = []
    last_end = 0
    for start, end in pairs:
        # Add text up to the start of the current span
        result.append(text[last_end:start])
        # Update last_end to the end of the current span, to skip its content
        last_end = end + 1
    # Add any remaining text after the last span
    result.append(text[last_end:])
    return ''.join(result)


def _map_brackets_to_text(text: str) -> Tuple[str, Dict[int, List[str]]]:
    """
    Map the content of brackets to the preceeding character in a string.

    :param text: The input string.
    :type text: str

    :return: A dictionary containing the start index of the bracket and the index of the preceeding character.
    :rtype: Dict[int, int]

    .. code-block:: python

        >>> _map_brackets_to_text("a[b]c")
        ('ac', {0: ['b']})

        >>> _map_brackets_to_text("a[b][x]c")
        ('ac', {0: ['b', 'x']})

        >>> _map_brackets_to_text("a[b]cd[1]e")
        ('acde', {0: ['b'], 2: ['1']})

        >>> _map_brackets_to_text("a[b[c]d]e")
        ('ae', {0: ['b[c]d']})

        >>> _map_brackets_to_text("[x]-a[b]c-[x]")
        ('-ac-', {-1: ['x'], 1: ['b'], 3: ['x']})

        >>> _map_brackets_to_text("[1]?[2]-a[b]c")
        ('?-ac', {-1: ['1'], 0: ['2'], 2: ['b']})

        >>> _map_brackets_to_text("[b]ac")
        ('ac', {-1: ['b']})

    """

    pairs = _map_brackets_to_preceding_char(text)  # Reuse the map_brackets function to find pairs

    bracket_contents = {}
    offset = 0  # Keep track of the offset caused by previous brackets

    for start, end in pairs:
        # Adjust the start index by subtracting the offset
        adjusted_start = start - offset
        preceding_char_index = adjusted_start - 1

        if preceding_char_index not in bracket_contents:
            bracket_contents[preceding_char_index] = []

        # Map the adjusted start index to the bracket content
        bracket_contents[preceding_char_index].append(text[start + 1:end])
        # Update the offset for the next iteration
        offset += end - start + 1

    sequence, spans = _strip_bracket_spans(text, pairs), bracket_contents

    return sequence, spans


def map_bracket_content_to_index(text: str, ignore_index_error: bool = False) -> Tuple[str, Dict[int, List[str]]]:
    """
    Map the content of brackets to the preceeding character in a string.

    :param text: The input string.
    :type text: str
    :return: A dictionary containing the start index of the bracket and the index of the preceeding character.
    :rtype: Dict[int, int]

    .. code-block:: python

        >>> map_bracket_content_to_index("a[b]c")
        ('ac', {0: ['b']})

        >>> map_bracket_content_to_index("a[b][x]c")
        ('ac', {0: ['b', 'x']})

        >>> map_bracket_content_to_index("a[b]cd[1]e")
        ('acde', {0: ['b'], 2: ['1']})

        >>> map_bracket_content_to_index("[1]?[2]-a[b]c")
        ('?-ac', {0: ['1'], 1: ['2'], 2: ['b']})

        >>> map_bracket_content_to_index("[b]ac")
        Traceback (most recent call last):
        ValueError: Invalid sequence

    """
    sequence, spans = _map_brackets_to_text(text)
    sequence, spans = _fix_sequential_brackets(sequence, spans)
    sequence, spans = _fix_preceding_brackets(sequence, spans, ignore_index_error)

    return sequence, spans


def _fix_sequential_brackets(text: str, bracket_contents: Dict[int, List[str]]) -> Tuple[str, Dict[int, List[str]]]:
    """
    Fix multiple modifications in a string.

    :param text: The input string.
    :type text: str
    :param bracket_contents: A dictionary containing the start index of the bracket and the index of the
    preceeding character.
    :type bracket_contents: Dict[int, int]

    :return: A dictionary containing the start index of the bracket and the index of the preceeding character.
    :rtype: Dict[int, int]

    .. code-block:: python

        >>> _fix_sequential_brackets('a^2cde', {0: ['b'], 4: ['1']})
        ('acde', {0: ['b', 'b'], 2: ['1']})

        >>> _fix_sequential_brackets('a^2cde', {0: ['a', 'b'], 4: ['1']})
        ('acde', {0: ['a', 'b', 'b'], 2: ['1']})

        >>> _fix_sequential_brackets('^2?-ac', {-1: ['1'], 2: ['2'], 4: ['b']})
        ('?-ac', {-1: ['1', '1'], 0: ['2'], 2: ['b']})

    """

    pattern = r'\^(\d+)'
    matches = re.finditer(pattern, text)

    # Process each match in reverse order (to not mess up the indices)
    for match in reversed(list(matches)):
        index = match.start()
        multiplier = int(match.group(1))  # The digit after '^' indicates how many times to duplicate the mod

        len_match = 1 + len(str(multiplier))

        # Check if there is a modification at the index before '^'
        if index - 1 in bracket_contents:
            last_mod = bracket_contents[index - 1][-1]
            # Duplicate the modification
            bracket_contents[index - 1].extend([last_mod] * (multiplier - 1))

            # Remove the '^n' from the text
            text = text[:index] + text[index + len(match.group(0)):]

        # Update the offset for the next iteration
        for key in list(bracket_contents.keys()):
            if key > index:
                bracket_contents[key - len_match] = bracket_contents[key]
                del bracket_contents[key]

    return text, bracket_contents


def _fix_preceding_brackets(sequence: str, mods: Dict[int, List[str]], ignore_index_error) -> Tuple[str, Dict[int, List[str]]]:
    """
    Fix the indices of the modifications in a string.

    :param sequence: The input string.
    :type sequence: str
    :param mods: A dictionary containing the start index of the bracket and the index of the preceeding character.
    :type mods: Dict[int, int]

    :return: A dictionary containing the start index of the bracket and the index of the preceeding character.
    :rtype: Dict[int, int]

    .. code-block:: python

        >>> _fix_preceding_brackets('?-ac', {-1: ['1', '1'], 0: ['2'], 2: ['b']}, False)
        ('?-ac', {0: ['1', '1'], 1: ['2'], 2: ['b']})

        >>> _fix_preceding_brackets('ac', {-1: ['b']}, False)
        Traceback (most recent call last):
        ValueError: Invalid sequence

        >>> _fix_preceding_brackets('ac', {-1: ['b']}, True)
        ('ac', {0: ['b']})

    """

    shift_index = 0
    for i, c in enumerate(sequence):
        if c not in '?-':
            shift_index = i
            break

    # update all dict keys below shift index
    new_mods = {}
    for k in list(mods.keys()):
        if k < shift_index:

            if shift_index == k + 1:
                if not ignore_index_error:
                    raise ValueError('Invalid sequence: A modification is on the wrong side!')

            new_mods[k + 1] = mods[k]
        else:
            new_mods[k] = mods[k]

    return sequence, new_mods
