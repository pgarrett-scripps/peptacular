from typing import Union, List, Tuple, Dict
import regex as re

from peptacular.constants import VALID_ION_TYPES, FORWARD_ION_TYPES


def validate_ion_type(ion_type: str) -> None:
    """
    Validates if the given ion type is valid.

    This function checks if the ion type provided is one of the valid ion types. If not, it raises a ValueError.

    :param ion_type: The type of ion to validate.
    :type ion_type: str

    :raises ValueError: If the ion type is not valid.

    .. code-block:: python

        >>> validate_ion_type("a")  # No error raised

        >>> validate_ion_type("z")  # No error raised

        >>> validate_ion_type("h")  # Raises ValueError
        Traceback (most recent call last):
        ValueError: Ion type h is invalid.
    """

    if ion_type not in VALID_ION_TYPES:
        raise ValueError(f"Ion type {ion_type} is invalid.")


def is_forward(ion_type: str) -> bool:
    """
        Checks if the provided ion type is a forward ion.

        This function first validates the ion type and then checks if it's a forward ion.

        :param ion_type: The type of ion to check.
        :type ion_type: str

        :return: True if it's a forward ion, False otherwise.
        :rtype: bool

        .. code-block:: python

            >>> is_forward("a")  # Returns False
            False

            >>> is_forward("y")  # Returns True
            True

    """

    validate_ion_type(ion_type)
    return ion_type not in FORWARD_ION_TYPES


def _are_parentheses_balanced(text: str, open_char='(', closed_char=')') -> bool:
    """
    Check if parentheses in the given text are balanced or not.

    :param text: The input string to check for balanced parentheses.
    :type text: str

    :return: True if parentheses are balanced, False otherwise.
    :rtype: bool

    .. code-block:: python

        >>> _are_parentheses_balanced("PEPTIDE(1.2345)")
        True

        >>> _are_parentheses_balanced("PEPTIDE(1.2345")
        False

        >>> _are_parentheses_balanced("PEPTIDE((1.2345")
        False

        >>> _are_parentheses_balanced("PEPTIDE((1.2345)")
        False

    """
    stack = []
    for i, char in enumerate(text):
        if char == open_char:
            stack.append(i)
        elif char == closed_char:
            if not stack:
                return False
            stack.pop()
    return len(stack) == 0


def validate_parentheses(sequence: str) -> None:
    """
    Validate the parentheses in the given sequence.

    :param sequence: The sequence to validate.
    :type sequence: str

    :raises ValueError: If the parentheses are not balanced.

    .. code-block:: python
        >>> validate_parentheses("PEPTIDE(1.2345)")

        >>> validate_parentheses("PEPTIDE(1.2345")
        Traceback (most recent call last):
        ValueError: Incorrect modification notation in peptide sequence: "PEPTIDE(1.2345".

        >>> validate_parentheses("PEPTIDE(1.2345)[Amide")
        Traceback (most recent call last):
        ValueError: Incorrect modification notation in peptide sequence: "PEPTIDE(1.2345)[Amide".

    """
    # Check parentheses balance
    if not _are_parentheses_balanced(sequence, '(', ')'):
        raise ValueError(f'Incorrect modification notation in peptide sequence: "{sequence}".')

    if not _are_parentheses_balanced(sequence, '[', ']'):
        raise ValueError(f'Incorrect modification notation in peptide sequence: "{sequence}".')


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


def identify_regex_indexes(input_str: str, regex_str: str, offset: int = 0) -> List[int]:
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

        >>> identify_regex_indexes("PEPTIDE", "P")
        [0, 2]

        >>> identify_regex_indexes("PEPTIDE", "E")
        [1, 6]

        # More complex regex
        >>> identify_regex_indexes("PEPTIDEP", "P[ST]")
        [2]

        >>> identify_regex_indexes("PPPP", "PP")
        [0, 1, 2]

    """

    return [i[0] for i in identify_regex_indexes_range(input_str, regex_str, offset)]


def identify_regex_indexes_range(input_str: str, regex_str: str, offset: int = 0) -> List[Tuple[int, int]]:
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

        >>> identify_regex_indexes_range("PEPTIDE", "P")
        [(0, 1), (2, 3)]

        >>> identify_regex_indexes_range("PEPTIDE", "E")
        [(1, 2), (6, 7)]

        # More complex regex
        >>> identify_regex_indexes_range("PEPTIDEP", "P[ST]")
        [(2, 4)]

        >>> identify_regex_indexes_range("PPPP", "PP")
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


def map_brackets(text: str) -> List[Tuple[int, int]]:
    """
    Find matching pairs of brackets in a string.
    
    :param text: The input string.
    :type text: str
    
    :return: A list of tuples containing the start and end indices of matching pairs of brackets.
    :rtype: List[Tuple[int, int]]
    
    .. code-block:: python
    
        >>> map_brackets("a[b]c")
        [(1, 3)]

        >>> map_brackets("a[b][x]c")
        [(1, 3), (4, 6)]

        >>> map_brackets("a[b[c]d]e")
        [(1, 7)]

        >>> map_brackets("a[b[c]d]e]]")
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
                # Handle case of unmatched closing bracket if necessary
                pass

    # No need to handle unmatched opening brackets as they are ignored for outermost pairs

    return pairs


def remove_spans_from_text(text: str, pairs: List[Tuple[int, int]]) -> str:
    """
    Remove spans from a string.

    :param text: The input string.
    :type text: str
    :param pairs: A list of tuples containing the start and end indices of spans to remove.
    :type pairs: List[Tuple[int, int]]

    :return: The input string with the spans removed.
    :rtype: str

    .. code-block:: python

        >>> remove_spans_from_text("a[b]c", [(1, 3)])
        'ac'

        >>> remove_spans_from_text("a[b][x]c", [(1, 3), (4, 6)])
        'ac'

        >>> remove_spans_from_text("[x]a[b]c[x]", [(0, 2), (4, 6), (8, 10)])
        'ac'

        >>> remove_spans_from_text("a[b[c]d]e", [(1, 7)])
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

        >>> map_brackets_to_text("a[b]c")
        ('ac', {0: ['b']})

        >>> map_brackets_to_text("a[b][x]c")
        ('ac', {0: ['b', 'x']})

        >>> map_brackets_to_text("a[b]cd[1]e")
        ('acde', {0: ['b'], 2: ['1']})

        >>> map_brackets_to_text("a[b[c]d]e")
        ('ae', {0: ['b[c]d']})

        >>> map_brackets_to_text("[x]-a[b]c-[x]")
        ('-ac-', {-1: ['x'], 1: ['b'], 3: ['x']})

    """

    pairs = map_brackets(text)  # Reuse the map_brackets function to find pairs

    bracket_contents = {}
    offset = 0  # Keep track of the offset caused by previous brackets

    for start, end in pairs:
        # Adjust the start index by subtracting the offset
        adjusted_start = start - offset
        preceding_char_index = adjusted_start - 1

        if preceding_char_index not in bracket_contents:
            bracket_contents[preceding_char_index] = []
        # Map the adjusted start index to the bracket content
        bracket_contents[preceding_char_index].append(text[start+1:end])
        # Update the offset for the next iteration
        offset += end - start + 1


    sequence, spans = remove_spans_from_text(text, pairs), bracket_contents

    return sequence, spans


def map_brackets_to_text(text: str) -> Tuple[str, Dict[int, List[str]]]:

    sequence, spans = _map_brackets_to_text(text)
    sequence, spans = fix_multiple_mods(sequence, spans)
    return sequence, spans

def fix_multiple_mods(text: str, bracket_contents: Dict[int, List[str]]) -> Tuple[str, Dict[int, List[str]]]:
    """
    Fix multiple modifications in a string.

    :param text: The input string.
    :type text: str
    :param bracket_contents: A dictionary containing the start index of the bracket and the index of the preceeding character.
    :type bracket_contents: Dict[int, int]

    :return: A dictionary containing the start index of the bracket and the index of the preceeding character.
    :rtype: Dict[int, int]

    .. code-block:: python

        >>> fix_multiple_mods('a^2cde', {0: ['b'], 4: ['1']})
        ('acde', {0: ['b', 'b'], 2: ['1']})

        >>> fix_multiple_mods('a^2cde', {0: ['a', 'b'], 4: ['1']})
        ('acde', {0: ['a', 'b', 'b'], 2: ['1']})

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
            bracket_contents[index - 1].extend([last_mod]*(multiplier - 1))

            # Remove the '^n' from the text
            text = text[:index] + text[index + len(match.group(0)):]


        # Update the offset for the next iteration
        for key in list(bracket_contents.keys()):
            if key > index:
                bracket_contents[key-len_match] = bracket_contents[key]
                del bracket_contents[key]

    return text, bracket_contents