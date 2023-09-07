from typing import Union, List, Tuple
import re

from peptacular.constants import VALID_ION_TYPES, FORWARD_IONS


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
    return ion_type not in FORWARD_IONS


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

    """

    regex_indexes = [match.start() + offset for match in re.finditer(regex_str, input_str)]
    return regex_indexes


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
