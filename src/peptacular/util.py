from peptacular.constants import VALID_ION_TYPES, FORWARD_IONS


def validate_ion_type(ion_type: str) -> None:
    """
    Validates if the given ion type is valid.

    This function checks if the ion type provided is one of the valid ion types. If not, it raises a ValueError.

    :param ion_type: The type of ion to validate.
    :type ion_type: str

    :raises ValueError: If the ion type is not valid.

    :Example:

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

        :Example:

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

    :Example:

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

    :return: None

    :raises ValueError: If the parentheses are not balanced.
    :Example:
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
