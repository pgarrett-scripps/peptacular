from typing import Union, List, Tuple
import regex as re

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

        >>> convert_type("1.0")
        1.0

        >>> convert_type("abc")
        'abc'
    """
    for conversion in (int, float):
        try:
            return conversion(val)
        except ValueError:
            continue
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

    regex_indexes = [match.start() + offset for match in re.finditer(regex_str, input_str, overlapped=True)]
    return regex_indexes
