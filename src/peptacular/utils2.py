from typing import Optional, Tuple, Union

from .constants import BACKWARD_ION_TYPES, FORWARD_ION_TYPES, INTERNAL_ION_TYPES


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


def _pop_ion_count(ion: str) -> Tuple[int, str]:
    """
    Parse the charge of an ion from a string.

    :param ion: The ion string to parse.
    :type ion: str

    :return: A tuple containing the charge and the remaining ion string.
    :rtype: Tuple[int, str]

    .. code-block:: python

        >>> _pop_ion_count('+H+')
        (1, 'H+')

        >>> _pop_ion_count('+2Na+')
        (2, 'Na+')

        >>> _pop_ion_count('2I-')
        (2, 'I-')

        >>> _pop_ion_count('+e-')
        (1, 'e-')

        >>> _pop_ion_count('-2Na+')
        (-2, 'Na+')

        >>> _pop_ion_count('+2Mg2+')
        (2, 'Mg2+')

    """

    count_str = ""

    charge = 1
    for i, c in enumerate(ion):
        if c == "-":
            charge = -1
        elif c == "+":
            pass
        elif c.isdigit():
            count_str += c
        else:
            cnt = int(count_str) if count_str else 1
            return cnt * charge, ion[i:]

    raise ValueError(f"Bad Ion Count: {ion}")


def _pop_ion_symbol(ion: str) -> Tuple[str, str]:
    """
    Parse the symbol of an ion from a string.

    :param ion: The ion string to parse.
    :type ion: str

    :return: A tuple containing the symbol and the remaining ion string.
    :rtype: Tuple[str, str]

    .. code-block:: python

        >>> _pop_ion_symbol('H+')
        ('H', '+')

        >>> _pop_ion_symbol('Na+')
        ('Na', '+')

        >>> _pop_ion_symbol('I-')
        ('I', '-')

        >>> _pop_ion_symbol('e-')
        ('e', '-')

        >>> _pop_ion_symbol('Na+')
        ('Na', '+')

        >>> _pop_ion_symbol('Mg2+')
        ('Mg', '2+')

    """

    symbol = ""
    for i, c in enumerate(ion):
        if c.isdigit() or c in ["+", "-"]:
            return symbol, ion[i:]
        symbol += c
    return symbol, ""


def _pop_ion_charge(ion: str) -> Tuple[int, str]:
    """
    Parse the charge of an ion from a string.

    :param ion: The ion string to parse.
    :type ion: str

    :return: A tuple containing the charge and the remaining ion string.
    :rtype: Tuple[int, str]

    .. code-block:: python

        >>> _pop_ion_charge('+')
        (1, '')

        >>> _pop_ion_charge('+')
        (1, '')

        >>> _pop_ion_charge('-')
        (-1, '')

        >>> _pop_ion_charge('-')
        (-1, '')

        >>> _pop_ion_charge('+')
        (1, '')

        >>> _pop_ion_charge('2+')
        (2, '')

        >>> _pop_ion_charge('2-')
        (-2, '')

    """
    count_str, sign = "", 1
    for c in ion:
        if c == "-":
            sign = -1
        elif c == "+":
            continue
        else:
            count_str += c
    cnt = int(count_str) if count_str else 1
    return cnt * sign, ""


def parse_ion_elements(ion: str) -> Tuple[int, str, int]:
    """
    Parse the count, element, and charge of an ion from a string.

    :param ion: The ion string to parse.
    :type ion: str

    :return: A tuple containing the count, element, and charge of the ion.
    :rtype: Tuple[int, str, int]

    Examples:
        >>> parse_ion_elements('+H+')
        (1, 'H', 1)

        >>> parse_ion_elements('+2Na+')
        (2, 'Na', 1)

        >>> parse_ion_elements('2I-')
        (2, 'I', -1)

        >>> parse_ion_elements('+e-')
        (1, 'e', -1)

        >>> parse_ion_elements('-2Na+')
        (-2, 'Na', 1)

        >>> parse_ion_elements('+2Mg2+')
        (2, 'Mg', 2)

        >>> parse_ion_elements('+2Mg2-')
        (2, 'Mg', -2)

    """
    count, ion = _pop_ion_count(ion)
    symbol, charge = _pop_ion_symbol(ion)
    charge, _ = _pop_ion_charge(charge)
    return count, symbol, charge


class ModLabler:
    """
    A class to label modifications in a peptide sequence.

    :param mod: The modification to label.
    :type mod: str

    :param label: The label for the modification.
    :type label: str
    """

    def __init__(self):
        self.index = 0
        self.letter_index = 0
        self.curr_label = None

    def increment_num(self):
        """Generate a new numeric label based on the current index."""
        self.index += 1

    def increment_letter(self):
        """
        Generate a new letter label based on the current letter index.

        :return: The new letter label.
        :rtype: str
        """
        self.letter_index += 1

    @property
    def letter(self):
        """Generate a new letter label based on the current letter index."""
        # Convert index to letter(s): 0->a, 1->b, ..., 25->z, 26->aa, 27->ab, etc.
        result = ""
        temp_index = self.letter_index
        while True:
            result = chr(ord("a") + (temp_index % 26)) + result
            temp_index //= 26
            if temp_index == 0:
                break
            temp_index -= 1

        return result

    @property
    def num(self):
        """
        Generate a new numeric label based on the current index.

        :return: The new numeric label.
        :rtype: str
        """
        return str(self.index + 1)

    @property
    def label(self):
        """
        Get the current label and increment the index.

        :return: The current label.
        :rtype: str
        """
        return f"{self.letter}{self.num}"


def get_label(
    ion_type: str, charge: int, number: str, loss: float, isotope: int
) -> str:
    """
    Returns the label of the fragment, e.g., b2, y3i, etc.

    :return: Label of the fragment.
    :rtype: str
    """

    return (
        f"{'+' * charge}"
        f"{ion_type}"
        f"{number}"
        f"{'(' + str(loss) + ')' if loss != 0.0 else ''}"
        f"{'*' * isotope if isotope > 0 else ''}"
    )


def get_number(ion_type: str, len_sequence: int, start: int, end: int) -> str:
    """
    Returns the number of the fragment, e.g., 2 for b2, 3 for y3, etc.

    :return: Number of the fragment.
    :rtype: str
    """

    if ion_type in FORWARD_ION_TYPES:
        number = end
    elif ion_type in BACKWARD_ION_TYPES:
        number = len_sequence - start
    elif ion_type in INTERNAL_ION_TYPES:
        number = f"{start}-{end}"
    elif ion_type == "i":
        number = start
    else:
        raise ValueError("Wrong Ion Type")

    return str(number)


def round_to_precision(value: float, precision: Optional[int] = None) -> float:
    """
    Round a float to a specified number of decimal places.

    :param value: The float value to round.
    :type value: float
    :param precision: The number of decimal places to round to.
    :type precision: int

    :return: The rounded float value.
    :rtype: float
    """
    if precision is not None:
        value = round(value, precision)
    return value
