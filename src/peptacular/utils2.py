from typing import Tuple, Union


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
