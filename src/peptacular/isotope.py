from typing import List, Union


def get_isotopes(sequence: str, n: Union[int, float]) -> List[float]:
    """
    Gets the abundance of N isotopes for a given sequence. If N is a float, the function will return the abundance of
    the first N isotopes up to the total abundance.

    :param sequence: The sequence.
    :type sequence: str
    :param n: The number of isotopes to return, or the total abundance.
    :type n: Union[int, float]

    :return: The abundance of N isotopes.
    :rtype: List[float]
    """
    pass


def get_isotope(sequence: str, n: int) -> float:
    """
    Gets the abundance of the Nth isotope for a given sequence.

    :param sequence: The sequence.
    :type sequence: str
    :param n: The isotope number.
    :type n: int

    :return: The abundance of the Nth isotope.
    :rtype: float
    """
    pass
