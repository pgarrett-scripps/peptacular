"""
This module provides functions to generate permutations, combinations, and products of sequences.
"""

from typing import Generator, Optional, Union

from ..proforma.annotation import ProFormaAnnotation
from .util import get_annotation_input


def permutations(
    sequence: Union[str, ProFormaAnnotation],
    size: Optional[int] = None,
    include_plus: bool = False,
    precision: Optional[int] = None,
) -> Generator[str, None, None]:
    """
    Generates all permutations of the input sequence. Terminal sequence are kept in place.

    :param sequence: The sequence to be permuted.
    :type sequence: str
    :param size: The size of the permutations.
    :type size: int

    :return: A list of all permutations of the input sequence.
    :rtype: List[str]

    .. code-block:: python

        >>> list(permutations('PET'))
        ['PET', 'PTE', 'EPT', 'ETP', 'TPE', 'TEP']

        >>> list(permutations('[3]-PET-[1]'))
        ['[3]-PET-[1]', '[3]-PTE-[1]', '[3]-EPT-[1]', '[3]-ETP-[1]', '[3]-TPE-[1]', '[3]-TEP-[1]']

        >>> list(permutations('PE[3.14]T'))
        ['PE[3.14]T', 'PTE[3.14]', 'E[3.14]PT', 'E[3.14]TP', 'TPE[3.14]', 'TE[3.14]P']

        >>> list(permutations('<13C>PET'))
        ['<13C>PET', '<13C>PTE', '<13C>EPT', '<13C>ETP', '<13C>TPE', '<13C>TEP']

    """
    annotation = get_annotation_input(sequence, copy=False)
    return (
        a.serialize(include_plus=include_plus, precision=precision)
        for a in annotation.permutations(size=size)
    )


def product(
    sequence: Union[str, ProFormaAnnotation],
    repeat: Union[int, None],
    include_plus: bool = False,
    precision: Optional[int] = None,
) -> Generator[str, None, None]:
    """
    Generates all sartesian products of the input sequence of a given size. Terminal sequence are kept in place.

    :param sequence: The sequence to be combined.
    :type sequence: str
    :param repeat: The size of the combinations to be generated.
    :type repeat: int

    :return: A list of all combinations of the input sequence of the given size.
    :rtype: List[str]

    .. code-block:: python

        >>> list(product('PET', 2))
        ['PP', 'PE', 'PT', 'EP', 'EE', 'ET', 'TP', 'TE', 'TT']

        >>> list(product('[3]-PET-[1]', 2))[:5]
        ['[3]-PP-[1]', '[3]-PE-[1]', '[3]-PT-[1]', '[3]-EP-[1]', '[3]-EE-[1]']

        >>> list(product('PE[3.14]T', 2))
        ['PP', 'PE[3.14]', 'PT', 'E[3.14]P', 'E[3.14]E[3.14]', 'E[3.14]T', 'TP', 'TE[3.14]', 'TT']

        >>> list(product('<13C>PET', 2))
        ['<13C>PP', '<13C>PE', '<13C>PT', '<13C>EP', '<13C>EE', '<13C>ET', '<13C>TP', '<13C>TE', '<13C>TT']

    """

    annotation = get_annotation_input(sequence=sequence, copy=False)
    return (
        a.serialize(include_plus=include_plus, precision=precision)
        for a in annotation.product(repeat=repeat)
    )


def combinations(
    sequence: Union[str, ProFormaAnnotation],
    size: Union[int, None],
    include_plus: bool = False,
    precision: Optional[int] = None,
) -> Generator[str, None, None]:
    """
    Generates all combinations of the input sequence of a given size. Terminal sequence are kept in place.

    :param sequence: The sequence to be combined.
    :type sequence: str

    :param size: The size of the combinations to be generated.
    :type size: int

    :return: A list of all combinations of the input sequence of the given size.
    :rtype: List[str]

    .. code-block:: python

        >>> list(combinations('PET', 2))
        ['PE', 'PT', 'ET']

        >>> list(combinations('[3]-PET-[1]', 2))
        ['[3]-PE-[1]', '[3]-PT-[1]', '[3]-ET-[1]']

        >>> list(combinations('PE[3.14]T', 2))
        ['PE[3.14]', 'PT', 'E[3.14]T']

        >>> list(combinations('<13C>PET', 2))
        ['<13C>PE', '<13C>PT', '<13C>ET']

    """
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return (
        a.serialize(include_plus=include_plus, precision=precision)
        for a in annotation.combinations(r=size)
    )


def combinations_with_replacement(
    sequence: Union[str, ProFormaAnnotation],
    size: Union[int, None],
    include_plus: bool = False,
    precision: Optional[int] = None,
) -> Generator[str, None, None]:
    """
    Generates all combinations with replacement of the input sequence of a given size. Terminal sequence are kept
    in place.

    :param sequence: The sequence to be combined.
    :type sequence: str

    :param size: The size of the combinations to be generated.
    :type size: int

    :return: A list of all combinations of the input sequence of the given size.
    :rtype: List[str]

    .. code-block:: python

        >>> list(combinations_with_replacement('PET', 2))
        ['PP', 'PE', 'PT', 'EE', 'ET', 'TT']

        >>> list(combinations_with_replacement('[3]-PET-[1]', 2))
        ['[3]-PP-[1]', '[3]-PE-[1]', '[3]-PT-[1]', '[3]-EE-[1]', '[3]-ET-[1]', '[3]-TT-[1]']

        >>> list(combinations_with_replacement('PE[3.14]T', 2))
        ['PP', 'PE[3.14]', 'PT', 'E[3.14]E[3.14]', 'E[3.14]T', 'TT']

        >>> list(combinations_with_replacement('<13C>PET', 2))
        ['<13C>PP', '<13C>PE', '<13C>PT', '<13C>EE', '<13C>ET', '<13C>TT']

    """
    annotation = get_annotation_input(sequence=sequence, copy=False)
    return (
        a.serialize(include_plus=include_plus, precision=precision)
        for a in annotation.combinations_with_replacement(r=size)
    )
