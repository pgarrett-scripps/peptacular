"""
This module provides functions to generate permutations, combinations, and products of sequences.
"""

from typing import List, Union

from peptacular.proforma.proforma_parser import ProFormaAnnotation
from peptacular.sequence.sequence_funcs import sequence_to_annotation


def permutations(sequence: Union[str, ProFormaAnnotation], size: int = None) -> List[str]:
    """
    Generates all permutations of the input sequence. Terminal sequence are kept in place.

    :param sequence: The sequence to be permuted.
    :type sequence: str
    :param size: The size of the permutations.
    :type size: int

    :return: A list of all permutations of the input sequence.
    :rtype: List[str]

    .. code-block:: python

        >>> permutations('PET')
        ['PET', 'PTE', 'EPT', 'ETP', 'TPE', 'TEP']

        >>> permutations('[3]-PET-[1]')
        ['[3]-PET-[1]', '[3]-PTE-[1]', '[3]-EPT-[1]', '[3]-ETP-[1]', '[3]-TPE-[1]', '[3]-TEP-[1]']

        >>> permutations('PE[3.14]T')
        ['PE[3.14]T', 'PTE[3.14]', 'E[3.14]PT', 'E[3.14]TP', 'TPE[3.14]', 'TE[3.14]P']

        >>> permutations('<13C>PET')
        ['<13C>PET', '<13C>PTE', '<13C>EPT', '<13C>ETP', '<13C>TPE', '<13C>TEP']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return [a.serialize() for a in annotation.permutations(size)]


def product(sequence: Union[str, ProFormaAnnotation], repeat: Union[int, None]) -> List[str]:
    """
    Generates all sartesian products of the input sequence of a given size. Terminal sequence are kept in place.

    :param sequence: The sequence to be combined.
    :type sequence: str
    :param repeat: The size of the combinations to be generated.
    :type repeat: int

    :return: A list of all combinations of the input sequence of the given size.
    :rtype: List[str]

    .. code-block:: python

        >>> product('PET', 2)
        ['PP', 'PE', 'PT', 'EP', 'EE', 'ET', 'TP', 'TE', 'TT']

        >>> product('[3]-PET-[1]', 2)[:5]
        ['[3]-PP-[1]', '[3]-PE-[1]', '[3]-PT-[1]', '[3]-EP-[1]', '[3]-EE-[1]']

        >>> product('PE[3.14]T', 2)
        ['PP', 'PE[3.14]', 'PT', 'E[3.14]P', 'E[3.14]E[3.14]', 'E[3.14]T', 'TP', 'TE[3.14]', 'TT']

        >>> product('<13C>PET', 2)
        ['<13C>PP', '<13C>PE', '<13C>PT', '<13C>EP', '<13C>EE', '<13C>ET', '<13C>TP', '<13C>TE', '<13C>TT']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return [a.serialize() for a in annotation.product(repeat)]


def combinations(sequence: Union[str, ProFormaAnnotation], size: Union[int, None]) -> List[str]:
    """
    Generates all combinations of the input sequence of a given size. Terminal sequence are kept in place.

    :param sequence: The sequence to be combined.
    :type sequence: str

    :param size: The size of the combinations to be generated.
    :type size: int

    :return: A list of all combinations of the input sequence of the given size.
    :rtype: List[str]

    .. code-block:: python

        >>> combinations('PET', 2)
        ['PE', 'PT', 'ET']

        >>> combinations('[3]-PET-[1]', 2)
        ['[3]-PE-[1]', '[3]-PT-[1]', '[3]-ET-[1]']

        >>> combinations('PE[3.14]T', 2)
        ['PE[3.14]', 'PT', 'E[3.14]T']

        >>> combinations('<13C>PET', 2)
        ['<13C>PE', '<13C>PT', '<13C>ET']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return [a.serialize() for a in annotation.combinations(size)]


def combinations_with_replacement(sequence: Union[str, ProFormaAnnotation], size: Union[int, None]) -> List[str]:
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

        >>> combinations_with_replacement('PET', 2)
        ['PP', 'PE', 'PT', 'EE', 'ET', 'TT']

        >>> combinations_with_replacement('[3]-PET-[1]', 2)
        ['[3]-PP-[1]', '[3]-PE-[1]', '[3]-PT-[1]', '[3]-EE-[1]', '[3]-ET-[1]', '[3]-TT-[1]']

        >>> combinations_with_replacement('PE[3.14]T', 2)
        ['PP', 'PE[3.14]', 'PT', 'E[3.14]E[3.14]', 'E[3.14]T', 'TT']

        >>> combinations_with_replacement('<13C>PET', 2)
        ['<13C>PP', '<13C>PE', '<13C>PT', '<13C>EE', '<13C>ET', '<13C>TT']

    """

    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    else:
        annotation = sequence

    return [a.serialize() for a in annotation.combinations_with_replacement(size)]
