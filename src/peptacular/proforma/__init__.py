"""
This module provides the main interface for the ProForma package.
"""

from typing import Literal, overload
from .dclasses import *
from .annotation import ProFormaAnnotation
from .multi_annotation import MultiProFormaAnnotation


@overload
def parse(sequence: str, sequence_type: Literal["single"]) -> ProFormaAnnotation: ...


@overload
def parse(
    sequence: str, sequence_type: Literal["multi"]
) -> MultiProFormaAnnotation: ...


@overload
def parse(
    sequence: str, sequence_type: Literal["both"]
) -> ProFormaAnnotation | MultiProFormaAnnotation: ...


@overload
def parse(
    sequence: str, sequence_type: Literal["single", "multi", "both"] = "single"
) -> ProFormaAnnotation: ...


def parse(
    sequence: str, sequence_type: Literal["single", "multi", "both"] = "single"
) -> ProFormaAnnotation | MultiProFormaAnnotation:
    """
    Parses a ProForma sequence string and returns its corresponding annotation object.

    Note that the function's behavior and the type of object returned depend on the structure of the input sequence.
    Single sequences result in ProFormaAnnotation objects, while multi-sequences result in MultiProFormaAnnotation
    objects.

    :param sequence: The sequence to parse.
    :type sequence: str

    :raises ProFormaFormatError: If the sequence is not valid.

    :return: Either a ProFormaAnnotation or a MultiProFormaAnnotation, based on the input
    :rtype: Union[ProFormaAnnotation, MultiProFormaAnnotation]

    .. python::

        Parsing a simple peptide sequence:
        >>> isinstance(parse('PEPTIDE'), ProFormaAnnotation)
        True

        Parsing a sequence with multiple peptides or complex modifications:
        >>> isinstance(parse('PEPTIDE+PEPTIDE'), MultiProFormaAnnotation)
        True


    """

    match sequence_type:
        case "single":
            return ProFormaAnnotation.parse(sequence)
        case "multi":
            return MultiProFormaAnnotation.parse(sequence)
        case "both":
            try:
                return ProFormaAnnotation.parse(sequence)
            except ValueError:
                return MultiProFormaAnnotation.parse(sequence)
        case _:
            raise ValueError(f"Unknown sequence type: {sequence_type}")


def serialize(
    annotation: ProFormaAnnotation | MultiProFormaAnnotation,
    include_plus: bool = False,
) -> str:
    """
    Serializes a ProForma annotation or multiple ProForma annotations into a single string representation.

    :param annotation: Either a ProFormaAnnotation or a MultiProFormaAnnotation.
    :type annotation: Union[ProFormaAnnotation, MultiProFormaAnnotation]

    :return: A string representation of the ProForma annotation.
    :rtype: str

    . python::

        Serializing a simple ProForma annotation:
        >>> serialize(ProFormaAnnotation(sequence='PEPTIDE'))
        'PEPTIDE'

        >>> pfa1 = ProFormaAnnotation(sequence='PEPTIDE')
        >>> pfa2 = ProFormaAnnotation(sequence='PEPTIDE')

        Serializing a MultiProFormaAnnotation with chimeric connections:
        >>> multi_annotation = MultiProFormaAnnotation([pfa1, pfa2], [False])
        >>> serialize(multi_annotation)
        'PEPTIDE+PEPTIDE'

        Serializing a MultiProFormaAnnotation with crosslink connections:
        >>> multi_annotation = MultiProFormaAnnotation([pfa1, pfa2], [True])
        >>> p = serialize(multi_annotation)
        >>> p == r'PEPTIDE\\\PEPTIDE'
        True

    """

    return annotation.serialize(include_plus)
