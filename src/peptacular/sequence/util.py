import warnings
from enum import StrEnum
from typing import List, Literal, Optional, Union, overload

from ..proforma.annotation import ProFormaAnnotation
from ..proforma.dclasses import MOD_VALUE_TYPES
from ..proforma.multi_annotation import (
    MultiProFormaAnnotation,
)


class SequenceType(StrEnum):
    SINGLE = "single"
    MULTI = "multi"
    BOTH = "both"


@overload
def sequence_to_annotation(
    sequence: str,
    sequence_type: Literal[SequenceType.SINGLE] = SequenceType.SINGLE,
) -> ProFormaAnnotation: ...


@overload
def sequence_to_annotation(
    sequence: str, sequence_type: Literal[SequenceType.MULTI]
) -> MultiProFormaAnnotation: ...


@overload
def sequence_to_annotation(
    sequence: str, sequence_type: Literal[SequenceType.BOTH]
) -> ProFormaAnnotation | MultiProFormaAnnotation: ...


def sequence_to_annotation(
    sequence: str, sequence_type: SequenceType = SequenceType.SINGLE
) -> ProFormaAnnotation | MultiProFormaAnnotation:
    """
    Parses a peptide sequence with modifications and returns a ProFormaAnnotation or MultiProFormaAnnotation object.

    :param sequence: The amino acid sequence.
    :type sequence: str
    :param sequence_type: The type of sequence to parse (SINGLE, MULTI, or BOTH).
    :type sequence_type: SequenceType
    :raises ValueError: If the input sequence contains multiple sequences when SINGLE is specified.
    :raises ProFormaFormatError: if the proforma sequence is not valid
    :return: A ProFormaAnnotation or MultiProFormaAnnotation object representing the input sequence.
    :rtype: Union[ProFormaAnnotation, MultiProFormaAnnotation]
    """
    if sequence_type == SequenceType.SINGLE:
        annotation = ProFormaAnnotation.parse(sequence)
    elif sequence_type == SequenceType.MULTI:
        annotation = MultiProFormaAnnotation.parse(sequence)
    else:  # SequenceType.BOTH
        try:
            annotation = ProFormaAnnotation.parse(sequence)
        except:
            annotation = MultiProFormaAnnotation.parse(sequence)
    return annotation


@overload
def get_annotation_input(
    sequence: str | ProFormaAnnotation | MultiProFormaAnnotation,
    sequence_type: Literal[SequenceType.SINGLE] = SequenceType.SINGLE,
    copy: bool = True,
) -> ProFormaAnnotation: ...


@overload
def get_annotation_input(
    sequence: str | ProFormaAnnotation | MultiProFormaAnnotation,
    sequence_type: Literal[SequenceType.MULTI],
    copy: bool = True,
) -> MultiProFormaAnnotation: ...


@overload
def get_annotation_input(
    sequence: str | ProFormaAnnotation | MultiProFormaAnnotation,
    sequence_type: Literal[SequenceType.BOTH],
    copy: bool = True,
) -> ProFormaAnnotation | MultiProFormaAnnotation: ...


def get_annotation_input(
    sequence: str | ProFormaAnnotation | MultiProFormaAnnotation,
    sequence_type: SequenceType = SequenceType.SINGLE,
    copy: bool = True,
) -> ProFormaAnnotation | MultiProFormaAnnotation:
    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence, sequence_type)
    elif isinstance(sequence, ProFormaAnnotation) and sequence_type in (
        SequenceType.SINGLE,
        SequenceType.BOTH,
    ):
        if copy:
            annotation = sequence.copy()
        else:
            annotation = sequence
    elif isinstance(sequence, MultiProFormaAnnotation) and sequence_type in (
        SequenceType.MULTI,
        SequenceType.BOTH,
    ):
        if copy:
            annotation = sequence.copy()
        else:
            annotation = sequence
    else:
        raise TypeError(f"Expected str or ProFormaAnnotation, got {type(sequence)}")
    return annotation


def override_annotation_properties(
    annotation: ProFormaAnnotation,
    charge: Optional[int] = None,
    charge_adducts: Optional[Union[MOD_VALUE_TYPES, List[MOD_VALUE_TYPES]]] = None,
    isotope_mods: Optional[Union[MOD_VALUE_TYPES, List[MOD_VALUE_TYPES]]] = None,
):
    if annotation.charge is not None and charge is not None:
        # warning
        warnings.warn(
            "Both 'annotation.charge' and 'charge' are provided. Using user provided 'charge'.",
            UserWarning,
        )
        annotation.charge = charge
    elif charge is not None:
        annotation.charge = charge
    else:
        pass

    if annotation.has_charge_adducts and charge_adducts is not None:
        # warning
        warnings.warn(
            "Both 'annotation.charge_adducts' and 'charge_adducts' are provided. Using user provided 'charge_adducts'.",
            UserWarning,
        )
        annotation.charge_adducts = charge_adducts  # type: ignore

    if annotation.has_isotope_mods and isotope_mods is not None:
        # warning
        warnings.warn(
            "Both 'annotation.isotope_mods' and 'isotope_mods' are provided. Using user provided 'isotope_mods'.",
            UserWarning,
        )
        annotation.isotope_mods = isotope_mods  # type: ignore


def is_sequence_valid(sequence: Union[str, ProFormaAnnotation]) -> bool:
    """
    Checks if the input sequence is a valid ProForma sequence.

    :param sequence: The sequence or ProFormaAnnotation object to be validated.
    :type sequence: Union[str, ProFormaAnnotation]

    :return: True if the sequence is a valid ProForma sequence, False otherwise.
    :rtype: bool

    """

    if isinstance(sequence, str):
        try:
            _ = sequence_to_annotation(sequence)
        except ValueError:
            return False
    return True
