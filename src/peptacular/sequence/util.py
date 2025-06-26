from typing import *
import warnings

from ..errors import ProFormaFormatError
from ..proforma_dataclasses import Mod
from ..proforma.annot import (
    ProFormaAnnotation,
)
from ..proforma.multi_annotation import MultiProFormaAnnotation
from ..proforma.parser import parse


def sequence_to_annotation(sequence: str) -> ProFormaAnnotation:
    """
    Parses a peptide sequence with modifications and returns a ProFormaAnnotation object.

    :param sequence: The amino acid sequence.
    :type sequence: str

    :raises ValueError: If the input sequence contains multiple sequences.
    :raises ProFormaFormatError: if the proforma sequence is not valid

    :return: A ProFormaAnnotation object representing the input sequence.
    :rtype: ProFormaAnnotation

    """
    annotation = parse(sequence)

    if isinstance(annotation, MultiProFormaAnnotation):
        raise ValueError(f"Invalid sequence: {sequence}")

    return annotation


def get_annotation_input(
    sequence: Union[str, ProFormaAnnotation], copy: bool = True
) -> ProFormaAnnotation:
    if isinstance(sequence, str):
        annotation = sequence_to_annotation(sequence)
    elif isinstance(sequence, ProFormaAnnotation):
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
    charge_adducts: Optional[List[Mod]] = None,
    isotope_mods: Optional[List[Mod]] = None,
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

    if annotation.charge_adducts is not None and charge_adducts is not None:
        # warning
        warnings.warn(
            "Both 'annotation.charge_adducts' and 'charge_adducts' are provided. Using user provided 'charge_adducts'.",
            UserWarning,
        )
        annotation.charge_adducts = charge_adducts
    elif charge_adducts is not None:
        annotation.charge_adducts = charge_adducts
    else:
        pass

    if annotation.isotope_mods is not None and isotope_mods is not None:
        # warning
        warnings.warn(
            "Both 'annotation.isotope_mods' and 'isotope_mods' are provided. Using user provided 'isotope_mods'.",
            UserWarning,
        )
        annotation.isotope_mods = isotope_mods
    elif isotope_mods is not None:
        annotation.isotope_mods = isotope_mods
    else:
        pass


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
        except ProFormaFormatError:
            return False
    return True
