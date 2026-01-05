import warnings
from typing import List, Optional, Union

from ..annotation import ProFormaAnnotation


def sequence_to_annotation(sequence: str) -> ProFormaAnnotation:
    """
    Parses a peptide sequence with modifications and returns a ProFormaAnnotation object.
    """
    return ProFormaAnnotation.parse(sequence)


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


def get_annotation_input(
    sequence: str | ProFormaAnnotation,
    copy: bool = True,
) -> ProFormaAnnotation:
    if isinstance(sequence, str):
        return sequence_to_annotation(sequence)
    elif isinstance(sequence, ProFormaAnnotation):  # type: ignore
        return sequence.copy() if copy else sequence
    raise TypeError("Input sequence must be a string or ProFormaAnnotation object.")


def override_annotation_properties(
    annotation: ProFormaAnnotation,
    charge: Optional[int] = None,
    charge_adducts: Optional[Union[str | int | float, List[str | int | float]]] = None,
    isotope_mods: Optional[Union[str | int | float, List[str | int | float]]] = None,
):
    if annotation.has_charge and charge is not None:
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

    if annotation.has_charge and charge_adducts is not None:
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
