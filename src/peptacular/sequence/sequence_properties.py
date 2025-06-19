"""
This section is largely based on the code from the `biopython` library, specifically the `Bio.SeqUtils` module.
"""

from typing import List, Literal, Union, Optional, Dict
import math

from .sequence_funcs import get_annotation_input, count_residues, percent_residues
from ..proforma.proforma_parser import ProFormaAnnotation
from .ProtParamData import *
from .weights import get_weights


def aromaticity(
    sequence: Union[str, ProFormaAnnotation],
    aromatic_residues: str = "YWF",
    precision: Optional[float] = None,
) -> float:
    """
    Calculates the aromaticity value of a protein according to Lobry, 1994.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: Union[str, ProFormaAnnotation]
    :param aromatic_residues: String containing the amino acids considered aromatic (default: "YWF")
    :type aromatic_residues: str
    :raises ValueError: If the input sequence contains multiple sequences
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The aromaticity value (relative frequency of Phe+Trp+Tyr)
    :rtype: float

    .. code-block:: python

        >>> aromaticity('PEPTIDE', precision=2)
        0.0
        >>> aromaticity('PEPTIDEYWY', precision=2)
        0.3
    """
    residue_perc = percent_residues(sequence, precision=2)
    aromaticity = sum(residue_perc.get(aa, 0) / 100 for aa in aromatic_residues)

    if precision is not None:
        aromaticity = round(aromaticity, precision)

    return aromaticity


def get_aa_value(
    aa: str,
    data: Dict[str, float],
    default: Literal["zero", "avg", "min", "max", "median", "error", "skip"],
    weight: float = 1,
    normalize: bool = False,
) -> Optional[float]:
    if default == "zero":
        default_value = 0.0
    elif default == "avg":
        default_value = sum(data.values()) / len(data)
    elif default == "min":
        default_value = min(data.values())
    elif default == "max":
        default_value = max(data.values())
    elif default == "median":
        sorted_values = sorted(data.values())
        mid = len(sorted_values) // 2
        if len(sorted_values) % 2 == 0:
            default_value = (sorted_values[mid - 1] + sorted_values[mid]) / 2
        else:
            default_value = sorted_values[mid]
    elif default == "error":
        default_value = "error"
    elif default == "skip":
        return "skip"
    else:
        raise ValueError(
            f"Invalid default value: {default}. Choose from 'zero', 'avg', 'min', 'max', 'median', 'error', or 'skip'."
        )

    # B -> Aspartic acid or Asparagine
    # J -> Leucine or Isoleucine
    # Z -> Glutamic acid or Glutamine
    # X -> Any amino acid (not a valid single-letter code)

    if aa in data:
        if normalize:
            # Normalize the value to a range of 0-1
            min_value = min(data.values())
            max_value = max(data.values())
            return ((data[aa] - min_value) / (max_value - min_value)) * weight
        else:
            # Return the raw value multiplied by the weight
            return data[aa] * weight
    elif aa == "B":
        # average of Aspartic acid and Asparagine
        return (
            (get_aa_value("D", data, default) + get_aa_value("N", data, default)) / 2
        ) * weight
    elif aa == "J":
        # average of Leucine and Isoleucine
        return (
            (get_aa_value("L", data, default) + get_aa_value("I", data, default)) / 2
        ) * weight
    elif aa == "Z":
        # average of Glutamic acid and Glutamine
        return (
            (get_aa_value("E", data, default) + get_aa_value("Q", data, default)) / 2
        ) * weight
    elif aa == "X":
        # Mean of all
        if normalize:
            # Normalize the average value to a range of 0-1
            min_value = min(data.values())
            max_value = max(data.values())
            return (
                (sum(data.values()) / len(data) - min_value) / (max_value - min_value)
            ) * weight
        else:
            # Return the average value multiplied by the weight
            return (sum(data.values()) / len(data)) * weight
    else:
        if default_value == "error":
            raise ValueError(
                f"Invalid amino acid: {aa}. No hydrophobicity value found."
            )
        elif default_value == "skip":
            return None
        return default_value * weight


def calc_property(
    sequence: Union[str, ProFormaAnnotation],
    scale_name: str,
    default: Literal["zero", "avg", "min", "max", "median", "error", "skip"] = "error",
    normalization: Literal["sum", "avg"] = "avg",
    normalize_scale: bool = False,
    weights: Union[
        list[float],
        Literal[
            "uniform",
            "linear",
            "exponential",
            "gaussian",
            "sigmoid",
            "cosine",
            "sinusoidal",
        ],
    ] = "uniform",
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    precision: Optional[float] = None,
) -> float:
    """
    Generic function to calculate the average property value of a protein sequence.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :param scale_name: The name of the scale to use
    :param scales_dict: Dictionary containing the available scales
    :param default: How to handle missing amino acids
    :param precision: Number of decimal places to round to
    :return: The average property value
    """
    if scale_name not in hydrophobicity_scales:
        raise ValueError(
            f"Invalid scale: {scale_name}. Choose from {list(hydrophobicity_scales.keys())}."
        )

    annotation = get_annotation_input(sequence, copy=False)

    weights = get_weights(
        length=len(annotation),
        weights=weights,
        min_weight=min_weight,
        max_weight=max_weight,
    )
    values = []
    for i, aa in enumerate(annotation.stripped_sequence):
        val = get_aa_value(
            aa=aa,
            data=hydrophobicity_scales[scale_name],
            default=default,
            weight=weights[i],
            normalize=normalize_scale,
        )
        if isinstance(val, (int, float)):
            values.append(val)

    if normalization == "sum":
        result = sum(values) if values else 0.0
    elif normalization == "avg":
        result = sum(values) / len(values) if values else 0.0
    else:
        raise ValueError(
            f"Invalid normalization method: {normalization}. Choose 'sum' or 'avg'."
        )

    if precision is not None:
        result = round(result, precision)

    return result


def calc_window_property(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "KyteDoolitle",
    window_size: int = 9,
    default: Literal["zero", "avg", "min", "max", "median", "error", "skip"] = "error",
    normalization: Literal["sum", "avg"] = "avg",
    normalize_scale: bool = False,
    weights: Union[
        list[float],
        Literal[
            "uniform",
            "linear",
            "exponential",
            "gaussian",
            "sigmoid",
            "cosine",
            "sinusoidal",
        ],
    ] = "uniform",
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    precision: Optional[float] = None,
) -> List[float]:

    annotation = get_annotation_input(sequence, copy=False)

    weights = get_weights(
        window_size, weights=weights, min_weight=min_weight, max_weight=max_weight
    )

    l = []
    for i, window_sequence in enumerate(annotation.sliding_windows(window_size)):
        if len(window_sequence) != window_size:
            raise ValueError(
                f"Window size {window_size} does not match sequence length {len(window_sequence)}."
            )

        # Calculate the property average for the current window
        window_average = calc_property(
            sequence=window_sequence,
            scale_name=scale,
            default=default,
            normalization=normalization,
            normalize_scale=normalize_scale,
            weights=weights,
            min_weight=min_weight,
            max_weight=max_weight,
            precision=precision,
        )
        l.append(window_average)

    return l


def hydrophobicity(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "KyteDoolitle",
    default: Literal["zero", "avg", "min", "max", "median", "error", "skip"] = "error",
    normalization: Literal["sum", "avg"] = "avg",
    precision: Optional[float] = None,
) -> float:
    """
    Calculates the hydrophobicity value of a protein based on Cowan's hydrophobicity indices.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: Union[str, ProFormaAnnotation]
    :param scale: The hydrophobicity scale to use (default: 'KyteDoolitle')
    :type scale: str
    :param default: How to handle missing amino acids (default: 'error')
    :type default: Literal['zero', 'avg', 'min', 'max', 'median', 'error', 'skip']
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :raises ValueError: If the input sequence contains multiple sequences
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The hydrophobicity value
    :rtype: float
    """
    return calc_property(
        sequence=sequence,
        scale=scale,
        default=default,
        normalization=normalization,
        precision=precision,
    )


def flexibility_vihinen(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "Vihinen",
    default: Literal["zero", "avg", "min", "max", "median", "error", "skip"] = "error",
    precision: Optional[float] = None,
) -> float:
    """
    Calculates the flexibility value of a protein based on normalized flexibility
    parameters (B-values) according to Vihinen et al., 1994.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: Union[str, ProFormaAnnotation]
    :param scale: The flexibility scale to use (default: 'Vihinen')
    :type scale: str
    :param default: How to handle missing amino acids (default: 'error')
    :type default: Literal['zero', 'avg', 'min', 'max', 'median', 'error', 'skip']
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :raises ValueError: If the input sequence contains multiple sequences or invalid scale
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The flexibility value
    :rtype: float
    """
    return calc_property(sequence, scale, flexibility_scales, default, precision)


def hydrophilicity(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "HoppWood",
    default: Literal["zero", "avg", "min", "max", "median", "error", "skip"] = "error",
    normalization: Literal["sum", "avg"] = "avg",
    precision: Optional[float] = None,
) -> float:
    """
    Calculates the hydrophilicity value of a protein based on Hopp & Wood
    hydrophilicity scale (Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981)).

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: Union[str, ProFormaAnnotation]
    :param scale: The hydrophilicity scale to use (default: 'HoppWood')
    :type scale: str
    :param default: How to handle missing amino acids (default: 'error')
    :type default: Literal['zero', 'avg', 'min', 'max', 'median', 'error', 'skip']
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :raises ValueError: If the input sequence contains multiple sequences or invalid scale
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The hydrophilicity value
    :rtype: float

    """
    return calc_property(
        sequence=sequence,
        scale=scale,
        default=default,
        normalization=normalization,
        precision=precision,
    )


def surface_accessibility(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "Emini",
    default: Literal["zero", "avg", "min", "max", "median", "error", "skip"] = "error",
    normalization: Literal["sum", "avg"] = "avg",
    precision: Optional[float] = None,
) -> float:
    """
    Calculates the surface accessibility value of a protein based on different scales.

    Available scales:
    - 'Emini': Emini Surface fractional probability (default)
    - 'Janin': Janin Interior to surface transfer energy scale

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: Union[str, ProFormaAnnotation]
    :param scale: The surface accessibility scale to use (default: 'Emini')
    :type scale: str
    :param default: How to handle missing amino acids (default: 'error')
    :type default: Literal['zero', 'avg', 'min', 'max', 'median', 'error', 'skip']
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :raises ValueError: If the input sequence contains multiple sequences or invalid scale
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The surface accessibility value
    :rtype: float
    """
    return calc_property(
        sequence=sequence,
        scale=scale,
        default=default,
        normalization=normalization,
        precision=precision,
    )


def charge_at_ph(
    sequence: Union[str, ProFormaAnnotation],
    pH: float = 7.0,
    precision: Optional[float] = None,
) -> float:
    """
    Calculate the charge of a protein at given pH using the Henderson-Hasselbalch equation.

    Uses updated amino acid pKa values with sequence-specific N-terminal and C-terminal pK values.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: Union[str, ProFormaAnnotation]
    :param pH: The pH at which to calculate the charge (default: 7.0)
    :type pH: float
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :raises ValueError: If the input sequence contains multiple sequences
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The net charge at the given pH
    :rtype: float
    """
    seq_str = get_annotation_input(sequence, copy=False).sequence

    # Count amino acids
    aa_counts = count_residues(sequence)

    # Get terminal residues
    nterm, cterm = seq_str[0], seq_str[-1]

    # Calculate positive charge (basic groups)
    positive_charge = 0.0

    # N-terminal charge
    nterm_pK = get_aa_value(aa=nterm, data=pk_nterminal, default="error")
    partial_charge = 1.0 / (10 ** (pH - nterm_pK) + 1.0)
    positive_charge += partial_charge

    # Side chain positive charges
    for aa in "KRH":
        count = float(aa_counts.get(aa, 0))
        if count > 0:
            pK = get_aa_value(aa=aa, data=pk_sidechain, default="error")
            partial_charge = 1.0 / (10 ** (pH - pK) + 1.0)
            positive_charge += count * partial_charge

    # Calculate negative charge (acidic groups)
    negative_charge = 0.0

    # C-terminal charge
    cterm_pK = get_aa_value(aa=cterm, data=pk_cterminal, default="error")
    partial_charge = 1.0 / (10 ** (cterm_pK - pH) + 1.0)
    negative_charge += partial_charge

    # Side chain negative charges
    for aa in "DECY":
        count = float(aa_counts.get(aa, 0))
        if count > 0:
            pK = get_aa_value(aa=aa, data=pk_sidechain, default="error")
            if pK > 0:  # Only calculate if pK exists (non-zero)
                partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
                negative_charge += count * partial_charge

    net_charge = positive_charge - negative_charge

    if precision is not None:
        net_charge = round(net_charge, precision)

    return net_charge


def pi(
    sequence: Union[str, ProFormaAnnotation], precision: Optional[float] = None
) -> float:
    """
    Calculate the isoelectric point (pI) of a protein using the bisection method.

    The isoelectric point is the pH at which the net charge of the protein is zero.
    Uses the Bjellqvist method with sequence-specific N-terminal and C-terminal pK values.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: Union[str, ProFormaAnnotation]
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :raises ValueError: If the input sequence contains multiple sequences
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The isoelectric point (pI)
    :rtype: float

    .. code-block:: python

        >>> pi('PEPTIDE', precision=2)
        4.05
        >>> pi('INGAR', precision=2)
        11.04
    """

    def _calculate_pi(
        pH: float = 7.775, min_: float = 4.05, max_: float = 12.0, tol_: float = 0.001
    ) -> float:
        """Recursive bisection method to find pI."""
        charge = charge_at_ph(sequence, pH)
        if max_ - min_ > tol_:
            if charge > 0.0:
                min_ = pH
            else:
                max_ = pH
            next_pH = (min_ + max_) / 2
            return _calculate_pi(next_pH, min_, max_, tol_)
        return pH

    isoelectric_point = _calculate_pi()

    if precision is not None:
        isoelectric_point = round(isoelectric_point, precision)

    return isoelectric_point
