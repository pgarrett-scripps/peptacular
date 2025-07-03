from typing import Dict, List, Literal, Optional, Union

from ..proforma.annot_properties import AGGREGATION_METHODS, MISSING_AA_HANDLING, WEIGHTING_SCHEMES

from .util import get_annotation_input

from ..proforma.annot import ProFormaAnnotation



def calc_property(
    sequence: Union[str, ProFormaAnnotation],
    scale: Union[str, Dict[str, float]],
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = False,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.calc_property(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
        precision=precision,
    )


def hydrophobicity(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "Kyte-Doolittle",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.hydrophobicity(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def flexibility(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "flexibility_vihinen",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.flexibility(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def hydrophilicity(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "hydrophilicity_hop_wood",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.hydrophilicity(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def surface_accessibility(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "surface_accessibility_vergoten",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.surface_accessibility(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def polarity(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "polarity_grantham",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.polarity(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def mutability(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "mutability",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.mutability(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def codons(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "codons",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "sum",
    normalize: bool = False,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.codons(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def bulkiness(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "bulkiness",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.bulkiness(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def recognition_factors(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "recognition_factors",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "sum",
    normalize: bool = False,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.recognition_factors(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def transmembrane_tendency(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "transmembrane_tendency",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.transmembrane_tendency(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def average_buried_area(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "average_buried_area",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.average_buried_area(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def hplc(
    sequence: Union[str, ProFormaAnnotation],
    scale: str = "hplc_meek_2_1",
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = True,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    precision: Optional[int] = None,
) -> float:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.hplc(
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        precision=precision,
    )


def calc_window_property(
    sequence: Union[str, ProFormaAnnotation],
    scale: Union[str, Dict[str, float]],
    window_size: int = 9,
    missing_aa_handling: MISSING_AA_HANDLING = "error",
    aggregation_method: AGGREGATION_METHODS = "avg",
    normalize: bool = False,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    precision: Optional[int] = None,
) -> List[float]:

    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.calc_window_property(
        scale=scale,
        window_size=window_size,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
        precision=precision,
    )


def charge_at_ph(
    sequence: Union[str, ProFormaAnnotation],
    pH: float = 7.0,
    precision: Optional[int] = None,
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

    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.charge_at_ph(
        pH=pH,
        precision=precision,
    )


def pi(
    sequence: Union[str, ProFormaAnnotation],
    precision: Optional[int] = None,
) -> float:
    """
    Calculate the isoelectric point (pI) of a protein sequence.

    Uses updated amino acid pKa values with sequence-specific N-terminal and C-terminal pK values.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: Union[str, ProFormaAnnotation]
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :raises ValueError: If the input sequence contains multiple sequences
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The isoelectric point (pI) of the sequence
    :rtype: float
    """

    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.pi(precision=precision)


def aa_property_percentage(
    sequence: Union[str, ProFormaAnnotation],
    residues: List[str],
    precision: Optional[int] = None,
) -> float:

    # Calculate the percentage of specified amino acid residues in a sequence.

    annotation = get_annotation_input(sequence=sequence, copy=False)
    return annotation.aa_property_percentage(
        residues=residues,
        precision=precision,
    )


def aromaticity(
    sequence: Union[str, ProFormaAnnotation],
    aromatic_residues: str = "YWF",
    precision: Optional[int] = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.aromaticity(
        aromatic_residues=aromatic_residues,
        precision=precision,
    )


def secondary_structure(
    sequence: Union[str, ProFormaAnnotation],
    scale: Literal["DeleageRoux", "Levitt", "ChouFasman"] = "DeleageRoux",
    precision: Optional[int] = None,
) -> Dict[str, Optional[float]]:
    # Calculate the secondary structure propensity of a peptide sequence.
    annotation = get_annotation_input(sequence=sequence, copy=True)
    return annotation.secondary_structure(
        scale=scale,
        precision=precision,
    )
