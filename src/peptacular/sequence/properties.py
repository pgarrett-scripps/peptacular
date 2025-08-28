from .util import get_annotation_input
from ..property.types import (
    AggregationMethodLiteral,
    MissingAAHandling,
    AggregationMethod,
    MissingAAHandlingLiteral,
    WeightingMethods,
    WeightingMethodsLiteral,
)
from ..property.properties import (
    HPLCScale,
    PhysicalPropertyScale,
    PolarityScale,
    SecondaryStructureMethod,
    SecondaryStructureType,
    SurfaceAccessibilityScale,
    HydrophobicityScale,
)
from ..property.core import (
    calc_property as _calc_property,
    calc_window_property as _calc_window_property,
    aa_property_percentage as _aa_property_percentage,
    charge_at_ph as _charge_at_ph,
    secondary_structure as _secondary_structure,
)
from ..utils2 import round_to_precision
from ..proforma.annotation import ProFormaAnnotation


def calc_property(
    sequence: str | ProFormaAnnotation,
    scale: str | dict[str, float],
    missing_aa_handling: (
        MissingAAHandlingLiteral | MissingAAHandling
    ) = MissingAAHandling.ERROR,
    aggregation_method: (
        AggregationMethodLiteral | AggregationMethod
    ) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (
        WeightingMethodsLiteral | WeightingMethods
    ) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=scale,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )
    return round_to_precision(val, precision)


def hydrophobicity(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=HydrophobicityScale.KYTE_DOOLITTLE,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def flexibility(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PhysicalPropertyScale.FLEXIBILITY_VIHINEN,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def hydrophilicity(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def surface_accessibility(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=SurfaceAccessibilityScale.VERGOTEN,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def polarity(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PolarityScale.GRANTHAM,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def mutability(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PhysicalPropertyScale.MUTABILITY,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def codons(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PhysicalPropertyScale.CODONS,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def bulkiness(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PhysicalPropertyScale.BULKINESS,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def recognition_factors(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PhysicalPropertyScale.RECOGNITION_FACTORS,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def transmembrane_tendency(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def average_buried_area(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=SurfaceAccessibilityScale.AVERAGE_BURIED_AREA,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def hplc(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=HPLCScale.MEEK_2_1,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def refractivity(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=PhysicalPropertyScale.REFRACTIVITY,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


def calc_window_property(
    sequence: str | ProFormaAnnotation,
    scale: str | dict[str, float],
    window_size: int = 9,
    missing_aa_handling: (
        MissingAAHandlingLiteral | MissingAAHandling
    ) = MissingAAHandling.ERROR,
    aggregation_method: (
        AggregationMethodLiteral | AggregationMethod
    ) = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: (
        WeightingMethodsLiteral | WeightingMethods
    ) = WeightingMethods.UNIFORM,
    min_weight: float = 0.1,
    max_weight: float = 1.0,
    precision: int | None = None,
) -> list[float]:
    annotation = get_annotation_input(sequence=sequence, copy=True)
    vals = _calc_window_property(
        sequence=annotation.stripped_sequence,
        scale=scale,
        window_size=window_size,
        missing_aa_handling=missing_aa_handling,
        aggregation_method=aggregation_method,
        normalize=normalize,
        weighting_scheme=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )
    if precision is not None:
        vals = [round_to_precision(v, precision) for v in vals]
    return vals


def charge_at_ph(
    sequence: str | ProFormaAnnotation,
    pH: float = 7.0,
    precision: int | None = None,
) -> float:
    """
    Calculate the charge of a protein at given pH using the Henderson-Hasselbalch equation.

    Uses updated amino acid pKa values with sequence-specific N-terminal and C-terminal pK values.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: str | ProFormaAnnotation
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
    val = _charge_at_ph(
        sequence=annotation.stripped_sequence,
        pH=pH,
    )
    return round_to_precision(val, precision)


def pi(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """
    Calculate the isoelectric point (pI) of a protein sequence.

    Uses updated amino acid pKa values with sequence-specific N-terminal and C-terminal pK values.

    :param sequence: The amino acid sequence or ProFormaAnnotation object
    :type sequence: str | ProFormaAnnotation
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :raises ValueError: If the input sequence contains multiple sequences
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The isoelectric point (pI) of the sequence
    :rtype: float
    """
    annotation = get_annotation_input(sequence=sequence, copy=False)

    def _calculate_pi(
        ph: float = 7.775,
        min_: float = 4.05,
        max_: float = 12.0,
        tol_: float = 0.001,
    ) -> float:
        """Recursive bisection method to find pI."""
        charge = _charge_at_ph(sequence=annotation.stripped_sequence, pH=ph)
        if max_ - min_ > tol_:
            if charge > 0.0:
                min_ = ph
            else:
                max_ = ph
            next_ph = (min_ + max_) / 2
            return _calculate_pi(next_ph, min_, max_, tol_)
        return ph

    isoelectric_point = _calculate_pi()
    return round_to_precision(isoelectric_point, precision)


def aa_property_percentage(
    sequence: str | ProFormaAnnotation,
    residues: list[str],
    precision: int | None = None,
) -> float:
    """
    Calculate the percentage of specified amino acid residues in a sequence.
    """
    annotation = get_annotation_input(sequence=sequence, copy=False)
    val = _aa_property_percentage(
        sequence=annotation.stripped_sequence,
        residues=residues,
    )
    return round_to_precision(val, precision)


def aromaticity(
    sequence: str | ProFormaAnnotation,
    aromatic_residues: list[str] = ["Y", "W", "F"],
    precision: int | None = None,
) -> float:
    """
    Calculate the aromaticity value of a protein according to Lobry, 1994.
    """
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _aa_property_percentage(
        sequence=annotation.stripped_sequence,
        residues=aromatic_residues,
    )
    return round_to_precision(val, precision)


def secondary_structure(
    sequence: str | ProFormaAnnotation,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    precision: int | None = None,
) -> dict[str, float]:
    """Calculate the secondary structure propensity of a peptide sequence."""
    annotation = get_annotation_input(sequence=sequence, copy=True)
    d = _secondary_structure(
        sequence=annotation.stripped_sequence,
        scale=scale,
    )
    if precision is not None:
        d = {k: round_to_precision(v, precision) for k, v in d.items()}
    return d


def alpha_helix_percent(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """
    Calculate the propensity for alpha helix formation using the Deleage-Roux scale.
    """
    d = secondary_structure(
        sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX, precision=precision
    )
    return d[SecondaryStructureType.ALPHA_HELIX]


def beta_sheet_percent(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """
    Calculate the propensity for beta sheet formation using the Deleage-Roux scale.
    """
    d = secondary_structure(
        sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX, precision=precision
    )
    return d[SecondaryStructureType.BETA_SHEET]


def beta_turn_percent(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """
    Calculate the propensity for beta turn formation using the Deleage-Roux scale.
    """
    d = secondary_structure(
        sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX, precision=precision
    )
    return d[SecondaryStructureType.BETA_TURN]


def coil_percent(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """
    Calculate the propensity for coil formation using the Deleage-Roux scale.
    """
    d = secondary_structure(
        sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX, precision=precision
    )
    return d[SecondaryStructureType.COIL]
