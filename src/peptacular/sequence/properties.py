from collections.abc import Sequence
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
from ..funcs import round_to_precision
from ..proforma.annotation import ProFormaAnnotation
from .parrallel import parallel_apply_internal
from typing import overload, Literal


def _calc_property_single(
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
    """Calculate property for a single sequence"""
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


@overload
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
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def calc_property(
    sequence: Sequence[str | ProFormaAnnotation],
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
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def calc_property(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
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
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate a physicochemical property for a sequence or list of sequences.

    Automatically uses parallel processing when a list of sequences is provided.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _calc_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
            precision=precision,
        )
    else:
        return _calc_property_single(
            sequence=sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
            precision=precision,
        )


# Helper function for simple property calculations
def _simple_property_single(
    sequence: str | ProFormaAnnotation,
    scale: str | dict[str, float],
    precision: int | None = None,
) -> float:
    """Calculate a simple property for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=True)
    val = _calc_property(
        sequence=annotation.stripped_sequence,
        scale=scale,
        missing_aa_handling=MissingAAHandling.ERROR,
        aggregation_method=AggregationMethod.AVG,
        normalize=True,
        weighting_scheme=WeightingMethods.UNIFORM,
    )
    return round_to_precision(val, precision)


@overload
def hydrophobicity(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def hydrophobicity(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def hydrophobicity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate hydrophobicity. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=HydrophobicityScale.KYTE_DOOLITTLE,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=HydrophobicityScale.KYTE_DOOLITTLE,
            precision=precision,
        )


@overload
def flexibility(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def flexibility(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def flexibility(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate flexibility. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.FLEXIBILITY_VIHINEN,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.FLEXIBILITY_VIHINEN,
            precision=precision,
        )


@overload
def hydrophilicity(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def hydrophilicity(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def hydrophilicity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate hydrophilicity. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD,
            precision=precision,
        )


@overload
def surface_accessibility(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def surface_accessibility(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def surface_accessibility(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate surface accessibility. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=SurfaceAccessibilityScale.VERGOTEN,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=SurfaceAccessibilityScale.VERGOTEN,
            precision=precision,
        )


@overload
def polarity(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def polarity(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def polarity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate polarity. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PolarityScale.GRANTHAM,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PolarityScale.GRANTHAM,
            precision=precision,
        )


@overload
def mutability(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def mutability(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def mutability(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate mutability. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.MUTABILITY,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.MUTABILITY,
            precision=precision,
        )


@overload
def codons(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def codons(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def codons(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate codons. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.CODONS,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.CODONS,
            precision=precision,
        )


@overload
def bulkiness(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def bulkiness(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def bulkiness(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate bulkiness. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.BULKINESS,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.BULKINESS,
            precision=precision,
        )


@overload
def recognition_factors(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def recognition_factors(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def recognition_factors(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate recognition factors. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.RECOGNITION_FACTORS,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.RECOGNITION_FACTORS,
            precision=precision,
        )


@overload
def transmembrane_tendency(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def transmembrane_tendency(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def transmembrane_tendency(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate transmembrane tendency. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY,
            precision=precision,
        )


@overload
def average_buried_area(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def average_buried_area(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def average_buried_area(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate average buried area. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=SurfaceAccessibilityScale.AVERAGE_BURIED_AREA,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=SurfaceAccessibilityScale.AVERAGE_BURIED_AREA,
            precision=precision,
        )


@overload
def hplc(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def hplc(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def hplc(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate HPLC retention. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=HPLCScale.MEEK_2_1,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=HPLCScale.MEEK_2_1,
            precision=precision,
        )


@overload
def refractivity(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def refractivity(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def refractivity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """Calculate refractivity. Supports parallel processing for lists."""
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _simple_property_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=PhysicalPropertyScale.REFRACTIVITY,
            precision=precision,
        )
    else:
        return _simple_property_single(
            sequence=sequence,
            scale=PhysicalPropertyScale.REFRACTIVITY,
            precision=precision,
        )


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


def _charge_at_ph_single(
    sequence: str | ProFormaAnnotation,
    pH: float = 7.0,
    precision: int | None = None,
) -> float:
    """Calculate charge at pH for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=False)
    val = _charge_at_ph(
        sequence=annotation.stripped_sequence,
        pH=pH,
    )
    return round_to_precision(val, precision)


@overload
def charge_at_ph(
    sequence: str | ProFormaAnnotation,
    pH: float = 7.0,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def charge_at_ph(
    sequence: Sequence[str | ProFormaAnnotation],
    pH: float = 7.0,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def charge_at_ph(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    pH: float = 7.0,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate the charge of a protein at given pH using the Henderson-Hasselbalch equation.

    Uses updated amino acid pKa values with sequence-specific N-terminal and C-terminal pK values.
    Supports parallel processing for lists of sequences.

    :param sequence: The amino acid sequence, ProFormaAnnotation object, or list of sequences
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param pH: The pH at which to calculate the charge (default: 7.0)
    :type pH: float
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None
    :raises ValueError: If the input sequence contains multiple sequences
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The net charge at the given pH, or list of charges
    :rtype: float | list[float]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _charge_at_ph_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            pH=pH,
            precision=precision,
        )
    else:
        return _charge_at_ph_single(
            sequence=sequence,
            pH=pH,
            precision=precision,
        )


def _pi_single(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """Calculate pI for a single sequence"""
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


@overload
def pi(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def pi(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def pi(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate the isoelectric point (pI) of a protein sequence.

    Uses updated amino acid pKa values with sequence-specific N-terminal and C-terminal pK values.
    Supports parallel processing for lists of sequences.

    :param sequence: The amino acid sequence, ProFormaAnnotation object, or list of sequences
    :type sequence: str | ProFormaAnnotation | list[str | ProFormaAnnotation]
    :param precision: Number of decimal places to round to (default: None)
    :type precision: Optional[float]
    :param n_workers: Number of worker processes (only for lists). If None, uses CPU count.
    :type n_workers: int | None
    :param chunksize: Number of items per chunk (only for lists). If None, auto-calculated.
    :type chunksize: int | None
    :param method: 'process', 'thread', or None (auto-detect). Default is None.
    :type method: Literal["process", "thread"] | None
    :raises ValueError: If the input sequence contains multiple sequences
    :raises ProFormaFormatError: If the proforma sequence is not valid
    :return: The isoelectric point (pI) of the sequence, or list of pI values
    :rtype: float | list[float]
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _pi_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            precision=precision,
        )
    else:
        return _pi_single(
            sequence=sequence,
            precision=precision,
        )


def _aa_property_percentage_single(
    sequence: str | ProFormaAnnotation,
    residues: list[str],
    precision: int | None = None,
) -> float:
    """Calculate AA percentage for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=False)
    val = _aa_property_percentage(
        sequence=annotation.stripped_sequence,
        residues=residues,
    )
    return round_to_precision(val, precision)


@overload
def aa_property_percentage(
    sequence: str | ProFormaAnnotation,
    residues: list[str],
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def aa_property_percentage(
    sequence: Sequence[str | ProFormaAnnotation],
    residues: list[str],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def aa_property_percentage(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    residues: list[str],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate the percentage of specified amino acid residues in a sequence.
    Supports parallel processing for lists of sequences.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _aa_property_percentage_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            residues=residues,
            precision=precision,
        )
    else:
        return _aa_property_percentage_single(
            sequence=sequence,
            residues=residues,
            precision=precision,
        )


@overload
def aromaticity(
    sequence: str | ProFormaAnnotation,
    aromatic_residues: list[str] = ["Y", "W", "F"],
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def aromaticity(
    sequence: Sequence[str | ProFormaAnnotation],
    aromatic_residues: list[str] = ["Y", "W", "F"],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def aromaticity(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    aromatic_residues: list[str] = ["Y", "W", "F"],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate the aromaticity value of a protein according to Lobry, 1994.
    Supports parallel processing for lists of sequences.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _aa_property_percentage_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            residues=aromatic_residues,
            precision=precision,
        )
    else:
        return _aa_property_percentage_single(
            sequence=sequence,
            residues=aromatic_residues,
            precision=precision,
        )


def _secondary_structure_single(
    sequence: str | ProFormaAnnotation,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    precision: int | None = None,
) -> dict[str, float]:
    """Calculate secondary structure for a single sequence"""
    annotation = get_annotation_input(sequence=sequence, copy=True)
    d = _secondary_structure(
        sequence=annotation.stripped_sequence,
        scale=scale,
    )
    if precision is not None:
        d = {k: round_to_precision(v, precision) for k, v in d.items()}
    return d


@overload
def secondary_structure(
    sequence: str | ProFormaAnnotation,
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> dict[str, float]: ...


@overload
def secondary_structure(
    sequence: Sequence[str | ProFormaAnnotation],
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[dict[str, float]]: ...


def secondary_structure(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> dict[str, float] | list[dict[str, float]]:
    """
    Calculate the secondary structure propensity of a peptide sequence.
    Supports parallel processing for lists of sequences.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _secondary_structure_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            scale=scale,
            precision=precision,
        )
    else:
        return _secondary_structure_single(
            sequence=sequence,
            scale=scale,
            precision=precision,
        )


def _alpha_helix_percent_single(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """Calculate alpha helix percent for a single sequence"""
    d = _secondary_structure_single(
        sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX, precision=precision
    )
    return d[SecondaryStructureType.ALPHA_HELIX]


@overload
def alpha_helix_percent(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def alpha_helix_percent(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def alpha_helix_percent(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate the propensity for alpha helix formation using the Deleage-Roux scale.
    Supports parallel processing for lists of sequences.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _alpha_helix_percent_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            precision=precision,
        )
    else:
        return _alpha_helix_percent_single(
            sequence=sequence,
            precision=precision,
        )


def _beta_sheet_percent_single(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """Calculate beta sheet percent for a single sequence"""
    d = _secondary_structure_single(
        sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX, precision=precision
    )
    return d[SecondaryStructureType.BETA_SHEET]


@overload
def beta_sheet_percent(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def beta_sheet_percent(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def beta_sheet_percent(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate the propensity for beta sheet formation using the Deleage-Roux scale.
    Supports parallel processing for lists of sequences.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _beta_sheet_percent_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            precision=precision,
        )
    else:
        return _beta_sheet_percent_single(
            sequence=sequence,
            precision=precision,
        )


def _beta_turn_percent_single(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """Calculate beta turn percent for a single sequence"""
    d = _secondary_structure_single(
        sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX, precision=precision
    )
    return d[SecondaryStructureType.BETA_TURN]


@overload
def beta_turn_percent(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def beta_turn_percent(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def beta_turn_percent(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate the propensity for beta turn formation using the Deleage-Roux scale.
    Supports parallel processing for lists of sequences.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _beta_turn_percent_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            precision=precision,
        )
    else:
        return _beta_turn_percent_single(
            sequence=sequence,
            precision=precision,
        )


def _coil_percent_single(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
) -> float:
    """Calculate coil percent for a single sequence"""
    d = _secondary_structure_single(
        sequence, scale=SecondaryStructureMethod.DELEAGE_ROUX, precision=precision
    )
    return d[SecondaryStructureType.COIL]


@overload
def coil_percent(
    sequence: str | ProFormaAnnotation,
    precision: int | None = None,
    n_workers: None = None,
    chunksize: None = None,
    method: Literal["process", "thread"] | None = None,
) -> float: ...


@overload
def coil_percent(
    sequence: Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> list[float]: ...


def coil_percent(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    precision: int | None = None,
    n_workers: int | None = None,
    chunksize: int | None = None,
    method: Literal["process", "thread"] | None = None,
) -> float | list[float]:
    """
    Calculate the propensity for coil formation using the Deleage-Roux scale.
    Supports parallel processing for lists of sequences.
    """
    if (
        isinstance(sequence, Sequence)
        and not isinstance(sequence, str)
        and not isinstance(sequence, ProFormaAnnotation)
    ):
        return parallel_apply_internal(
            _coil_percent_single,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
            precision=precision,
        )
    else:
        return _coil_percent_single(
            sequence=sequence,
            precision=precision,
        )
