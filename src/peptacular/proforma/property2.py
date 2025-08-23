from __future__ import annotations
from collections.abc import Generator, Iterable
from typing import Dict, List, Optional, Protocol, Union, Literal

# Import the property calculation function and data
from ..property.weights import get_weights
from ..property.data import (
    all_property_scales,
    pk_nterminal,
    pk_cterminal,
    pk_sidechain,
    secondary_structure_scales_by_name,
)
from ..utils2 import round_to_precision

from enum import Enum
import statistics


class MissingAAHandling(str, Enum):
    """Strategy for handling missing amino acid values"""

    ZERO = "zero"
    AVG = "avg"
    MIN = "min"
    MAX = "max"
    MEDIAN = "median"
    ERROR = "error"
    SKIP = "skip"

class AggregationMethod(str, Enum):
    """Strategy for aggregating amino acid values"""

    SUM = "sum"
    AVG = "avg"

def _calculate_median(values: list[float]) -> float:
    """Calculate median of a list of values"""
    return statistics.median(values)


def _get_default_value(
    handling: MissingAAHandling, aa_data: dict[str, float], aa: str
) -> float | None:
    """Calculate default value based on missing AA handling strategy"""
    match handling:
        case MissingAAHandling.ZERO:
            return 0.0
        case MissingAAHandling.AVG:
            return sum(aa_data.values()) / len(aa_data)
        case MissingAAHandling.MIN:
            return min(aa_data.values())
        case MissingAAHandling.MAX:
            return max(aa_data.values())
        case MissingAAHandling.MEDIAN:
            return _calculate_median(list(aa_data.values()))
        case MissingAAHandling.ERROR:
            raise ValueError(f"Invalid amino acid: {aa}")
        case MissingAAHandling.SKIP:
            return None
        case _:
            valid_options = ", ".join([f"'{opt.value}'" for opt in MissingAAHandling])
            raise ValueError(
                f"Invalid missing_aa_handling: {handling}. Choose from {valid_options}"
            )


def _normalize_value(value: float, aa_data: Dict[str, float]) -> float:
    """Normalize value to 0-1 range based on min/max in aa_data"""
    min_value = min(aa_data.values())
    max_value = max(aa_data.values())
    if max_value == min_value:
        return 0.0  # Avoid division by zero
    return (value - min_value) / (max_value - min_value)


def _apply_weighting_and_normalization(
    value: float, aa_data: dict[str, float], weighting_scheme: float, normalize: bool
) -> float:
    """Apply normalization and weighting to a value"""
    if normalize:
        value = _normalize_value(value, aa_data)
    return value * weighting_scheme


def _get_ambiguous_aa_value(
    aa: str,
    constituent_aas: list[str],
    aa_data: dict[str, float],
    missing_aa_handling: MissingAAHandling,
    weighting_scheme: float,
    normalize: bool,
) -> float:
    """Get average value for ambiguous amino acids (B, J, Z)"""
    values: list[float] = []
    for constituent_aa in constituent_aas:
        val = _get_aa_value(constituent_aa, aa_data, missing_aa_handling, 1.0, False)
        if val is not None:
            values.append(val)

    if not values:
        raise ValueError(f"No valid values found for ambiguous amino acid {aa}")

    avg_value = sum(values) / len(values)
    return _apply_weighting_and_normalization(
        avg_value, aa_data, weighting_scheme, normalize
    )


def _generate_string_sliding_windows(
    sequence: str,
    window_size: int,
    reverse: bool = False,
) -> Generator[str, None, None]:
    """
    Generate sliding windows of a plain amino acid sequence.
    
    Args:
        sequence: The amino acid sequence string
        window_size: Size of each window
        reverse: If True, generate windows from right to left
        
    Yields:
        String subsequences representing each window
        
    Raises:
        ValueError: If window_size is invalid or sequence is empty
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")
    
    if window_size <= 0:
        raise ValueError("Window size must be positive")
    
    seq_len = len(sequence)
    if window_size > seq_len:
        raise ValueError(
            f"Window size {window_size} cannot be greater than sequence length {seq_len}."
        )

    if reverse:
        # Generate windows from right to left
        for start in range(seq_len - window_size, -1, -1):
            stop = start + window_size
            yield sequence[start:stop]
    else:
        # Generate windows from left to right
        for start in range(seq_len - window_size + 1):
            stop = start + window_size
            yield sequence[start:stop]

def _get_aa_value(
    aa: str,
    aa_data: dict[str, float],
    missing_aa_handling: MissingAAHandling,
    weighting_scheme: float = 1.0,
    normalize: bool = False,
) -> float | None:
    """
    Get amino acid value with configurable handling for missing values.

    Args:
        aa: Single letter amino acid code
        aa_data: Dictionary mapping amino acid codes to values
        missing_aa_handling: Strategy for handling missing amino acids
        weighting_scheme: Multiplier applied to final value
        normalize: Whether to normalize value to 0-1 range

    Returns:
        Amino acid value or None if skipped

    Raises:
        ValueError: If amino acid is invalid and handling is ERROR
    """

    # Handle standard amino acids
    if aa in aa_data:
        value = aa_data[aa]
        return _apply_weighting_and_normalization(
            value, aa_data, weighting_scheme, normalize
        )

    # Handle ambiguous amino acids
    ambiguous_mappings = {
        "B": ["D", "N"],  # Aspartic acid or Asparagine
        "J": ["L", "I"],  # Leucine or Isoleucine
        "Z": ["E", "Q"],  # Glutamic acid or Glutamine
    }

    if aa in ambiguous_mappings:
        return _get_ambiguous_aa_value(
            aa,
            ambiguous_mappings[aa],
            aa_data,
            missing_aa_handling,
            weighting_scheme,
            normalize,
        )

    # Handle X (any amino acid) - use average of all
    if aa == "X":
        avg_value = sum(aa_data.values()) / len(aa_data)
        return _apply_weighting_and_normalization(
            avg_value, aa_data, weighting_scheme, normalize
        )

    # Handle unknown amino acids
    default_value = _get_default_value(missing_aa_handling, aa_data, aa)
    if default_value is None:
        return None

    return _apply_weighting_and_normalization(
        default_value, aa_data, weighting_scheme, normalize
    )


# Type aliases
WEIGHTING_SCHEMES = Union[
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
]


def _calc_property(
    sequence: str,
    scale: str | dict[str, float],
    missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
    aggregation_method: AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> float:
    """
    Generic function to calculate the average property value of a protein sequence.
    """
    weighting_scheme = get_weights(
        length=len(sequence),
        weights=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )

    if isinstance(scale, str):
        if scale not in all_property_scales:
            raise ValueError(
                f"Scale '{scale}' not found in available property scales: {list(all_property_scales.keys())}"
            )
        aa_data = all_property_scales[scale]
    else:
        aa_data = scale

    values: List[float] = []
    for i, aa in enumerate(sequence):
        val = _get_aa_value(
            aa=aa,
            aa_data=aa_data,
            missing_aa_handling=missing_aa_handling,
            weighting_scheme=weighting_scheme[i],
            normalize=normalize,
        )
        if isinstance(val, (int, float)):
            values.append(val)

    if aggregation_method == AggregationMethod.SUM:
        result = sum(values) if values else 0.0
    elif aggregation_method == AggregationMethod.AVG:
        result = sum(values) / len(values) if values else 0.0
    else:
        raise ValueError(
            f"Invalid aggregation method: {aggregation_method}. Choose '{AggregationMethod.SUM}' or '{AggregationMethod.AVG}'."
        )

    return result


def _count_residues(
    sequence: str,
) -> dict[str, int]:
    """
    Counts the occurrences of each amino acid residue in the sequence.
    """
    counts: dict[str, int] = {}
    for aa in sequence:
        counts[aa] = counts.get(aa, 0) + 1
    return counts


def _percent_residues(
    sequence: str,
) -> dict[str, float]:
    """
    Calculates the percentage of each amino acid residue in the sequence.
    """
    total = len(sequence)
    if total == 0:
        return {}

    counts = _count_residues(sequence)
    return {
        aa: count / total for aa, count in counts.items()
    }


def _aa_property_percentage(
    sequence: str,
    residues: Iterable[str],
) -> float:
    """
    Calculates the percentage of specified amino acids in the sequence.
    """

    residue_perc = _percent_residues(sequence)
    val = sum(residue_perc.get(aa, 0) for aa in residues)

    return val



def _calc_window_property(
    sequence: str,
    scale: str | dict[str, float],
    window_size: int = 9,
    missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
    aggregation_method: AggregationMethod = AggregationMethod.AVG,
    normalize: bool = False,
    weighting_scheme: WEIGHTING_SCHEMES = "uniform",
    min_weight: float = 0.1,
    max_weight: float = 1.0,
) -> List[float]:
    """Calculate property values over sliding windows"""
    weighting_scheme = get_weights(
        window_size,
        weights=weighting_scheme,
        min_weight=min_weight,
        max_weight=max_weight,
    )

    results: List[float] = []
    for window_sequence in _generate_string_sliding_windows(sequence, window_size):
        if len(window_sequence) != window_size:
            raise ValueError(
                f"Window size {window_size} does not match sequence length {len(window_sequence)}."
            )

        # Calculate the property average for the current window
        window_value = _calc_property(
            sequence=window_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )
        results.append(window_value)

    return results



def _charge_at_ph(
    sequence: str,
    pH: float = 7.0,
) -> float:
    """Calculate net charge at given pH"""
    # Count amino acids
    aa_counts = _count_residues(sequence=sequence)

    # Get terminal residues
    if not sequence:
        return 0.0

    nterm, cterm = sequence[0], sequence[-1]

    # Calculate positive charge (basic groups)
    positive_charge = 0.0

    # N-terminal charge
    nterm_pK = _get_aa_value(
        aa=nterm, aa_data=pk_nterminal, missing_aa_handling=MissingAAHandling.ERROR
    )

    if not isinstance(nterm_pK, (int, float)):
        raise ValueError(
            f"Invalid pK value for N-terminal amino acid '{nterm}': {nterm_pK}"
        )

    partial_charge: float = 1.0 / (10 ** (pH - nterm_pK) + 1.0)
    positive_charge += partial_charge

    # Side chain positive charges
    for aa in "KRH":
        count = float(aa_counts.get(aa, 0))
        if count > 0:
            pK = _get_aa_value(
                aa=aa, aa_data=pk_sidechain, missing_aa_handling=MissingAAHandling.ERROR
            )

            if not isinstance(pK, (int, float)) or pK <= 0:
                raise ValueError(
                    f"Invalid pK value for side chain amino acid '{aa}': {pK}"
                )

            partial_charge = 1.0 / (10 ** (pH - pK) + 1.0)
            positive_charge += count * partial_charge

    # Calculate negative charge (acidic groups)
    negative_charge = 0.0

    # C-terminal charge
    cterm_pK = _get_aa_value(
        aa=cterm, aa_data=pk_cterminal, missing_aa_handling=MissingAAHandling.ERROR
    )

    if not isinstance(cterm_pK, (int, float)):
        raise ValueError(
            f"Invalid pK value for C-terminal amino acid '{cterm}': {cterm_pK}"
        )

    partial_charge = 1.0 / (10 ** (cterm_pK - pH) + 1.0)
    negative_charge += partial_charge

    # Side chain negative charges
    for aa in "DECY":
        count = float(aa_counts.get(aa, 0))
        if count > 0:
            pK = _get_aa_value(
                aa=aa, aa_data=pk_sidechain, missing_aa_handling=MissingAAHandling.ERROR
            )

            if not isinstance(pK, (int, float)) or pK <= 0:
                raise ValueError(
                    f"Invalid pK value for side chain amino acid '{aa}': {pK}"
                )

            if pK > 0:  # Only calculate if pK exists (non-zero)
                partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
                negative_charge += count * partial_charge

    net_charge = positive_charge - negative_charge
    return net_charge


class SequenceProtocol(Protocol):
    """Protocol defining the interface for objects with sequences"""

    @property
    def stripped_sequence(self) -> str:
        """The sequence without modifications (for property calculations)"""
        ...


class SequencePropertyMixin:
    """Mixin providing named property calculation methods"""

    def calc_property(
        self: SequenceProtocol,
        scale: str | dict[str, float],
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = False,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        min_weight: float = 0.1,
        max_weight: float = 1.0,
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
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
        self: SequenceProtocol,
        scale: str = "Kyte-Doolittle",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def flexibility(
        self: SequenceProtocol,
        scale: str = "flexibility_vihinen",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def hydrophilicity(
        self: SequenceProtocol,
        scale: str = "hydrophilicity_hop_wood",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def surface_accessibility(
        self: SequenceProtocol,
        scale: str = "surface_accessibility_vergoten",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def polarity(
        self: SequenceProtocol,
        scale: str = "polarity_grantham",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def mutability(
        self: SequenceProtocol,
        scale: str = "mutability",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def codons(
        self: SequenceProtocol,
        scale: str = "codons",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.SUM,
        normalize: bool = False,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def bulkiness(
        self: SequenceProtocol,
        scale: str = "bulkiness",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def recognition_factors(
        self: SequenceProtocol,
        scale: str = "recognition_factors",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.SUM,
        normalize: bool = False,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def transmembrane_tendency(
        self: SequenceProtocol,
        scale: str = "transmembrane_tendency",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def average_buried_area(
        self: SequenceProtocol,
        scale: str = "average_buried_area",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def hplc(
        self: SequenceProtocol,
        scale: str = "hplc_meek_2_1",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def refractivity(
        self: SequenceProtocol,
        scale: str = "refractivity",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = _calc_property(
            sequence=self.stripped_sequence,
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
        )
        return round_to_precision(val, precision)

    def aromaticity(
        self: SequenceProtocol,
        aromatic_residues: str = "YWF",
        precision: int | None = None,
    ) -> float:
        """
        Calculate the aromaticity value of a protein according to Lobry, 1994.
        """
        val = _aa_property_percentage(
            sequence=self.stripped_sequence,
            residues=list(aromatic_residues),
        )
        return round_to_precision(val, precision)

    def calc_window_property(
        self: SequenceProtocol,
        scale: str | dict[str, float],
        window_size: int = 9,
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = False,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        min_weight: float = 0.1,
        max_weight: float = 1.0,
        precision: int | None = None,
    ) -> list[float]:
        """Calculate property values over sliding windows"""
        vals = _calc_window_property(
            sequence=self.stripped_sequence,
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


    def aa_property_percentage(
        self: SequenceProtocol,
        residues: Iterable[str],
        precision: int | None = None,
    ) -> float:
        """
        Calculate the percentage of specified amino acids in the sequence.
        """
        val = _aa_property_percentage(
            sequence=self.stripped_sequence,
            residues=list(residues),
        )
        return round_to_precision(val, precision)

    def charge_at_ph(
        self: SequenceProtocol,
        pH: float = 7.0,
        precision: int | None = None,
    ) -> float:
        """Calculate net charge at given pH"""
        # Count amino acids
        val = _charge_at_ph(
            sequence=self.stripped_sequence,
            pH=pH,
        )
        return round_to_precision(val, precision)

    def pi(self: SequenceProtocol, precision: int | None = None) -> float:
        """Calculate isoelectric point using bisection method"""

        def _calculate_pi(
            ph: float = 7.775,
            min_: float = 4.05,
            max_: float = 12.0,
            tol_: float = 0.001,
        ) -> float:
            """Recursive bisection method to find pI."""
            charge = _charge_at_ph(sequence=self.stripped_sequence, pH=ph)
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

    def secondary_structure(
        self: SequenceProtocol,
        scale: Literal["DeleageRoux", "Levitt", "ChouFasman"] = "DeleageRoux",
        precision: int | None = None,
    ) -> Dict[str, Optional[float]]:
        """Calculate secondary structure propensities"""
        if scale not in secondary_structure_scales_by_name:
            raise ValueError(
                f"Scale '{scale}' not found in available secondary structure scales: {list(secondary_structure_scales_by_name.keys())}"
            )

        d: Dict[str, float] = {}
        for structure_scale_name, structure_scale in secondary_structure_scales_by_name[
            scale
        ].items():
            val = _calc_property(
                sequence=self.stripped_sequence,
                scale=structure_scale,
                missing_aa_handling=MissingAAHandling.ERROR,
                aggregation_method=AggregationMethod.AVG,
                normalize=True,
                weighting_scheme="uniform",
            )
            d[structure_scale_name] = round_to_precision(val, precision)

        # Normalize the values to sum to 1
        total = sum(d.values())
        if total > 0:
            for key in d:
                d[key] /= total

        if "coil" not in d:
            d["coil"] = None  # type: ignore

        return d  # type: ignore
