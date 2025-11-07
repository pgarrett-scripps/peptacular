"""Mixin class providing property calculation methods for sequences."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from peptacular.property.data import AROMATIC_AMINO_ACIDS

from ..funcs import round_to_precision
from .core import (
    aa_property_percentage,
    calc_property,
    calc_window_property,
    charge_at_ph,
    generate_sliding_window_features,
    secondary_structure,
)
from .properties import (
    HPLCScale,
    HydrophobicityScale,
    PhysicalPropertyScale,
    PolarityScale,
    SecondaryStructureMethod,
    SecondaryStructureType,
    SurfaceAccessibilityScale,
)
from .types import (
    AggregationMethod,
    AggregationMethodLiteral,
    MissingAAHandling,
    MissingAAHandlingLiteral,
    SequenceProtocol,
    WeightingMethods,
    WeightingMethodsLiteral,
)


class SequencePropertyMixin:
    """Mixin providing named property calculation methods"""

    __slots__ = ()

    def calc_property(
        self: SequenceProtocol,
        scale: str | dict[str, float],
        missing_aa_handling: (
            MissingAAHandlingLiteral | MissingAAHandling
        ) = MissingAAHandling.ERROR,
        aggregation_method: (
            AggregationMethodLiteral | AggregationMethod
        ) = AggregationMethod.AVG,
        normalize: bool = False,
        weighting_scheme: (
            WeightingMethodsLiteral | WeightingMethods | Sequence[float]
        ) = WeightingMethods.UNIFORM,
        min_weight: float = 0.1,
        max_weight: float = 1.0,
        precision: int | None = None,
    ) -> float:
        val = calc_property(
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

    @property
    def hydrophobicity(self: SequenceProtocol) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=HydrophobicityScale.KYTE_DOOLITTLE,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def flexibility(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PhysicalPropertyScale.FLEXIBILITY_VIHINEN,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def hydrophilicity(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PhysicalPropertyScale.HYDROPHILICITY_HOP_WOOD,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def surface_accessibility(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=SurfaceAccessibilityScale.VERGOTEN,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def polarity(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PolarityScale.GRANTHAM,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def mutability(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PhysicalPropertyScale.MUTABILITY,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def codons(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PhysicalPropertyScale.CODONS,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def bulkiness(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PhysicalPropertyScale.BULKINESS,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def recognition_factors(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PhysicalPropertyScale.RECOGNITION_FACTORS,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def transmembrane_tendency(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PhysicalPropertyScale.TRANSMEMBRANE_TENDENCY,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def average_buried_area(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=SurfaceAccessibilityScale.AVERAGE_BURIED_AREA,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def hplc(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=HPLCScale.MEEK_2_1,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def refractivity(
        self: SequenceProtocol,
    ) -> float:
        val = calc_property(
            sequence=self.stripped_sequence,
            scale=PhysicalPropertyScale.REFRACTIVITY,
            missing_aa_handling=MissingAAHandling.ERROR,
            aggregation_method=AggregationMethod.AVG,
            normalize=True,
            weighting_scheme=WeightingMethods.UNIFORM,
        )
        return val

    @property
    def aromaticity(
        self: SequenceProtocol,
        aromatic_residues: Iterable[str] = AROMATIC_AMINO_ACIDS,
        precision: int | None = None,
    ) -> float:
        """
        Calculate the aromaticity value of a protein according to Lobry, 1994.
        """
        val = aa_property_percentage(
            sequence=self.stripped_sequence,
            residues=list(aromatic_residues),
        )
        return round_to_precision(val, precision)

    def calc_window_property(
        self: SequenceProtocol,
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
            WeightingMethodsLiteral | WeightingMethods | Sequence[float]
        ) = WeightingMethods.UNIFORM,
        min_weight: float = 0.1,
        max_weight: float = 1.0,
        precision: int | None = None,
    ) -> list[float]:
        """Calculate property values over sliding windows"""
        vals = calc_window_property(
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
        val = aa_property_percentage(
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
        val = charge_at_ph(
            sequence=self.stripped_sequence,
            pH=pH,
        )
        return round_to_precision(val, precision)

    @property
    def pi(self: SequenceProtocol) -> float:
        """Calculate isoelectric point using bisection method"""

        def _calculate_pi(
            ph: float = 7.775,
            min_: float = 4.05,
            max_: float = 12.0,
            tol_: float = 0.001,
        ) -> float:
            """Recursive bisection method to find pI."""
            charge = charge_at_ph(sequence=self.stripped_sequence, pH=ph)
            if max_ - min_ > tol_:
                if charge > 0.0:
                    min_ = ph
                else:
                    max_ = ph
                next_ph = (min_ + max_) / 2
                return _calculate_pi(next_ph, min_, max_, tol_)
            return ph

        isoelectric_point = _calculate_pi()
        return isoelectric_point

    def secondary_structure(
        self: SequenceProtocol,
        scale: str = SecondaryStructureMethod.DELEAGE_ROUX,
        precision: int | None = None,
    ) -> dict[str, float]:
        """Calculate secondary structure propensities"""

        d = secondary_structure(
            sequence=self.stripped_sequence,
            scale=scale,
        )
        if precision is not None:
            d = {k: round_to_precision(v, precision) for k, v in d.items()}
        return d

    @property
    def alpha_helix_percent(self: SequenceProtocol) -> float:
        """
        Calculate the propensity for alpha helix formation using the Deleage-Roux scale.
        """
        d = secondary_structure(
            sequence=self.stripped_sequence,
            scale=SecondaryStructureMethod.DELEAGE_ROUX,
        )
        return d[SecondaryStructureType.ALPHA_HELIX]

    @property
    def beta_sheet_percent(self: SequenceProtocol) -> float:
        """
        Calculate the propensity for beta sheet formation using the Deleage-Roux scale.
        """
        d = secondary_structure(
            sequence=self.stripped_sequence,
            scale=SecondaryStructureMethod.DELEAGE_ROUX,
        )
        return d[SecondaryStructureType.BETA_SHEET]

    @property
    def beta_turn_percent(self: SequenceProtocol) -> float:
        """
        Calculate the propensity for beta turn formation using the Deleage-Roux scale.
        """
        d = secondary_structure(
            sequence=self.stripped_sequence,
            scale=SecondaryStructureMethod.DELEAGE_ROUX,
        )
        return d[SecondaryStructureType.BETA_TURN]

    @property
    def coil_percent(self: SequenceProtocol) -> float:
        """
        Calculate the propensity for coil formation using the Deleage-Roux scale.
        """
        d = secondary_structure(
            sequence=self.stripped_sequence,
            scale=SecondaryStructureMethod.DELEAGE_ROUX,
        )
        return d[SecondaryStructureType.COIL]

    def generate_sliding_window_features(
        self: SequenceProtocol,
        scale: str | dict[str, float],
        num_windows: int = 5,
        aa_overlap: int = 0,
        missing_aa_handling: (
            MissingAAHandlingLiteral | MissingAAHandling
        ) = MissingAAHandling.AVG,
        aggregation_method: (
            AggregationMethodLiteral | AggregationMethod
        ) = AggregationMethod.AVG,
        precision: int | None = None,
    ) -> list[float]:
        """Generate sliding window features for the sequence."""
        vals = generate_sliding_window_features(
            sequence=self.stripped_sequence,
            scale=scale,
            num_windows=num_windows,
            aa_overlap=aa_overlap,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
        )
        if precision is not None:
            vals = [round_to_precision(v, precision) for v in vals]
        return vals
