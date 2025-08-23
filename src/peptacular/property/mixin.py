"""Mixin class providing property calculation methods for sequences."""

from __future__ import annotations
from collections.abc import Iterable
from typing import Dict, Optional, Literal

from .data import secondary_structure_scales_by_name
from .types import (
    MissingAAHandling,
    AggregationMethod, 
    WEIGHTING_SCHEMES,
    SequenceProtocol
)
from .core import (
    calc_property,
    calc_window_property,
    aa_property_percentage,
    charge_at_ph,
)
from ..utils2 import round_to_precision



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

    def hydrophobicity(
        self: SequenceProtocol,
        scale: str = "Kyte-Doolittle",
        missing_aa_handling: MissingAAHandling = MissingAAHandling.ERROR,
        aggregation_method: AggregationMethod = AggregationMethod.AVG,
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        precision: int | None = None,
    ) -> float:
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = calc_property(
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
        val = aa_property_percentage(
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

    def pi(self: SequenceProtocol, precision: int | None = None) -> float:
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
            val = calc_property(
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
