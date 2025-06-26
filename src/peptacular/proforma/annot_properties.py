"""
ProForma Parser Module
"""

from typing import *

from .annot_manipulation import ProFormaAnnotationManipulation
from .utils import get_aa_value
from ..weights import get_weights
from .property_data import (
    all_property_scales,
    pk_nterminal,
    pk_cterminal,
    pk_sidechain,
    secondary_structure_scales_by_name,
)


class ProFormaAnnotationProperty(ProFormaAnnotationManipulation):
    """
    A ProForma annotation with additional methods for serialization and parsing.
    This class is used to represent a ProForma annotation with a sequence, modifications, and other properties.
    """

    def calc_property(
        self,
        scale: Union[str, Dict[str, float]],
        missing_aa_handling: Literal[
            "zero", "avg", "min", "max", "median", "error", "skip"
        ] = "error",
        aggregation_method: Literal["sum", "avg"] = "avg",
        normalize: bool = False,
        weighting_scheme: Union[
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

        weighting_scheme = get_weights(
            length=len(self),
            weights=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )

        if isinstance(scale, str):

            if scale not in all_property_scales:
                raise ValueError(
                    f"Scale '{scale}' not found in available property scales: {all_property_scales.keys()}"
                )

            aa_data = all_property_scales[scale]
        else:
            aa_data = scale

        values = []
        for i, aa in enumerate(self.stripped_sequence):
            val = get_aa_value(
                aa=aa,
                aa_data=aa_data,
                missing_aa_handling=missing_aa_handling,
                weighting_scheme=weighting_scheme[i],
                normalize=normalize,
            )
            if isinstance(val, (int, float)):
                values.append(val)

        if aggregation_method == "sum":
            result = sum(values) if values else 0.0
        elif aggregation_method == "avg":
            result = sum(values) / len(values) if values else 0.0
        else:
            raise ValueError(
                f"Invalid normalization method: {aggregation_method}. Choose 'sum' or 'avg'."
            )

        if precision is not None:
            result = round(result, precision)

        return result

    def hydrophobicity(
        self,
        scale: str = "Kyte-Doolittle",
        missing_aa_handling: str = "error",
        aggregation_method: str = "avg",
        normalize: bool = True,
        weighting_scheme: str = "uniform",
        precision: Optional[int] = None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def flexibility(
        self,
        scale="flexibility_vihinen",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def hydrophilicity(
        self,
        scale="hydrophilicity_hop_wood",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def surface_accessibility(
        self,
        scale="surface_accessibility_vergoten",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def polarity(
        self,
        scale="polarity_grantham",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def mutability(
        self,
        scale="mutability",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def codons(
        self,
        scale="codons",
        missing_aa_handling="error",
        aggregation_method="sum",
        normalize=False,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def bulkiness(
        self,
        scale="bulkiness",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def recognition_factors(
        self,
        scale="recognition_factors",
        missing_aa_handling="error",
        aggregation_method="sum",
        normalize=False,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def transmembrane_tendency(
        self,
        scale="transmembrane_tendency",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def average_buried_area(
        self,
        scale="average_buried_area",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def hplc(
        self,
        scale="hplc_meek_2_1",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def refractivity(
        self,
        scale="refractivity",
        missing_aa_handling="error",
        aggregation_method="avg",
        normalize=True,
        weighting_scheme="uniform",
        precision=None,
    ) -> float:
        return self.calc_property(
            scale=scale,
            missing_aa_handling=missing_aa_handling,
            aggregation_method=aggregation_method,
            normalize=normalize,
            weighting_scheme=weighting_scheme,
            precision=precision,
        )

    def calc_window_property(
        self,
        scale: Union[str, Dict[str, float]],
        window_size: int = 9,
        missing_aa_handling: Literal[
            "zero", "avg", "min", "max", "median", "error", "skip"
        ] = "error",
        aggregation_method: Literal["sum", "avg"] = "avg",
        normalize: bool = False,
        weighting_scheme: Union[
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

        weighting_scheme = get_weights(
            window_size,
            weights=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )

        l = []
        for i, window_sequence in enumerate(self.sliding_windows(window_size)):
            if len(window_sequence) != window_size:
                raise ValueError(
                    f"Window size {window_size} does not match sequence length {len(window_sequence)}."
                )

            # Calculate the property average for the current window
            window_average = window_sequence.calc_property(
                scale=scale,
                missing_aa_handling=missing_aa_handling,
                aggregation_method=aggregation_method,
                normalize=normalize,
                weighting_scheme=weighting_scheme,
                min_weight=min_weight,
                max_weight=max_weight,
                precision=precision,
            )
            l.append(window_average)

        return l

    def charge_at_ph(
        self,
        pH: float = 7.0,
        precision: Optional[float] = None,
    ) -> float:
        # Count amino acids
        aa_counts = self.count_residues()

        # Get terminal residues
        nterm, cterm = self.sequence[0], self.sequence[-1]

        # Calculate positive charge (basic groups)
        positive_charge = 0.0

        # N-terminal charge
        nterm_pK = get_aa_value(
            aa=nterm, aa_data=pk_nterminal, missing_aa_handling="error"
        )
        partial_charge = 1.0 / (10 ** (pH - nterm_pK) + 1.0)
        positive_charge += partial_charge

        # Side chain positive charges
        for aa in "KRH":
            count = float(aa_counts.get(aa, 0))
            if count > 0:
                pK = get_aa_value(
                    aa=aa, aa_data=pk_sidechain, missing_aa_handling="error"
                )
                partial_charge = 1.0 / (10 ** (pH - pK) + 1.0)
                positive_charge += count * partial_charge

        # Calculate negative charge (acidic groups)
        negative_charge = 0.0

        # C-terminal charge
        cterm_pK = get_aa_value(
            aa=cterm, aa_data=pk_cterminal, missing_aa_handling="error"
        )
        partial_charge = 1.0 / (10 ** (cterm_pK - pH) + 1.0)
        negative_charge += partial_charge

        # Side chain negative charges
        for aa in "DECY":
            count = float(aa_counts.get(aa, 0))
            if count > 0:
                pK = get_aa_value(
                    aa=aa, aa_data=pk_sidechain, missing_aa_handling="error"
                )
                if pK > 0:  # Only calculate if pK exists (non-zero)
                    partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
                    negative_charge += count * partial_charge

        net_charge = positive_charge - negative_charge

        if precision is not None:
            net_charge = round(net_charge, precision)

        return net_charge

    def pi(self, precision: Optional[float] = None) -> float:
        def _calculate_pi(
            ph: float = 7.775,
            min_: float = 4.05,
            max_: float = 12.0,
            tol_: float = 0.001,
        ) -> float:
            """Recursive bisection method to find pI."""
            charge = self.charge_at_ph(ph)
            if max_ - min_ > tol_:
                if charge > 0.0:
                    min_ = ph
                else:
                    max_ = ph
                next_ph = (min_ + max_) / 2
                return _calculate_pi(next_ph, min_, max_, tol_)
            return ph

        isoelectric_point = _calculate_pi()

        if precision is not None:
            isoelectric_point = round(isoelectric_point, precision)

        return isoelectric_point

    def aa_property_percentage(
        self,
        residues: Iterable[str],
        precision: Optional[float] = None,
    ) -> Dict[str, float]:
        residue_perc = self.percent_residues(precision=None)
        val = sum(residue_perc.get(aa, 0) / 100 for aa in residues)

        if precision is not None:
            val = round(val, precision)

        return val

    def secondary_structure(
        self,
        scale: Literal["DeleageRoux", "Levitt", "ChouFasman"] = "DeleageRoux",
        precision: Optional[int] = None,
    ) -> Dict[str, float]:
        if scale not in secondary_structure_scales_by_name:
            raise ValueError(
                f"Scale '{scale}' not found in available secondary structure scales: {secondary_structure_scales_by_name.keys()}"
            )

        d = {}
        for structure_scale_name, structure_scale in secondary_structure_scales_by_name[
            scale
        ].items():
            d[structure_scale_name] = self.calc_property(
                scale=structure_scale,
                missing_aa_handling="error",
                aggregation_method="avg",
                normalize=True,
                weighting_scheme="uniform",
                precision=precision,
            )

        # normalize the values to sum to 1
        total = sum(d.values())
        if total > 0:
            for key in d:
                d[key] /= total

        if "coil" not in d:
            d["coil"] = None

        return d

    def aromaticity(
        self,
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
        """
        return self.aa_property_percentage(
            residues=list(aromatic_residues), precision=precision
        )
