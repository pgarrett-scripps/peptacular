"""
ProForma Parser Module
"""

from typing import (
    Dict,
    Iterable,
    List,
    Optional,
    Union,
    Literal,
)

from .annot_digestion import ProFormaAnnotationDigestion
from .utils import get_aa_value
from ..weights import get_weights
from .property_data import (
    all_property_scales,
    pk_nterminal,
    pk_cterminal,
    pk_sidechain,
    secondary_structure_scales_by_name,
)
from ..utils2 import round_to_precision

AGGREGATION_METHODS = Literal["sum", "avg"]
MISSING_AA_HANDLING = Literal["zero", "avg", "min", "max", "median", "error", "skip"]
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


class ProFormaAnnotationProperty(ProFormaAnnotationDigestion):
    """
    A ProForma annotation with additional methods for serialization and parsing.
    This class is used to represent a ProForma annotation with a sequence, modifications, and other properties.
    """

    def calc_property(
        self,
        scale: Union[str, Dict[str, float]],
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = False,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
        min_weight: float = 0.1,
        max_weight: float = 1.0,
        precision: Optional[int] = None,
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

        values: List[float] = []
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

        return round_to_precision(result, precision)

    def hydrophobicity(
        self,
        scale: str = "Kyte-Doolittle",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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
        scale: str = "flexibility_vihinen",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def hydrophilicity(
        self,
        scale: str = "hydrophilicity_hop_wood",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def surface_accessibility(
        self,
        scale: str = "surface_accessibility_vergoten",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def polarity(
        self,
        scale: str = "polarity_grantham",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def mutability(
        self,
        scale: str = "mutability",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def codons(
        self,
        scale: str = "codons",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "sum",
        normalize: bool = False,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def bulkiness(
        self,
        scale: str = "bulkiness",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def recognition_factors(
        self,
        scale: str = "recognition_factors",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "sum",
        normalize: bool = False,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def transmembrane_tendency(
        self,
        scale: str = "transmembrane_tendency",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def average_buried_area(
        self,
        scale: str = "average_buried_area",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def hplc(
        self,
        scale: str = "hplc_meek_2_1",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def refractivity(
        self,
        scale: str = "refractivity",
        missing_aa_handling: MISSING_AA_HANDLING = "error",
        aggregation_method: AGGREGATION_METHODS = "avg",
        normalize: bool = True,
        weighting_scheme: WEIGHTING_SCHEMES = "uniform",
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

    def calc_window_property(
        self,
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

        weighting_scheme = get_weights(
            window_size,
            weights=weighting_scheme,
            min_weight=min_weight,
            max_weight=max_weight,
        )

        l: List[float] = []
        for window_sequence in self.sliding_windows(window_size):
            if len(window_sequence) != window_size:
                raise ValueError(
                    f"Window size {window_size} does not match sequence length {len(window_sequence)}."
                )

            # Calculate the property average for the current window
            window_value = window_sequence.calc_property(
                scale=scale,
                missing_aa_handling=missing_aa_handling,
                aggregation_method=aggregation_method,
                normalize=normalize,
                weighting_scheme=weighting_scheme,
                min_weight=min_weight,
                max_weight=max_weight,
                precision=precision,
            )
            l.append(window_value)

        return l

    def charge_at_ph(
        self,
        pH: float = 7.0,
        precision: Optional[int] = None,
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

        if not isinstance(nterm_pK, (int, float)):
            raise ValueError(f"Invalid pK value for N-terminal amino acid '{nterm}': {nterm_pK}")

        partial_charge: float = 1.0 / (10 ** (pH - nterm_pK) + 1.0)
        positive_charge += partial_charge

        # Side chain positive charges
        for aa in "KRH":
            count = float(aa_counts.get(aa, 0))
            if count > 0:
                pK = get_aa_value(
                    aa=aa, aa_data=pk_sidechain, missing_aa_handling="error"
                )

                if not isinstance(pK, (int, float)) or pK <= 0:
                    # Skip if pK is not a valid number or is non-positive
                    raise ValueError(
                        f"Invalid pK value for side chain amino acid '{aa}': {pK}"
                    )

                partial_charge: float = 1.0 / (10 ** (pH - pK) + 1.0)
                positive_charge += count * partial_charge

        # Calculate negative charge (acidic groups)
        negative_charge = 0.0

        # C-terminal charge
        cterm_pK = get_aa_value(
            aa=cterm, aa_data=pk_cterminal, missing_aa_handling="error"
        )

        if not isinstance(cterm_pK, (int, float)):
            raise ValueError(f"Invalid pK value for C-terminal amino acid '{cterm}': {cterm_pK}")
        
        partial_charge = 1.0 / (10 ** (cterm_pK - pH) + 1.0)
        negative_charge += partial_charge

        # Side chain negative charges
        for aa in "DECY":
            count = float(aa_counts.get(aa, 0))
            if count > 0:
                pK = get_aa_value(
                    aa=aa, aa_data=pk_sidechain, missing_aa_handling="error"
                )

                if not isinstance(pK, (int, float)) or pK <= 0:
                    # Skip if pK is not a valid number or is non-positive
                    raise ValueError(
                        f"Invalid pK value for side chain amino acid '{aa}': {pK}"
                    )

                if pK > 0:  # Only calculate if pK exists (non-zero)
                    partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
                    negative_charge += count * partial_charge

        net_charge = positive_charge - negative_charge

        return round_to_precision(net_charge, precision)

    def pi(self, precision: Optional[int] = None) -> float:
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

        return round_to_precision(isoelectric_point, precision)

    def aa_property_percentage(
        self,
        residues: Iterable[str],
        precision: Optional[int] = None,
    ) -> float:
        """
        Calculates the percentage of specified amino acids in the sequence.
        """

        residue_perc = self.percent_residues(precision=None)
        val = sum(residue_perc.get(aa, 0) / 100 for aa in residues)

        return round_to_precision(val, precision)

    def secondary_structure(
        self,
        scale: Literal["DeleageRoux", "Levitt", "ChouFasman"] = "DeleageRoux",
        precision: Optional[int] = None,
    ) -> Dict[str, Optional[float]]:
        if scale not in secondary_structure_scales_by_name:
            raise ValueError(
                f"Scale '{scale}' not found in available secondary structure scales: {secondary_structure_scales_by_name.keys()}"
            )

        d: Dict[str, float] = {}
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
            d["coil"] = None # type: ignore

        return d # type: ignore

    def aromaticity(
        self,
        aromatic_residues: str = "YWF",
        precision: Optional[int] = None,
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
