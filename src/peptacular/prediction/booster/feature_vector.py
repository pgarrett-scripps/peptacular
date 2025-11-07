from collections.abc import Sequence
from dataclasses import dataclass
from typing import overload, Callable

import numpy as np

from ...proforma.dclasses import Mod
from ...util import parse_static_mods
from ...constants import AMINO_ACIDS, ParrallelMethod, ParrallelMethodLiteral
from ...mass_calc import mod_mass
from ...proforma.annotation import ProFormaAnnotation
from ...property.properties import (
    ChargeScale,
    HPLCScale,
    HydrophobicityScale,
    PhysicalPropertyScale,
    SurfaceAccessibilityScale,
    SecondaryStructureScale,
)
from ...property.core import calc_property
from ...property.types import MissingAAHandling, AggregationMethod
from ...sequence import get_annotation_input, parallel_apply_internal


# Constants
AA_TO_INDEX = {aa: idx for idx, aa in enumerate(sorted(list(AMINO_ACIDS)))}
NUM_AMINO_ACIDS = len(AA_TO_INDEX)
NUM_WINDOWS = 5
WINDOW_OVERLAP = 1
WINDOW_FEATURE_SIZE = NUM_WINDOWS + 1  # +1 for full peptide value


@dataclass
class FeatureSpec:
    """Specification for a single feature or feature group."""
    name: str
    size: int
    extractor: Callable[[ProFormaAnnotation, np.ndarray, int], None]


# Define all basic features (non-windowed)
BASIC_FEATURES = [
    FeatureSpec("aa_counts", NUM_AMINO_ACIDS, 
                lambda a, v, p: _extract_aa_counts(a, v, p)),
    FeatureSpec("first_aa", NUM_AMINO_ACIDS, 
                lambda a, v, p: _extract_first_aa(a, v, p)),
    FeatureSpec("last_aa", NUM_AMINO_ACIDS, 
                lambda a, v, p: _extract_last_aa(a, v, p)),
    FeatureSpec("aa_mods", NUM_AMINO_ACIDS, 
                lambda a, v, p: _extract_aa_mods(a, v, p)),
    FeatureSpec("nterm_mod", 1, 
                lambda a, v, p: _extract_nterm_mod(a, v, p)),
    FeatureSpec("cterm_mod", 1, 
                lambda a, v, p: _extract_cterm_mod(a, v, p)),
    FeatureSpec("length", 1, 
                lambda a, v, p: _extract_length(a, v, p)),
    FeatureSpec("mass", 1, 
                lambda a, v, p: _extract_mass(a, v, p)),
    FeatureSpec("charge", 1, 
                lambda a, v, p: _extract_charge(a, v, p)),
    FeatureSpec("nterm-pkj", 1, 
                lambda a, v, p: _extract_nterm_pkj(a, v, p)),
    FeatureSpec("cterm-pkj", 1, 
                lambda a, v, p: _extract_cterm_pkj(a, v, p)),
]

RT_BASE_FEATURES = BASIC_FEATURES


# Define all window features
@dataclass
class WindowFeatureSpec:
    """Specification for sliding window features."""
    name: str
    scale: any  # The scale object to use


WINDOW_FEATURES = [
    WindowFeatureSpec("hydrophobicity", HydrophobicityScale.PARKER),
    WindowFeatureSpec("charge", ChargeScale.PK_SIDECHAIN),
    WindowFeatureSpec("surface_access", SurfaceAccessibilityScale.JANIN),
    WindowFeatureSpec("flexibility", PhysicalPropertyScale.FLEXIBILITY_VIHINEN),
    WindowFeatureSpec("bulkiness", PhysicalPropertyScale.BULKINESS),
    WindowFeatureSpec("hplc", HPLCScale.BROWNE),
    WindowFeatureSpec("alpha_helix", SecondaryStructureScale.LEVITT_ALPHA_HELIX),
    WindowFeatureSpec("beta_sheet", SecondaryStructureScale.LEVITT_BETA_SHEET),
    WindowFeatureSpec("beta_turn", SecondaryStructureScale.LEVITT_BETA_TURN)
]

RT_WINDOW_FEATURES = [WINDOW_FEATURES[0]]  # Only hydrophobicity for RT model


# Extractor functions for basic features
def _extract_aa_counts(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    for aa in annot.stripped_sequence:
        vector[pos + AA_TO_INDEX[aa]] += 1


def _extract_first_aa(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    vector[pos + AA_TO_INDEX[annot.stripped_sequence[0]]] = 1.0


def _extract_last_aa(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    vector[pos + AA_TO_INDEX[annot.stripped_sequence[-1]]] = 1.0


def _extract_aa_mods(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    if annot.has_internal_mods:
        for idx, mods in annot.get_internal_mod_dict().items():
            aa = annot.stripped_sequence[idx]
            vector[pos + AA_TO_INDEX[aa]] += sum(mod_mass(m) for m in mods)

    if annot.has_static_mods:
        static_mod_dict: dict[str, list[Mod]] = parse_static_mods(annot.get_static_mod_list().data)
        static_mod_masses: dict[str, float] = {mod_str: mod_mass(mods) for mod_str, mods in static_mod_dict.items()}
        for aa in static_mod_masses:
            if aa in AA_TO_INDEX:
                vector[pos + AA_TO_INDEX[aa]] += static_mod_masses[aa]



def _extract_nterm_mod(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    if annot.has_nterm_mods:
        vector[pos] = sum(mod_mass(m) for m in annot.nterm_mods)


def _extract_cterm_mod(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    if annot.has_cterm_mods:
        vector[pos] = sum(mod_mass(m) for m in annot.cterm_mods)


def _extract_length(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    vector[pos] = len(annot.stripped_sequence)


def _extract_mass(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    vector[pos] = annot.mass(charge=0)


def _extract_charge(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    vector[pos] = annot.charge or 0

def _extract_nterm_pkj(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    vector[pos] = calc_property(
        sequence=annot.stripped_sequence[0],
        scale=ChargeScale.PK_NTERMINAL,
        missing_aa_handling=MissingAAHandling.AVG,
        aggregation_method=AggregationMethod.AVG,
        normalize=False,
    )

def _extract_cterm_pkj(annot: ProFormaAnnotation, vector: np.ndarray, pos: int):
    vector[pos] = calc_property(
        sequence=annot.stripped_sequence[-1],
        scale=ChargeScale.PK_CTERMINAL,
        missing_aa_handling=MissingAAHandling.AVG,
        aggregation_method=AggregationMethod.AVG,
        normalize=False,
    )


class PeptideFeatures:
    """Feature vector with named property access."""
    
    def __init__(self, feature_vector: np.ndarray):
        self._vector = feature_vector
        self._offsets = self._calculate_offsets()
        
        # Validate size
        expected_size = self._calculate_total_size()
        if self._vector.shape[0] != expected_size:
            raise ValueError(f"Expected {expected_size} features, got {self._vector.shape[0]}")
    
    @staticmethod
    def _calculate_offsets() -> dict[str, int]:
        """Calculate starting position for each feature."""
        offsets = {}
        pos = 0
        
        # Basic features
        for spec in BASIC_FEATURES:
            offsets[spec.name] = pos
            pos += spec.size
        
        # Window features
        for spec in WINDOW_FEATURES:
            offsets[spec.name] = pos
            pos += WINDOW_FEATURE_SIZE
        
        return offsets
    
    @staticmethod
    def _calculate_total_size() -> int:
        """Calculate total feature vector size."""
        basic_size = sum(spec.size for spec in BASIC_FEATURES)
        window_size = len(WINDOW_FEATURES) * WINDOW_FEATURE_SIZE
        return basic_size + window_size
    
    def _get_slice(self, name: str, size: int) -> np.ndarray:
        """Get feature slice by name."""
        start = self._offsets[name]
        return self._vector[start:start + size]
    
    def _get_scalar(self, name: str) -> float:
        """Get scalar feature by name."""
        return float(self._vector[self._offsets[name]])
    
    @property
    def vector(self) -> np.ndarray:
        return self._vector
    
    # Basic features
    @property
    def aa_counts(self) -> np.ndarray:
        return self._get_slice("aa_counts", NUM_AMINO_ACIDS)
    
    @property
    def first_aa_encoding(self) -> np.ndarray:
        return self._get_slice("first_aa", NUM_AMINO_ACIDS)
    
    @property
    def last_aa_encoding(self) -> np.ndarray:
        return self._get_slice("last_aa", NUM_AMINO_ACIDS)
    
    @property
    def aa_mod_sums(self) -> np.ndarray:
        return self._get_slice("aa_mods", NUM_AMINO_ACIDS)
    
    @property
    def nterm_mod_mass(self) -> float:
        return self._get_scalar("nterm_mod")
    
    @property
    def cterm_mod_mass(self) -> float:
        return self._get_scalar("cterm_mod")
    
    @property
    def peptide_length(self) -> float:
        return self._get_scalar("length")
    
    @property
    def peptide_mass(self) -> float:
        return self._get_scalar("mass")
    
    @property
    def charge_state(self) -> float:
        return self._get_scalar("charge")
    
    # Window features
    @property
    def hydrophobicity_windows(self) -> np.ndarray:
        return self._get_slice("hydrophobicity", WINDOW_FEATURE_SIZE)
    
    @property
    def charge_windows(self) -> np.ndarray:
        return self._get_slice("charge", WINDOW_FEATURE_SIZE)
    
    @property
    def surface_accessibility_windows(self) -> np.ndarray:
        return self._get_slice("surface_access", WINDOW_FEATURE_SIZE)
    
    @property
    def flexibility_windows(self) -> np.ndarray:
        return self._get_slice("flexibility", WINDOW_FEATURE_SIZE)
    
    @property
    def bulkiness_windows(self) -> np.ndarray:
        return self._get_slice("bulkiness", WINDOW_FEATURE_SIZE)
    
    # Helpers
    @property
    def first_aa(self) -> str:
        idx = int(np.argmax(self.first_aa_encoding))
        return next(aa for aa, i in AA_TO_INDEX.items() if i == idx)
    
    @property
    def last_aa(self) -> str:
        idx = int(np.argmax(self.last_aa_encoding))
        return next(aa for aa, i in AA_TO_INDEX.items() if i == idx)
    
    @property
    def aa_count_dict(self) -> dict[str, int]:
        return {aa: int(self.aa_counts[i]) for aa, i in AA_TO_INDEX.items() if self.aa_counts[i] > 0}
    
    @classmethod
    def get_feature_names(cls) -> list[str]:
        """Get all feature names."""
        names: list[str] = []
        aa_list = sorted(AA_TO_INDEX.keys())
        
        # Basic AA features (multi-dimensional)
        for spec in BASIC_FEATURES:
            if spec.size == NUM_AMINO_ACIDS:
                names.extend(f"{spec.name}_{aa}" for aa in aa_list)
            else:
                names.append(spec.name)
        
        # Window features
        for spec in WINDOW_FEATURES:
            names.append(f"{spec.name}_full")
            names.extend(f"{spec.name}_window_{i+1}" for i in range(NUM_WINDOWS))
        
        return names
    
    def __repr__(self) -> str:
        return (f"PeptideFeatures(length={self.peptide_length:.0f}, mass={self.peptide_mass:.2f}, "
                f"first_aa={self.first_aa}, last_aa={self.last_aa})")


def _single_generate_feature_vector(sequence: str | ProFormaAnnotation) -> PeptideFeatures:
    """Generate feature vector for a single sequence."""
    annot = get_annotation_input(sequence=sequence, copy=False)
    
    # Allocate vector
    total_size = PeptideFeatures._calculate_total_size()
    vector = np.zeros(total_size, dtype=np.float32)
    
    pos = 0
    
    # Extract basic features
    for spec in BASIC_FEATURES:
        spec.extractor(annot, vector, pos)
        pos += spec.size
    
    # Extract window features
    for spec in WINDOW_FEATURES:
        # First calculate full peptide value
        full_value = calc_property(
            sequence=annot.stripped_sequence,
            scale=spec.scale,
            missing_aa_handling=MissingAAHandling.AVG,
            aggregation_method=AggregationMethod.AVG,
            normalize=False,
        )
        vector[pos] = full_value
        pos += 1
        
        # Then calculate sliding windows
        windows = annot.generate_sliding_window_features(spec.scale, NUM_WINDOWS, WINDOW_OVERLAP)
        vector[pos:pos + NUM_WINDOWS] = windows
        pos += NUM_WINDOWS
    
    return PeptideFeatures(vector)


@overload
def generate_feature_vector(
    sequence: str | ProFormaAnnotation,
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> PeptideFeatures: ...


@overload
def generate_feature_vector(
    sequence: Sequence[str | ProFormaAnnotation],
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> list[PeptideFeatures]: ...


def generate_feature_vector(
    sequence: str | ProFormaAnnotation | Sequence[str | ProFormaAnnotation],
    n_workers: None = None,
    chunksize: None = None,
    method: ParrallelMethod | ParrallelMethodLiteral | None = None,
) -> PeptideFeatures | list[PeptideFeatures]:
    """Generate feature vector(s) for given sequence(s)."""
    if isinstance(sequence, Sequence) and not isinstance(sequence, (str, ProFormaAnnotation)):
        return parallel_apply_internal(
            _single_generate_feature_vector,
            sequence,
            n_workers=n_workers,
            chunksize=chunksize,
            method=method,
        )
    return _single_generate_feature_vector(sequence)