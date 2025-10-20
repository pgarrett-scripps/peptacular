"""
Data structure classes and constants used throughout the package.
"""

from typing import NamedTuple, Protocol, TypeVar, Generic
from dataclasses import dataclass
from .score import score_isotopic_pattern, IsotopicPatternScore
from ...isotope import estimate_isotopic_distribution, IsotopeLookup
from ...constants import PROTON_MASS, C13_NEUTRON_MASS


class PeakLike(Protocol):
    """Protocol for any peak-like object with mz and intensity."""

    mz: float
    intensity: float


# TypeVar bound to PeakLike protocol
P = TypeVar("P", bound=PeakLike)


class IsotopeGap(NamedTuple):
    """
    Stores information about isotope spacing between two peaks.
    """

    offset: float
    charge: int
    mz_error: float
    ppm_error: float


@dataclass
class SpectrumPeak:
    """Peak with m/z, intensity, and optional charge."""

    mz: float
    intensity: float
    charge: int | None = None

    @property
    def neutral_mass(self) -> float | None:
        """Calculate neutral mass if charge is known."""
        if self.charge is None or self.charge == 0:
            return None
        return self.mz * abs(self.charge) - self.charge * PROTON_MASS


@dataclass
class DeconvolutedPeak(Generic[P]):
    """
    Represents the result of deconvolution on an isotopic envelope.
    Generic over peak type P.
    """

    peaks: list[P]
    charge: int | None
    score: IsotopicPatternScore | None = None

    @property
    def base_peak(self) -> P:
        """Returns the monoisotopic peak."""
        return self.peaks[0]

    @property
    def largest_peak(self) -> P:
        """Returns the largest peak."""
        return max(self.peaks, key=lambda x: x.intensity)

    @property
    def mz_window(self) -> tuple[float, float]:
        """Returns the m/z window of the isotopic envelope."""
        return self.peaks[0].mz, self.peaks[-1].mz

    @property
    def intensity_window(self) -> tuple[float, float]:
        """Returns the intensity window of the isotopic envelope."""
        return self.peaks[0].intensity, self.peaks[-1].intensity

    @property
    def total_intensity(self) -> float:
        """Returns the total intensity of the isotopic envelope."""
        return sum(p.intensity for p in self.peaks)

    @property
    def num_peaks(self) -> int:
        """Returns the number of peaks in the isotopic envelope."""
        return len(self.peaks)

    @property
    def base_peak_neutral_mass(self) -> float | None:
        """Returns the neutral (uncharged) mass if charge is known."""
        if self.charge is None or self.base_peak.mz is None:
            return None
        return self.base_peak.mz * self.charge - self.charge * PROTON_MASS

    @property
    def largest_peak_neutral_mass(self) -> float | None:
        """Returns the neutral (uncharged) mass if charge is known."""
        if self.charge is None or self.largest_peak.mz is None:
            return None
        return self.largest_peak.mz * self.charge - self.charge * PROTON_MASS

    @property
    def isotope_gaps_mz(self) -> list[float]:
        """Returns a list of m/z gaps between peaks."""
        gaps: list[float] = []
        for i in range(1, len(self.peaks)):
            gap = self.peaks[i].mz - self.peaks[i - 1].mz
            gaps.append(gap)
        return gaps

    @property
    def isotope_gaps_neutral_mass(self) -> list[float]:
        """Returns a list of neutral mass gaps between peaks."""
        isotope_gaps = self.isotope_gaps_mz
        gaps: list[float] = []
        for gap in isotope_gaps:
            gaps.append(gap * self.charge)
        return gaps

    @property
    def isotope_gaps_ppm_error(self) -> list[float]:
        """Returns a list of ppm errors associated with isotope gaps."""
        gap_errors: list[float] = []
        for i in range(1, len(self.peaks)):
            theo_mz = self.peaks[i - 1].mz + (C13_NEUTRON_MASS / self.charge)
            expt_mz = self.peaks[i].mz
            gap_errors.append((expt_mz - theo_mz) / theo_mz * 1e6)
        return gap_errors

    @property
    def combined_score(self) -> float | None:
        """Returns the combined score from the scoring metrics."""
        return self.score.combined_score if self.score is not None else None

    def calculate_score(self, isotope_lookup: IsotopeLookup | None = None) -> None:
        """Score the fit of the observed isotopic envelope to the theoretical distribution."""
        if self.charge is None or self.largest_peak_neutral_mass is None:
            return None

        if not self.peaks or len(self.peaks) <= 1:
            return None

        theoretical_dist = (
            isotope_lookup.get_isotope_pattern(self.largest_peak_neutral_mass)
            if isotope_lookup is not None
            else estimate_isotopic_distribution(
                self.largest_peak_neutral_mass,
                max_isotopes=20,
                output_masses_for_neutron_offset=False,
                min_abundance_threshold=0.005,
                use_neutron_count=True,
                is_abundance_sum=True,
            )
        )

        observed: list[float] = [p.intensity for p in self.peaks]
        theo: list[float] = [abundance for _, abundance in theoretical_dist]
        score = score_isotopic_pattern(
            observed, theo, min_intensity=0.0, offset_range=0
        )
        self.score = score

    def __repr__(self) -> str:
        return f"DeconvolutedPeak(peaks={len(self.peaks)}, charge={int(self.charge) if self.charge is not None else 'None'}, score={self.score.combined_score if self.score is not None else 'None'})"

    def __str__(self) -> str:
        return f"DeconvolutedPeak: charge={int(self.charge) if self.charge is not None else 'None'}, peaks={len(self.peaks)}, score={self.score.combined_score if self.score is not None else 'None'}"


@dataclass
class GraphNode(Generic[P]):
    """
    Node wrapper for graph nodes, tagging each node with
    the original peak data and an index.
    """

    value: P
    index: int | None = None
    seen: bool = False

    def __str__(self) -> str:
        return f"GraphNode: {self.value} @ index: {self.index}"


@dataclass
class GraphEdge:
    """
    Edge wrapper for graph edges, carrying isotope gap info.
    """

    value: IsotopeGap
    index: int | None = None

    def __str__(self) -> str:
        return f"GraphEdge: {self.value} @ index: {self.index}"


class Graph(Generic[P]):
    """
    Pure Python graph implementation, generic over peak type P.
    P must satisfy the PeakLike protocol (have mz and intensity attributes).
    """

    def __init__(self):
        self._nodes: dict[int, GraphNode[P]] = {}
        self._edges: dict[int, tuple[int, int, GraphEdge]] = {}
        self._adjacency: dict[int, dict[int, list[int]]] = {}
        self._next_node_idx: int = 0
        self._next_edge_idx: int = 0

    def add_node(self, data: GraphNode[P]) -> int:
        """Add a single node and return its index."""
        idx = self._next_node_idx
        data.index = idx
        self._nodes[idx] = data
        self._adjacency[idx] = {}
        self._next_node_idx += 1
        return idx

    def add_nodes_from(self, node_data_list: list[GraphNode[P]]) -> list[int]:
        """Add multiple nodes from a list and return their indices."""
        indices: list[int] = []
        for data in node_data_list:
            indices.append(self.add_node(data))
        return indices

    def add_edge(self, source_idx: int, target_idx: int, edge_data: GraphEdge) -> int:
        """Add an edge between two nodes."""
        edge_idx = self._next_edge_idx
        edge_data.index = edge_idx
        self._edges[edge_idx] = (source_idx, target_idx, edge_data)

        # Add to adjacency lists
        if target_idx not in self._adjacency[source_idx]:
            self._adjacency[source_idx][target_idx] = []
        self._adjacency[source_idx][target_idx].append(edge_idx)

        if source_idx not in self._adjacency[target_idx]:
            self._adjacency[target_idx][source_idx] = []
        self._adjacency[target_idx][source_idx].append(edge_idx)

        self._next_edge_idx += 1
        return edge_idx

    def neighbors(self, node_idx: int) -> list[int]:
        """Return list of neighbor node indices."""
        return list(self._adjacency[node_idx].keys())

    def get_edge_data(
        self, source_idx: int, target_idx: int
    ) -> GraphEdge | list[GraphEdge] | None:
        """Get edge data between two nodes."""
        if target_idx not in self._adjacency[source_idx]:
            return None

        edge_indices = self._adjacency[source_idx][target_idx]
        if len(edge_indices) == 1:
            return self._edges[edge_indices[0]][2]
        else:
            return [self._edges[edge_idx][2] for edge_idx in edge_indices]

    def node_indices(self) -> list[int]:
        """Return all node indices."""
        return list(self._nodes.keys())

    def edge_index_map(self) -> dict[int, tuple[int, int, GraphEdge]]:
        """Return mapping of edge indices to (source, target, edge_data)."""
        return self._edges.copy()

    def __getitem__(self, node_idx: int) -> GraphNode[P]:
        """Get node data by index."""
        return self._nodes[node_idx]

    def __setitem__(self, node_idx: int, data: GraphNode[P]) -> None:
        """Set node data by index."""
        self._nodes[node_idx] = data

    def __repr__(self) -> str:
        nodes_str = {idx: str(node) for idx, node in self._nodes.items()}
        edges_str = {
            idx: (src, tgt, str(edge)) for idx, (src, tgt, edge) in self._edges.items()
        }
        return f"Graph(\n Nodes: {nodes_str},\n Edges: {edges_str}\n)"
