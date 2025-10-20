"""
Graph construction and subgraph separation utilities.
"""

from typing import Literal, Callable, TypeVar, Generic
from .dclass import GraphNode, GraphEdge, Graph, IsotopeGap, PeakLike, SpectrumPeak


# Generic type variable for peak types
P = TypeVar("P", bound=PeakLike)


def get_tolerance(
    mz: float, tolerance: float, tolerance_type: Literal["ppm", "da"]
) -> float:
    """
    Returns the absolute tolerance in Da, given a base m/z and either
    a parts-per-million (ppm) or Dalton (da) setting.
    """
    if tolerance_type == "ppm":
        return mz * tolerance / 1e6
    elif tolerance_type == "da":
        return tolerance
    else:
        raise ValueError(f"Unknown tolerance_type: {tolerance_type}")


def construct_graph(
    peaks: list[P],
    tolerance: float,
    tolerance_type: Literal["ppm", "da"],
    charge_range: tuple[int, int],
    isotope_mass: float,
    edge_filter: Callable[[P, P, int, float], bool] | None = None,
) -> Graph[P]:
    """
    Constructs a Graph from a list of peak-like objects.
    """
    graph: Graph[P] = Graph()

    # Add all peaks as nodes
    graph.add_nodes_from([GraphNode(peak) for peak in peaks])

    # Assign node indices
    for index in graph.node_indices():
        graph[index].index = index

    # Determine the smallest and largest offset for the given charge range
    min_isotope_offset = isotope_mass / charge_range[1]
    max_isotope_offset = isotope_mass / charge_range[0]

    # Precompute valid offsets for each charge
    valid_isotope_offsets: list[tuple[int, float]] = []
    for charge in range(charge_range[0], charge_range[1] + 1):
        valid_isotope_offsets.append((charge, isotope_mass / charge))

    # FIXED: Create a sorted list with original indices
    # Each element is (original_index, peak)
    indexed_peaks = [(i, peak) for i, peak in enumerate(peaks)]
    indexed_peaks.sort(key=lambda x: x[1].mz)

    # Loop over peaks and connect them if they match an isotope offset
    for idx_i in range(len(indexed_peaks)):
        for idx_j in range(idx_i + 1, len(indexed_peaks)):
            orig_i, peak_i = indexed_peaks[idx_i]
            orig_j, peak_j = indexed_peaks[idx_j]

            mz_i = peak_i.mz
            mz_j = peak_j.mz
            mz_diff = abs(mz_j - mz_i)

            tol = get_tolerance(mz_i, tolerance, tolerance_type)

            # If difference is too small or too large, skip early
            if mz_diff < min_isotope_offset - tol:
                continue
            if mz_diff > max_isotope_offset + tol:
                break

            # Check potential matches for each possible charge
            for charge, offset in valid_isotope_offsets:
                if abs(mz_diff - offset) <= tol:
                    # Apply custom filter if provided
                    if edge_filter is not None:
                        if not edge_filter(peak_i, peak_j, charge, mz_diff):
                            continue

                    ppm_error = (mz_diff - offset) / mz_i * 1e6
                    peak_edge = GraphEdge(
                        IsotopeGap(
                            offset=offset,
                            charge=charge,
                            mz_error=mz_diff - offset,
                            ppm_error=ppm_error,
                        )
                    )
                    # FIXED: Use original indices, not sorted positions
                    graph.add_edge(orig_i, orig_j, peak_edge)

    # Assign edge indices for easy reference
    for index, data in graph.edge_index_map().items():
        data[2].index = index

    return graph
