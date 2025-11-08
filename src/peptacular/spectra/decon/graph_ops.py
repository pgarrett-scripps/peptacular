"""
Graph construction and subgraph separation utilities.
"""

from typing import Callable, Literal, TypeVar

from .dclass import Graph, GraphEdge, GraphNode, IsotopeGap, PeakLike

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

    # Precompute bounds and offsets
    min_charge, max_charge = charge_range
    min_isotope_offset = isotope_mass / max_charge
    max_isotope_offset = isotope_mass / min_charge

    # Precompute valid offsets for each charge (already good)
    valid_isotope_offsets: list[tuple[int, float]] = [
        (charge, isotope_mass / charge) for charge in range(min_charge, max_charge + 1)
    ]

    # Create indexed peaks sorted by m/z
    indexed_peaks = [(i, peak.mz, peak) for i, peak in enumerate(peaks)]
    indexed_peaks.sort(key=lambda x: x[1])  # Sort by mz (index 1)

    n = len(indexed_peaks)

    # Precompute tolerance if it's DA (constant)
    use_constant_tol = tolerance_type == "da"
    constant_tol = tolerance if use_constant_tol else 0.0

    # Main loop - key optimization: early break conditions
    for idx_i in range(n):
        orig_i, mz_i, peak_i = indexed_peaks[idx_i]

        # Precompute tolerance once per peak_i if PPM
        tol_i = (
            constant_tol
            if use_constant_tol
            else get_tolerance(mz_i, tolerance, tolerance_type)
        )

        # Precompute search bounds
        min_search_mz = mz_i + min_isotope_offset - tol_i
        max_search_mz = mz_i + max_isotope_offset + tol_i

        for idx_j in range(idx_i + 1, n):
            orig_j, mz_j, peak_j = indexed_peaks[idx_j]

            # Early exit: peaks are sorted, so if we're past max, we're done
            if mz_j > max_search_mz:
                break

            # Skip if below minimum
            if mz_j < min_search_mz:
                continue

            mz_diff = mz_j - mz_i  # No need for abs() since sorted

            # Check potential matches for each possible charge
            for charge, offset in valid_isotope_offsets:
                diff_from_offset = abs(mz_diff - offset)
                if diff_from_offset <= tol_i:
                    # Apply custom filter if provided
                    if edge_filter is not None:
                        if not edge_filter(peak_i, peak_j, charge, mz_diff):
                            continue

                    # Compute error once
                    ppm_error = (mz_diff - offset) / mz_i * 1e6

                    peak_edge = GraphEdge(
                        IsotopeGap(
                            offset=offset,
                            charge=charge,
                            mz_error=mz_diff - offset,
                            ppm_error=ppm_error,
                        )
                    )

                    graph.add_edge(orig_i, orig_j, peak_edge)

    # Assign edge indices
    for index, data in graph.edge_index_map().items():
        data[2].index = index

    return graph
