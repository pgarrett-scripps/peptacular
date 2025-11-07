"""
Functions to traverse (navigate) left or right in the graph based on
peak intensity and m/z relationships.
"""

from .dclass import Graph, PeakLike


def navigate_left(
    graph: Graph[PeakLike],
    start_node_idx: int,
    charge: int,
    max_intensity_change: float,
) -> list[int]:
    """
    Navigates left in m/z space (toward smaller m/z), collecting connected peaks
    of the same charge while ensuring each next peak's intensity is not too large
    a drop from the current peak.
    """
    current_node = start_node_idx
    path = [current_node]

    while True:
        neighbors = graph.neighbors(current_node)
        valid_next: list[tuple[float, int]] = []
        current_mz = graph[current_node].value.mz

        for n in neighbors:
            edge_data_list = graph.get_edge_data(current_node, n)
            if edge_data_list is None:
                continue
            for edge_data in (
                edge_data_list if isinstance(edge_data_list, list) else [edge_data_list]
            ):
                if (
                    edge_data.value.charge == charge
                    and graph[n].value.intensity < graph[current_node].value.intensity
                    and graph[n].value.mz < current_mz
                    and not graph[n].seen
                    and (
                        (graph[current_node].value.intensity - graph[n].value.intensity)
                        / graph[current_node].value.intensity
                    )
                    < max_intensity_change
                ):
                    valid_next.append((graph[n].value.intensity, n))

        # If no valid next peaks, break
        if not valid_next:
            break

        # Select the neighbor with the largest intensity
        valid_next.sort(key=lambda x: x[0], reverse=True)
        next_node = valid_next[0][1]

        path.append(next_node)
        current_node = next_node

    return path


def navigate_right(
    graph: Graph[PeakLike],
    start_node_idx: int,
    charge: int,
    max_intensity_change: float,
) -> list[int]:
    """
    Navigates right in m/z space (toward larger m/z), collecting connected peaks
    of the same charge while ensuring each next peak's intensity is not too large
    a drop from the current peak.
    """
    current_node = start_node_idx
    path = [current_node]

    while True:
        neighbors = graph.neighbors(current_node)
        valid_next: list[tuple[float, int]] = []
        current_mz = graph[current_node].value.mz

        for n in neighbors:
            edge_data_list = graph.get_edge_data(current_node, n)
            if edge_data_list is None:
                continue
            for edge_data in (
                edge_data_list if isinstance(edge_data_list, list) else [edge_data_list]
            ):
                if (
                    edge_data.value.charge == charge
                    and graph[n].value.intensity < graph[current_node].value.intensity
                    and graph[n].value.mz > current_mz
                    and not graph[n].seen
                    and (
                        (graph[current_node].value.intensity - graph[n].value.intensity)
                        / graph[current_node].value.intensity
                    )
                    < max_intensity_change
                ):
                    valid_next.append((graph[n].value.intensity, n))

        # If no valid next peaks, break
        if not valid_next:
            break

        valid_next.sort(key=lambda x: x[0], reverse=True)
        next_node = valid_next[0][1]

        path.append(next_node)
        current_node = next_node

    return path
