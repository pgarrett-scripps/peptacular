# Enhanced spectrum plotting with DataFrame-based approach
import re
from collections.abc import Sequence
from typing import Any, Literal

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..fragment.types import Fragment
from ..proforma.annotation import ProFormaAnnotation
from ..score import FragmentMatch, MatchMode, Scorer, ToleranceType
from ..sequence.util import get_annotation_input

ION_TYPE_COLORS = {
    "a": "#3BC936",
    "b": "#2E30C5",
    "c": "#45B7D1",
    "x": "#96CEB4",
    "y": "#DD4A30",
    "z": "#FF9FF3",
    "i": "#9E0C8B",
    "unmatched": "#808080",
}

CHARGE_SYMBOLS = ["circle", "square", "diamond", "star", "cross"]


def _get_charge_symbol_mapping(
    fragment_matches: Sequence[FragmentMatch],
) -> dict[int, str]:
    """Map fragment charges to marker symbols."""
    charges = sorted(
        {match.fragment.charge for match in fragment_matches if match.fragment}
    )
    return {
        charge: CHARGE_SYMBOLS[min(i, len(CHARGE_SYMBOLS) - 1)]
        for i, charge in enumerate(charges)
    }


def _create_base_dataframe(
    mz_spectra: Sequence[float],
    intensity_spectra: Sequence[float],
    normalize: bool,
    baseline_offset: float,
) -> pd.DataFrame:
    """Create base spectrum DataFrame with standard columns."""
    df = pd.DataFrame(
        {
            "mz": mz_spectra,
            "intensity": intensity_spectra,
            "matched": False,
            "fragment_label": "",
            "ion_type": "",
            "charge": 0,
            "loss": 0.0,
            "error_ppm": 0.0,
            "theoretical_mz": 0.0,
            "isotope": 0,
            "color": ION_TYPE_COLORS["unmatched"],
            "marker_symbol": None,
            "line_dash": "solid",
            "show_marker": False,
            "show_label": False,
        }
    )

    if normalize and intensity_spectra:
        max_intensity = max(intensity_spectra)
        if max_intensity > 0:
            df["intensity"] = (df["intensity"] / max_intensity) * 100

    df["intensity"] += baseline_offset
    return df


def _update_matched_peaks(
    df: pd.DataFrame,
    fragment_matches: Sequence[FragmentMatch],
    charge_symbols: dict[int, str],
) -> None:
    """Update DataFrame with fragment match information."""
    for match in fragment_matches:
        if match.fragment is None:
            continue

        mask = (abs(df["mz"] - match.mz) < 1e-6) & (
            abs(df["intensity"] - match.intensity) < 1e-6
        )
        indices = df[mask].index

        if not len(indices):
            continue

        idx = indices[0]
        fragment = match.fragment
        ion_base = fragment.ion_type[0].lower()

        df.loc[
            idx,
            [
                "matched",
                "fragment_label",
                "ion_type",
                "charge",
                "loss",
                "error_ppm",
                "theoretical_mz",
                "isotope",
            ],
        ] = [
            True,
            fragment.label,
            fragment.ion_type,
            fragment.charge,
            fragment.loss,
            match.error_ppm,
            match.theo_mz,
            fragment.isotope,
        ]

        df.loc[idx, "color"] = ION_TYPE_COLORS.get(
            ion_base, ION_TYPE_COLORS["unmatched"]
        )
        df.loc[idx, "line_dash"] = "dash" if fragment.loss != 0 else "solid"
        df.loc[idx, "marker_symbol"] = charge_symbols.get(fragment.charge, "circle")
        df.loc[idx, ["show_marker", "show_label"]] = fragment.isotope == 0


def prepare_spectrum_dataframe(
    mz_spectra: Sequence[float],
    intensity_spectra: Sequence[float],
    scorer: Scorer,
    normalize: bool = False,
    baseline_offset: float = 0.0,
) -> pd.DataFrame:
    """Prepare DataFrame with spectrum data and fragment annotations."""
    df = _create_base_dataframe(
        mz_spectra, intensity_spectra, normalize, baseline_offset
    )
    charge_symbols = _get_charge_symbol_mapping(scorer.fragment_matches)
    _update_matched_peaks(df, scorer.fragment_matches, charge_symbols)
    return df


def _create_hover_template(row: pd.Series) -> str:
    """Create hover template based on peak match status."""
    if row["matched"]:
        return (
            f"Fragment: {row['fragment_label']}<br>"
            f"m/z: {row['theoretical_mz']:.4f} ({row['mz']:.4f})<br>"
            f"Error: {row['error_ppm']:.1f} ppm<br>"
            f"Intensity: {row['intensity']:.0f}<extra></extra>"
        )
    return f"m/z: {row['mz']:.4f}<br>Intensity: {row['intensity']:.0f}<extra></extra>"


def _add_spectrum_traces(
    fig: go.Figure,
    df: pd.DataFrame,
    line_width: float,
    label_size: int,
    font_family: str,
    baseline_offset: float,
    subplot_row: int,
    subplot_col: int,
) -> None:
    """Add spectrum lines, markers, and labels to figure."""
    add_trace_kwargs = {"row": subplot_row, "col": subplot_col} if subplot_row else {}

    # Calculate a fixed label offset based on the overall intensity range
    if not df.empty:
        intensity_range = df["intensity"].max() - baseline_offset
        label_offset = intensity_range * 0.01
    else:
        label_offset = 0

    for _, row in df.iterrows():
        hover_template = _create_hover_template(row)

        # Stem line
        fig.add_trace(
            go.Scatter(
                x=[row["mz"], row["mz"]],
                y=[baseline_offset, row["intensity"]],
                mode="lines",
                line=dict(width=line_width, color=row["color"], dash=row["line_dash"]),
                showlegend=False,
                hovertemplate=hover_template,
            ),
            **add_trace_kwargs,
        )

        # Marker
        if row["show_marker"]:
            # Make markers hollow for fragments with losses
            if row["loss"] != 0:
                marker_color = "white"
                marker_line_color = row["color"]
                marker_line_width = 1
            else:
                marker_color = row["color"]
                marker_line_color = row["color"]
                marker_line_width = 1

            fig.add_trace(
                go.Scatter(
                    x=[row["mz"]],
                    y=[row["intensity"]],
                    mode="markers",
                    marker=dict(
                        size=5,
                        color=marker_color,
                        symbol=row["marker_symbol"],
                        line=dict(width=marker_line_width, color=marker_line_color),
                    ),
                    showlegend=False,
                    hovertemplate=hover_template,
                ),
                **add_trace_kwargs,
            )

        # Label
        if row["show_label"]:
            annotation_kwargs = {
                "x": row["mz"],
                "y": row["intensity"]
                + label_offset,  # Use fixed offset instead of scaled offset
                "text": row["fragment_label"],
                "textangle": -90,
                "showarrow": False,
                "font": dict(size=label_size, color=row["color"], family=font_family),
                "xanchor": "center",
                "yanchor": "bottom",
            }
            if subplot_row:
                annotation_kwargs.update({"row": subplot_row, "col": subplot_col})

            fig.add_annotation(**annotation_kwargs)


def _get_charge_legend_labels(charge_symbols: dict[int, str]) -> dict[int, str]:
    """Generate legend labels for charge states, handling symbol collisions."""
    if not charge_symbols:
        return {}

    # Find which symbols are used by multiple charges
    symbol_to_charges = {}
    for charge, symbol in charge_symbols.items():
        if symbol not in symbol_to_charges:
            symbol_to_charges[symbol] = []
        symbol_to_charges[symbol].append(charge)

    # Generate labels
    labels = {}
    for charge, symbol in charge_symbols.items():
        charges_with_same_symbol = sorted(symbol_to_charges[symbol])

        if len(charges_with_same_symbol) > 1:
            # Multiple charges share this symbol
            min_charge = min(charges_with_same_symbol)
            if charge == min_charge:
                # This is the lowest charge with this symbol - add + to indicate "and higher"
                labels[charge] = f"{charge}+"
            else:
                # Higher charges with same symbol - don't show in legend
                labels[charge] = None
        else:
            # Only one charge uses this symbol - no + needed
            labels[charge] = str(charge)

    return labels


def _add_legends(fig: go.Figure, df: pd.DataFrame) -> None:
    """Add legend entries for ion types and charge states with proper grouping."""
    matched_df = df[df["matched"]]
    if matched_df.empty:
        return

    # Ion types - first entry gets the group title
    ion_types = sorted(set(matched_df["ion_type"].str[0].str.lower()))
    for i, ion_type in enumerate(ion_types):
        if ion_type in ION_TYPE_COLORS:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(size=10, color=ION_TYPE_COLORS[ion_type]),
                    name=ion_type,
                    showlegend=True,
                    legendgroup="ion_types",
                    legendgrouptitle_text="Ion" if i == 0 else None,
                )
            )

    # Charge states - get symbol mapping and legend labels
    charge_symbols = {
        row["charge"]: row["marker_symbol"]
        for _, row in matched_df.iterrows()
        if row["charge"] > 0
    }
    charge_labels = _get_charge_legend_labels(charge_symbols)

    legend_entries_added = 0
    for charge in sorted(charge_symbols.keys()):
        label = charge_labels.get(charge)
        if label is not None:  # Only add legend entry if label is not None
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker=dict(
                        size=10,
                        color="black",
                        symbol=charge_symbols[charge],
                        line=dict(width=1, color="white"),
                    ),
                    name=label,
                    showlegend=True,
                    legendgroup="charge_states",
                    legendgrouptitle_text=(
                        "Charge" if legend_entries_added == 0 else None
                    ),
                )
            )
            legend_entries_added += 1


def _get_legend_config_for_subplot(subplot_config: dict, has_matches: bool) -> dict:
    """Get legend configuration based on subplot layout and whether there are matches."""
    if not has_matches:
        return dict(showlegend=False)

    positions = subplot_config["positions"]

    # Position legend in upper left of spectrum plot
    if "spectrum" in positions:
        spectrum_pos = positions["spectrum"]
        if spectrum_pos == (1, 1):  # Single row with coverage
            return dict(
                yanchor="top",
                y=0.95,
                xanchor="left",
                x=0.02,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="lightgray",
                borderwidth=1,
            )
        elif spectrum_pos == (2, 1):  # Multi-row layout
            return dict(
                yanchor="top",
                y=0.65,
                xanchor="left",
                x=0.02,  # Adjusted for lower subplot
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="lightgray",
                borderwidth=1,
            )

    # Default positioning
    return dict(
        yanchor="top",
        y=0.95,
        xanchor="left",
        x=0.02,
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="lightgray",
        borderwidth=1,
    )


def _calculate_plot_ranges(
    df: pd.DataFrame, baseline_offset: float
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Calculate appropriate x and y ranges for the plot."""
    if df.empty:
        return (0, 1000), (0, 1000)

    x_margin = (df["mz"].max() - df["mz"].min()) * 0.05
    y_margin = (df["intensity"].max() - baseline_offset) * 0.1

    x_range = (df["mz"].min() - x_margin, df["mz"].max() + x_margin)
    y_range = (baseline_offset - y_margin, df["intensity"].max() + y_margin)

    return x_range, y_range


def _calculate_marker_sizes(
    intensities: Sequence[float], min_size: float, max_size: float
) -> list[float]:
    """Calculate marker sizes scaled to intensity values."""
    if len(intensities) == 0:
        return []

    intensities = np.array(intensities)
    if len(intensities) == 1:
        return [max_size]

    min_intensity = intensities.min()
    max_intensity = intensities.max()

    if max_intensity == min_intensity:
        return [max_size] * len(intensities)

    # Scale intensities to marker size range
    normalized = (intensities - min_intensity) / (max_intensity - min_intensity)
    marker_sizes = min_size + normalized * (max_size - min_size)

    return marker_sizes.tolist()


def _add_error_subplot(
    fig: go.Figure,
    df: pd.DataFrame,
    subplot_row: int,
    subplot_col: int,
    min_marker_size: float,
    max_marker_size: float,
) -> None:
    """Add mass error subplot for matched peaks."""
    matched_df = df[df["matched"]]
    if matched_df.empty:
        return

    # Calculate marker sizes based on intensity
    marker_sizes = _calculate_marker_sizes(
        matched_df["intensity"], min_marker_size, max_marker_size
    )

    # Apply same hollow marker logic for error plot
    marker_colors = []
    marker_line_colors = []
    marker_line_widths = []

    for _, row in matched_df.iterrows():
        if row["loss"] != 0:
            marker_colors.append("white")
            marker_line_colors.append(row["color"])
            marker_line_widths.append(1)
        else:
            marker_colors.append(row["color"])
            marker_line_colors.append(row["color"])
            marker_line_widths.append(1)

    fig.add_trace(
        go.Scatter(
            x=matched_df["mz"],
            y=matched_df["error_ppm"],
            mode="markers",
            marker=dict(
                size=marker_sizes,
                color=marker_colors,
                symbol=matched_df["marker_symbol"],
                line=dict(width=marker_line_widths, color=marker_line_colors),
            ),
            name="Mass Error",
            showlegend=False,
            hovertemplate="Fragment: %{text}<br>m/z: %{x:.4f}<br>Error: %{y:.1f} ppm<extra></extra>",
            text=matched_df["fragment_label"],
        ),
        row=subplot_row,
        col=subplot_col,
    )


def _add_coverage_subplot_vertical(
    fig: go.Figure,
    df: pd.DataFrame,
    peptide_sequence: str,
    font_family: str,
    axis_title_size: int,
    tick_size: int,
    show_grid: bool,
    subplot_row: int,
    subplot_col: int,
    min_marker_size: float,
    max_marker_size: float,
) -> None:
    """Add vertical fragment coverage subplot showing matched fragments by position and ion type."""
    matched_df = df[(df["matched"]) & (df["isotope"] == 0)]
    if matched_df.empty or not peptide_sequence:
        return

    present_ion_types = sorted(set(matched_df["ion_type"].str[0].str.lower()))
    all_ion_types = ["a", "b", "c", "x", "y", "z"]
    ion_type_order = [ion for ion in all_ion_types if ion in present_ion_types]

    if not ion_type_order:
        return

    ion_type_positions = {ion: i for i, ion in enumerate(ion_type_order)}

    # Extract fragment positions and ion types
    coverage_data = []
    for _, row in matched_df.iterrows():
        fragment_label = row["fragment_label"]
        ion_type = row["ion_type"][0].lower() if row["ion_type"] else ""

        if ion_type not in ion_type_positions:
            continue

        try:
            match = re.search(r"[a-z](\d+)", fragment_label.lower())
            if not match:
                continue
            position = int(match.group(1))

            if ion_type in ["a", "b", "c"]:
                aa_position = position
            else:
                aa_position = len(peptide_sequence) - position + 1

            if 1 <= aa_position <= len(peptide_sequence):
                coverage_data.append(
                    {
                        "aa_position": aa_position,
                        "ion_type": ion_type,
                        "ion_type_x": ion_type_positions[ion_type],
                        "color": row["color"],
                        "symbol": row["marker_symbol"],
                        "label": fragment_label,
                        "charge": row["charge"],
                        "error_ppm": row["error_ppm"],
                        "has_loss": row["loss"] != 0,
                        "intensity": row["intensity"],
                    }
                )
        except (ValueError, IndexError, AttributeError):
            continue

    if not coverage_data:
        return

    # Calculate marker sizes based on intensity
    intensities = [point["intensity"] for point in coverage_data]
    marker_sizes = _calculate_marker_sizes(
        intensities, min_marker_size, max_marker_size
    )

    # Apply jitter to overlapping points
    position_groups = {}
    for i, point in enumerate(coverage_data):
        point["marker_size"] = marker_sizes[i]
        key = (point["aa_position"], point["ion_type_x"])
        if key not in position_groups:
            position_groups[key] = []
        position_groups[key].append(point)

    jitter_amount = 0.15
    np.random.seed(42)

    jittered_data = []
    for key, points in position_groups.items():
        if len(points) == 1:
            jittered_data.extend(points)
        else:
            n_points = len(points)
            jitter_offsets = np.linspace(-jitter_amount, jitter_amount, n_points)

            for i, point in enumerate(points):
                jittered_point = point.copy()
                jittered_point["ion_type_x_jittered"] = (
                    point["ion_type_x"] + jitter_offsets[i]
                )
                jittered_data.append(jittered_point)

    # Add coverage points
    for point in jittered_data:
        x_pos = point.get("ion_type_x_jittered", point["ion_type_x"])

        # Apply hollow marker logic for coverage plot
        if point["has_loss"]:
            marker_color = "white"
            marker_line_color = point["color"]
            marker_line_width = 2
        else:
            marker_color = point["color"]
            marker_line_color = point["color"]
            marker_line_width = 2

        fig.add_trace(
            go.Scatter(
                x=[x_pos],
                y=[point["aa_position"]],
                mode="markers",
                marker=dict(
                    size=point["marker_size"],
                    color=marker_color,
                    symbol=point["symbol"],
                    line=dict(width=marker_line_width, color=marker_line_color),
                ),
                name=f"{point['ion_type']} ions",
                showlegend=False,
                hovertemplate=(
                    f"Fragment: {point['label']}<br>"
                    f"Position: {point['aa_position']}<br>"
                    f"Ion Type: {point['ion_type']}<br>"
                    f"Charge: {point['charge']}<br>"
                    f"Error: {point['error_ppm']:.1f} ppm<extra></extra>"
                ),
            ),
            row=subplot_row,
            col=subplot_col,
        )

    # Configure axes
    aa_positions = list(range(1, len(peptide_sequence) + 1))
    aa_labels = [aa for aa in peptide_sequence]

    fig.update_xaxes(
        tickfont=dict(size=tick_size, family=font_family),
        tickvals=list(range(len(ion_type_order))),
        ticktext=ion_type_order,
        range=(-0.5, len(ion_type_order) - 0.5),
        showgrid=show_grid,
        gridcolor="lightgray",
        row=subplot_row,
        col=subplot_col,
        zeroline=False,
    )
    fig.update_yaxes(
        tickfont=dict(size=tick_size - 2, family=font_family),
        tickvals=aa_positions,
        ticktext=aa_labels,
        range=(len(peptide_sequence) + 0.5, 0.5),
        showgrid=show_grid,
        gridcolor="lightgray",
        row=subplot_row,
        col=subplot_col,
        zeroline=False,
    )


def _get_subplot_config(show_error_plot: bool, has_coverage: bool):
    """Get subplot configuration based on enabled plots."""
    if not show_error_plot and not has_coverage:
        return None

    if has_coverage:
        if show_error_plot:
            return {
                "specs": [[{"rowspan": 1}, {"rowspan": 2}], [{"rowspan": 1}, None]],
                "subplot_titles": (None, "Coverage", None),
                "row_heights": [0.3, 0.7],
                "column_widths": [0.75, 0.25],
                "rows": 2,
                "cols": 2,
                "positions": {"error": (1, 1), "spectrum": (2, 1), "coverage": (1, 2)},
            }
        else:
            return {
                "specs": [[{"rowspan": 1}, {"rowspan": 1}]],
                "subplot_titles": (None, "Fragment Coverage"),
                "row_heights": None,
                "column_widths": [0.75, 0.25],
                "rows": 1,
                "cols": 2,
                "positions": {"spectrum": (1, 1), "coverage": (1, 2)},
            }
    else:
        return {
            "specs": [[{}], [{}]],
            "subplot_titles": (None, None),
            "row_heights": [0.25, 0.75],
            "column_widths": None,
            "rows": 2,
            "cols": 1,
            "positions": {"error": (1, 1), "spectrum": (2, 1)},
        }


def _update_axis_properties(
    fig: go.Figure,
    row: int,
    col: int,
    axis_type: str,
    title: str,
    title_size: int,
    tick_size: int,
    font_family: str,
    range_vals: tuple,
    show_grid: bool,
):
    """Update axis properties for subplot or main plot."""
    axis_dict = dict(
        title=dict(text=title, font=dict(size=title_size, family=font_family)),
        tickfont=dict(size=tick_size, family=font_family),
        range=range_vals,
        showgrid=show_grid,
        gridcolor="lightgray",
        zeroline=False,
    )

    if axis_type == "x":
        fig.update_xaxes(axis_dict, row=row, col=col)
    else:
        fig.update_yaxes(axis_dict, row=row, col=col)


def plot_spectrum_from_dataframe(
    df: pd.DataFrame,
    title: str | None = None,
    fig_width: int = 1000,
    fig_height: int = 600,
    line_width: float = 2.0,
    label_size: int = 10,
    baseline_offset: float = 0.0,
    font_family: str = "Arial",
    axis_title_size: int = 14,
    tick_size: int = 12,
    show_grid: bool = True,
    background_color: str = "white",
    normalize: bool = False,
    show_error_plot: bool = False,
    show_coverage_plot: bool = False,
    peptide_sequence: str = "",
    min_marker_size: float = 4.0,
    max_marker_size: float = 8.0,
) -> go.Figure:
    """Plot spectrum using prepared DataFrame."""
    has_coverage = show_coverage_plot and bool(peptide_sequence)
    has_matches = not df[df["matched"]].empty
    subplot_config = _get_subplot_config(show_error_plot, has_coverage)

    if subplot_config:
        fig = make_subplots(
            rows=subplot_config["rows"],
            cols=subplot_config["cols"],
            specs=subplot_config["specs"],
            row_heights=subplot_config["row_heights"],
            column_widths=subplot_config["column_widths"],
            vertical_spacing=0.05,
            horizontal_spacing=0.05,
            subplot_titles=subplot_config["subplot_titles"],
            shared_xaxes="columns",
        )

        positions = subplot_config["positions"]

        # Add subplots
        if show_error_plot:
            _add_error_subplot(
                fig,
                df,
                positions["error"][0],
                positions["error"][1],
                min_marker_size,
                max_marker_size,
            )

        if not df.empty:
            spectrum_pos = positions["spectrum"]
            _add_spectrum_traces(
                fig,
                df,
                line_width,
                label_size,
                font_family,
                baseline_offset,
                spectrum_pos[0],
                spectrum_pos[1],
            )

        if has_coverage:
            coverage_pos = positions["coverage"]
            _add_coverage_subplot_vertical(
                fig,
                df,
                peptide_sequence,
                font_family,
                axis_title_size,
                tick_size,
                show_grid,
                coverage_pos[0],
                coverage_pos[1],
                min_marker_size,
                max_marker_size,
            )

        _add_legends(fig, df)

        # Update layout and axes
        x_range, y_range = _calculate_plot_ranges(df, baseline_offset)
        matched_df = df[df["matched"]]

        if show_error_plot and not matched_df.empty:
            error_margin = abs(matched_df["error_ppm"]).max() * 0.1
            error_range = (
                matched_df["error_ppm"].min() - error_margin,
                matched_df["error_ppm"].max() + error_margin,
            )
            error_pos = positions["error"]
            _update_axis_properties(
                fig,
                error_pos[0],
                error_pos[1],
                "x",
                "",
                axis_title_size,
                tick_size,
                font_family,
                x_range,
                show_grid,
            )
            _update_axis_properties(
                fig,
                error_pos[0],
                error_pos[1],
                "y",
                "Error (ppm)",
                axis_title_size,
                tick_size,
                font_family,
                error_range,
                show_grid,
            )

        spectrum_pos = positions["spectrum"]
        _update_axis_properties(
            fig,
            spectrum_pos[0],
            spectrum_pos[1],
            "x",
            "m/z",
            axis_title_size,
            tick_size,
            font_family,
            x_range,
            show_grid,
        )
        _update_axis_properties(
            fig,
            spectrum_pos[0],
            spectrum_pos[1],
            "y",
            "Intensity" + (" (%)" if normalize else ""),
            axis_title_size,
            tick_size,
            font_family,
            y_range,
            show_grid,
        )

        # Get legend configuration based on subplot layout
        legend_config = _get_legend_config_for_subplot(subplot_config, has_matches)

        fig.update_layout(
            title=dict(
                text=title or "Mass Spectrum",
                font=dict(size=axis_title_size + 4, family=font_family),
            ),
            width=fig_width,
            height=fig_height,
            plot_bgcolor=background_color,
            font=dict(family=font_family),
            margin=dict(l=60, r=40, t=80, b=60),
            legend=legend_config,
        )

    else:
        # Single plot
        fig = go.Figure()

        if not df.empty:
            _add_spectrum_traces(
                fig,
                df,
                line_width,
                label_size,
                font_family,
                baseline_offset,
                None,
                None,
            )
            _add_legends(fig, df)

        x_range, y_range = _calculate_plot_ranges(df, baseline_offset)

        # Legend configuration for single plot
        legend_config = (
            dict(
                yanchor="top",
                y=0.95,
                xanchor="left",
                x=0.02,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="lightgray",
                borderwidth=1,
            )
            if has_matches
            else dict(showlegend=False)
        )

        fig.update_layout(
            title=dict(
                text=title or "Mass Spectrum",
                font=dict(size=axis_title_size + 4, family=font_family),
            ),
            xaxis=dict(
                title=dict(
                    text="m/z", font=dict(size=axis_title_size, family=font_family)
                ),
                tickfont=dict(size=tick_size, family=font_family),
                range=x_range,
                showgrid=show_grid,
                gridcolor="lightgray",
            ),
            yaxis=dict(
                title=dict(
                    text="Intensity" + (" (%)" if normalize else ""),
                    font=dict(size=axis_title_size, family=font_family),
                ),
                tickfont=dict(size=tick_size, family=font_family),
                range=y_range,
                showgrid=show_grid,
                gridcolor="lightgray",
                zeroline=False,
            ),
            width=fig_width,
            height=fig_height,
            plot_bgcolor=background_color,
            font=dict(family=font_family),
            margin=dict(l=60, r=40, t=60, b=60),
            legend=legend_config,
        )

    return fig


def _apply_intensity_filter(
    data_pairs: list[tuple[float, float]], threshold: float, is_max: bool
) -> list[tuple[float, float]]:
    """Apply single intensity threshold filter."""
    if threshold is None:
        return data_pairs

    if threshold < 0:
        raise ValueError(
            f"{'Maximum' if is_max else 'Minimum'} intensity threshold must be non-negative."
        )

    if threshold < 1 and data_pairs:
        max_intensity = max(intensity for _, intensity in data_pairs)
        threshold = threshold * max_intensity

    if is_max:
        return [
            (mz, intensity) for mz, intensity in data_pairs if intensity <= threshold
        ]
    return [(mz, intensity) for mz, intensity in data_pairs if intensity >= threshold]


def _apply_mz_filter(
    data_pairs: list[tuple[float, float]], bound: float, is_upper: bool
) -> list[tuple[float, float]]:
    """Apply single m/z range filter."""
    if bound is None:
        return data_pairs

    if is_upper:
        return [(mz, intensity) for mz, intensity in data_pairs if mz <= bound]
    return [(mz, intensity) for mz, intensity in data_pairs if mz >= bound]


def apply_spectral_filters(
    mz_list: Sequence[float],
    intensity_list: Sequence[float],
    intensity_threshold: tuple[float | None, float | None],
    mz_range: tuple[float | None, float | None],
) -> tuple[list[float], list[float]]:
    """Apply intensity and m/z filtering to spectral data."""
    data_pairs = list(zip(mz_list, intensity_list))

    intensity_min, intensity_max = intensity_threshold
    data_pairs = _apply_intensity_filter(data_pairs, intensity_min, False)
    data_pairs = _apply_intensity_filter(data_pairs, intensity_max, True)

    mz_lower, mz_upper = mz_range
    data_pairs = _apply_mz_filter(data_pairs, mz_lower, False)
    data_pairs = _apply_mz_filter(data_pairs, mz_upper, True)

    if not data_pairs:
        return [], []

    mz_filtered, intensity_filtered = zip(*data_pairs)
    return list(mz_filtered), list(intensity_filtered)


def prepare_spectrum_dataframe_from_matches(
    mz_spectra: Sequence[float],
    intensity_spectra: Sequence[float],
    fragment_matches: Sequence[FragmentMatch],
    normalize: bool = False,
    baseline_offset: float = 0.0,
) -> pd.DataFrame:
    """
    Prepare DataFrame with spectrum data and fragment match annotations.

    Parameters:
    -----------
    mz_spectra, intensity_spectra : Sequence[float]
        The spectrum m/z and intensity values.
    fragment_matches : Sequence[FragmentMatch]
        Pre-computed fragment matches to annotate.
    normalize : bool, default False
        Whether to normalize intensities to 100%.
    baseline_offset : float, default 0.0
        Offset to apply to all intensities.

    Returns:
    --------
    pd.DataFrame
        DataFrame with spectrum data and fragment annotations.
    """
    df = _create_base_dataframe(
        mz_spectra, intensity_spectra, normalize, baseline_offset
    )
    charge_symbols = _get_charge_symbol_mapping(fragment_matches)
    _update_matched_peaks(df, fragment_matches, charge_symbols)
    return df


def plot_annotated_spectra_from_dataframe(
    df: pd.DataFrame,
    peptide_sequence: str = "",
    title: str | None = None,
    show_error_plot: bool = True,
    show_coverage_plot: bool = True,
    min_marker_size: float = 4.0,
    max_marker_size: float = 8.0,
    **kwargs: dict[str, Any],
) -> go.Figure:
    """
    Plot an annotated mass spectrum from a prepared DataFrame.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame prepared by prepare_spectrum_dataframe_from_matches().
    peptide_sequence : str, default ""
        The peptide sequence for coverage plot.
    title : str | None, default None
        Plot title. If None, uses "Mass Spectrum".
    show_error_plot : bool, default True
        Whether to show mass error subplot.
    show_coverage_plot : bool, default True
        Whether to show fragment coverage subplot.
    min_marker_size : float, default 4.0
        Minimum marker size for intensity scaling.
    max_marker_size : float, default 6.0
        Maximum marker size for intensity scaling.
    **kwargs
        Additional plotting parameters (fig_width, fig_height, normalize, etc.).

    Returns:
    --------
    go.Figure
        The annotated plotly figure.
    """
    # Calculate statistics from DataFrame
    matched_df = df[df["matched"]]
    num_matches = len(matched_df)

    # Calculate matched intensity percentage
    if not df.empty and not matched_df.empty:
        total_intensity = df["intensity"].sum()
        matched_intensity = matched_df["intensity"].sum()
        matched_intensity_pct = (
            (matched_intensity / total_intensity) * 100 if total_intensity > 0 else 0
        )
    else:
        matched_intensity_pct = 0

    # Create enhanced title
    base_title = title or "Mass Spectrum"
    if peptide_sequence:
        full_title = f"{base_title}: {peptide_sequence}<br><sub>{num_matches} matches, {matched_intensity_pct:.1f}% intensity explained</sub>"
    else:
        full_title = f"{base_title}<br><sub>{num_matches} matches, {matched_intensity_pct:.1f}% intensity explained</sub>"

    return plot_spectrum_from_dataframe(
        df,
        title=full_title,
        show_error_plot=show_error_plot,
        show_coverage_plot=show_coverage_plot,
        peptide_sequence=peptide_sequence,
        min_marker_size=min_marker_size,
        max_marker_size=max_marker_size,
        **kwargs,
    )


def plot_annotated_spectra(
    peptide: str | ProFormaAnnotation,
    mz_spectra: Sequence[float],
    intensity_spectra: Sequence[float],
    fragments: Sequence[Fragment],
    tolerance_value: float = 50,
    tolerance_type: ToleranceType = ToleranceType.PPM,
    match_mode: MatchMode = MatchMode.CLOSEST,
    filter_fragments_with_iso_gap: bool = True,
    filter_fragments_without_mono: bool = True,
    remove_duplicate_matches: bool = True,
    filter_losses_without_base: bool = True,
    intensity_threshold: tuple[float | None, float | None] = (None, None),
    mz_range: tuple[float | None, float | None] = (None, None),
    show_error_plot: bool = True,
    show_coverage_plot: bool = True,
    min_marker_size: float = 4.0,
    max_marker_size: float = 8.0,
    **kwargs,
) -> go.Figure:
    """
    Plot an annotated mass spectrum with fragment ion assignments.

    Creates a spectrum plot with fragment annotations, including ion type coloring,
    charge state markers, and neutral loss indicators. Optionally includes mass
    error and fragment coverage subplots.

    Parameters:
    -----------
    peptide : str | ProFormaAnnotation
        The peptide sequence or ProForma annotation.
    mz_spectra, intensity_spectra : Sequence[float]
        The spectrum m/z and intensity values.
    fragments : Sequence[Fragment]
        Fragment ions to annotate.
    tolerance_value : float, default 50
        Tolerance for matching fragments to peaks.
    tolerance_type : ToleranceType, default PPM
        Tolerance type (PPM or TH).
    match_mode : MatchMode, default CLOSEST
        Matching strategy (CLOSEST, LARGEST, or ALL).
    filter_* : bool
        Various fragment filtering options.
    intensity_threshold : tuple[float | None, float | None], default (None, None)
        Min/max intensity thresholds. Values < 1 are treated as fractions.
    mz_range : tuple[float | None, float | None], default (None, None)
        Min/max m/z display range.
    show_error_plot : bool, default False
        Whether to show mass error subplot above the spectrum.
    show_coverage_plot : bool, default False
        Whether to show fragment coverage subplot below the spectrum.
    min_marker_size : float, default 4.0
        Minimum marker size for intensity scaling in error and coverage plots.
    max_marker_size : float, default 6.0
        Maximum marker size for intensity scaling in error and coverage plots.
    **kwargs
        Additional plotting parameters (title, normalize, etc.).

    Returns:
    --------
    go.Figure
        The annotated plotly figure.
    """
    annotation = get_annotation_input(peptide)

    # Apply spectral filtering
    mz_filtered, intensity_filtered = apply_spectral_filters(
        list(mz_spectra), list(intensity_spectra), intensity_threshold, mz_range
    )

    # Configure scorer to get fragment matches
    scorer = Scorer(
        experimental_spectra=(mz_filtered, intensity_filtered),
        fragments=fragments,
        tolerance_type=tolerance_type,
        tolerance=tolerance_value,
        match_mode=match_mode,
        filter_fragments_with_iso_gap=filter_fragments_with_iso_gap,
        filter_fragments_without_mono=filter_fragments_without_mono,
        remove_duplicate_matches=remove_duplicate_matches,
        filter_losses_without_base=filter_losses_without_base,
    )

    # Create DataFrame from fragment matches
    df = prepare_spectrum_dataframe_from_matches(
        mz_filtered,
        intensity_filtered,
        scorer.fragment_matches,
        kwargs.get("normalize", False),
        kwargs.get("baseline_offset", 0.0),
    )

    # Get sequence string for display
    sequence_str = (
        annotation.sequence if hasattr(annotation, "sequence") else str(annotation)
    )

    # Plot using the DataFrame-based function
    return plot_annotated_spectra_from_dataframe(
        df,
        peptide_sequence=sequence_str,
        show_error_plot=show_error_plot,
        show_coverage_plot=show_coverage_plot,
        min_marker_size=min_marker_size,
        max_marker_size=max_marker_size,
        **kwargs,
    )


def protein_coverage_plot(
    protein_sequence: str,
    psm_coverage_array: Sequence[float],
    peptide_coverage_array: Sequence[float],
    log_scale: bool = False,
) -> go.Figure:
    """
    Create side-by-side protein coverage visualizations for PSM and Peptide coverage.

    Args:
        protein_sequence: String of amino acid sequence
        psm_coverage_array: Array of integers representing PSM coverage
        peptide_coverage_array: Array of integers representing peptide coverage
        log_scale: If True, apply log10(x+1) transformation to coverage values

    Returns:
        Plotly figure with two side-by-side coverage heatmaps
    """

    # Validate inputs
    if len(protein_sequence) != len(psm_coverage_array) or len(protein_sequence) != len(
        peptide_coverage_array
    ):
        raise ValueError(
            "Protein sequence and both coverage arrays must be same length"
        )

    # Apply log scaling if requested
    if log_scale:
        psm_values = [np.log10(x + 1) if x >= 0 else -1 for x in psm_coverage_array]
        peptide_values = [
            np.log10(x + 1) if x >= 0 else -1 for x in peptide_coverage_array
        ]
        scale_label = "Log10(Count + 1)"
    else:
        psm_values = psm_coverage_array.copy()
        peptide_values = peptide_coverage_array.copy()
        scale_label = "Count"

    # Constants
    RESIDUES_PER_ROW = 30

    # Calculate number of rows needed
    num_residues = len(protein_sequence)
    num_rows = (num_residues + RESIDUES_PER_ROW - 1) // RESIDUES_PER_ROW

    # Get max values for consistent color scaling (exclude zeros for better scaling)
    max_psm_coverage = (
        max([x for x in psm_values if x > 0]) if any(x > 0 for x in psm_values) else 1
    )
    max_peptide_coverage = (
        max([x for x in peptide_values if x > 0])
        if any(x > 0 for x in peptide_values)
        else 1
    )
    global_max = max(max_psm_coverage, max_peptide_coverage)

    def prepare_heatmap_data(coverage_array, original_array, coverage_type):
        """Prepare data for a single heatmap with hover data"""
        heatmap_data = []
        hover_data = []

        for row in range(num_rows):
            start_idx = row * RESIDUES_PER_ROW
            end_idx = min(start_idx + RESIDUES_PER_ROW, num_residues)

            row_coverage = coverage_array[start_idx:end_idx].copy()
            row_original = original_array[start_idx:end_idx].copy()

            # Create hover text for this row
            hover_row = []
            for i in range(len(row_coverage)):
                seq_idx = start_idx + i
                if seq_idx < len(protein_sequence):
                    aa = protein_sequence[seq_idx]
                    pos = seq_idx + 1  # 1-indexed position
                    orig_cov = row_original[i]
                    hover_row.append(
                        f"Position: {pos}<br>Amino Acid: {aa}<br>{coverage_type}: {orig_cov}"
                    )
                else:
                    hover_row.append("")

            # Pad row to RESIDUES_PER_ROW if needed
            while len(row_coverage) < RESIDUES_PER_ROW:
                row_coverage.append(-1)  # Use -1 for empty cells
                hover_row.append("")

            heatmap_data.append(row_coverage)
            hover_data.append(hover_row)

        return heatmap_data, hover_data

    # Prepare data for both heatmaps
    psm_heat_data, psm_hover_data = prepare_heatmap_data(
        psm_values, psm_coverage_array, "PSM Coverage"
    )
    peptide_heat_data, peptide_hover_data = prepare_heatmap_data(
        peptide_values, peptide_coverage_array, "Peptide Coverage"
    )

    # Create subplots
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("PSM Coverage", "Peptide Coverage"),
        horizontal_spacing=0.1,
    )

    # Add PSM coverage heatmap with hover data
    fig.add_trace(
        go.Heatmap(
            z=psm_heat_data,
            customdata=psm_hover_data,
            hovertemplate="%{customdata}<extra></extra>",
            colorscale="Reds",
            colorbar=dict(
                title=dict(text=f"PSM {scale_label}", side="bottom"),
                orientation="h",
                x=0.225,
                y=-0.25,
                len=0.35,
                thickness=20,
            ),
            showscale=True,
            zmin=0,
            zmax=global_max,
            name="PSM",
        ),
        row=1,
        col=1,
    )

    # Add Peptide coverage heatmap with hover data
    fig.add_trace(
        go.Heatmap(
            z=peptide_heat_data,
            customdata=peptide_hover_data,
            hovertemplate="%{customdata}<extra></extra>",
            colorscale="Blues",
            colorbar=dict(
                title=dict(text=f"Peptide {scale_label}", side="bottom"),
                orientation="h",
                x=0.775,
                y=-0.25,
                len=0.35,
                thickness=20,
            ),
            showscale=True,
            zmin=0,
            zmax=global_max,
            name="Peptide",
        ),
        row=1,
        col=2,
    )

    # Add text annotations for amino acids
    for row in range(num_rows):
        for col_pos in range(RESIDUES_PER_ROW):
            seq_idx = row * RESIDUES_PER_ROW + col_pos

            if seq_idx < len(protein_sequence):
                aa = protein_sequence[seq_idx]

                # PSM side annotations
                psm_coverage = psm_coverage_array[seq_idx]
                psm_color = "#CCCCCC" if psm_coverage == 0 else "black"

                fig.add_annotation(
                    x=col_pos,
                    y=row,
                    text=aa,
                    showarrow=False,
                    font=dict(color=psm_color, size=14, family="monospace"),
                    xref="x",
                    yref="y",
                )

                # Peptide side annotations
                peptide_coverage = peptide_coverage_array[seq_idx]
                peptide_color = "#CCCCCC" if peptide_coverage == 0 else "black"

                fig.add_annotation(
                    x=col_pos,
                    y=row,
                    text=aa,
                    showarrow=False,
                    font=dict(color=peptide_color, size=14, family="monospace"),
                    xref="x2",
                    yref="y2",
                )

    # Calculate height based on number of rows
    plot_height = max(600, num_rows * 40 + 200)

    # Update layout
    fig.update_layout(
        title=dict(
            text=f"Protein Coverage Maps ({num_residues} residues){' - Log Scale' if log_scale else ''}",
            x=0.5,
            xanchor="center",
            font=dict(size=16),
        ),
        height=plot_height,
        width=1800,
        showlegend=False,
        margin=dict(l=20, r=20, t=80, b=180),
        font=dict(family="Arial", size=12),
    )

    # Remove all axis labels and ticks for both subplots
    for col in [1, 2]:
        fig.update_xaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            title_text="",
            row=1,
            col=col,
        )

        fig.update_yaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            title_text="",
            autorange="reversed",
            row=1,
            col=col,
        )

    return fig
