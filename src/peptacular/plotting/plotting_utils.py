from collections.abc import Sequence
from enum import StrEnum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd
    import plotly.graph_objects as go


class IonTypeColor(StrEnum):
    """Colors for different ion types in spectral plots."""

    A = "#3BC936"
    B = "#2E30C5"
    C = "#45B7D1"
    X = "#96CEB4"
    Y = "#DD4A30"
    Z = "#FF9FF3"
    I = "#9E0C8B"
    UNMATCHED = "#808080"

    @classmethod
    def get(cls, ion_type: str, default: str | None = None) -> str:
        """Get color for ion type (case-insensitive)."""
        try:
            return cls[ion_type.upper()].value
        except KeyError:
            return default or cls.UNMATCHED.value


class ChargeSymbol(StrEnum):
    """Marker symbols for different charge states."""

    CHARGE_1 = "circle"
    CHARGE_2 = "square"
    CHARGE_3 = "diamond"
    CHARGE_4 = "star"
    CHARGE_5 = "cross"

    @classmethod
    def get_by_index(cls, index: int) -> str:
        """Get symbol by index, returning the last symbol if index exceeds available symbols."""
        symbols = list(cls)
        return symbols[min(index, len(symbols) - 1)].value


def _get_charge_symbol_mapping(
    charges: Sequence[int],
) -> dict[int, str]:
    """Map fragment charges to marker symbols."""
    return {
        charge: ChargeSymbol.get_by_index(i) for i, charge in enumerate(sorted(charges))
    }


def _create_base_dataframe(
    mz_spectra: Sequence[float],
    intensity_spectra: Sequence[float],
    normalize: bool,
    baseline_offset: float,
) -> pd.DataFrame:
    """Create base spectrum DataFrame with standard columns."""
    import pandas as pd

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
            "color": IonTypeColor.UNMATCHED,
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


def _calculate_marker_sizes(
    intensities: Sequence[float], min_size: float, max_size: float
) -> list[float]:
    """Calculate marker sizes scaled to intensity values."""
    import numpy as np

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
