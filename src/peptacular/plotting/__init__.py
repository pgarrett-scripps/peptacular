try:
    import plotly
    import pandas as pd

    PLOTTING_AVAILABLE = True

    from .plotting_core import (
        plot_annotated_spectra_from_dataframe,
        apply_spectral_filters,
        prepare_spectrum_dataframe_from_matches,
        protein_coverage_plot,
    )

except ImportError:
    PLOTTING_AVAILABLE = False

    def plotly_not_available(*args, **kwargs):
        raise ImportError(
            "Plotly is not installed. Install with: pip install peptacular[plotly]"
        )

    # Create dummy functions for all expected exports
    plot_annotated_spectra_from_dataframe = plotly_not_available
    apply_spectral_filters = plotly_not_available
    prepare_spectrum_dataframe_from_matches = plotly_not_available
    protein_coverage_plot = plotly_not_available

__all__ = [
    "PLOTTING_AVAILABLE",
    "plot_annotated_spectra_from_dataframe",
    "apply_spectral_filters",
    "prepare_spectrum_dataframe_from_matches",
    "protein_coverage_plot",
]
