"""
Streamlit UI components for peptacular.
Only available if streamlit is installed.
"""

from typing import Any

try:
    import pandas as pd
    import streamlit as st

    STREAMLIT_AVAILABLE = True
    from .streamlit_core import annotated_spectrum_plot

except ImportError:
    STREAMLIT_AVAILABLE = False

    def streamlit_not_available(*args: Any, **kwargs: Any):
        raise ImportError(
            "Streamlit is not installed. Install with: "
            "pip install peptacular[streamlit]"
        )

    # Create dummy functions for all expected exports
    annotated_spectrum_plot = streamlit_not_available

__all__ = ["STREAMLIT_AVAILABLE", "annotated_spectrum_plot"]
