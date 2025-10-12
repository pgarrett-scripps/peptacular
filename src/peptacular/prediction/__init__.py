"""
Streamlit UI components for peptacular.
Only available if streamlit is installed.
"""

try:
    import ms2pip

    MS2PIP_AVAILABLE = True
    from .psm_utils_core import predict_msms_spectra

except ImportError:
    MS2PIP_AVAILABLE = False

    def ms2pip_not_available(*args, **kwargs):  # type: ignore
        raise ImportError(
            "Ms2pip is not installed. Install with: pip install peptacular[Ms2pip]"
        )

    # Create dummy functions for all expected exports
    predict_msms_spectra = ms2pip_not_available

__all__ = ["MS2PIP_AVAILABLE", "predict_msms_spectra"]
