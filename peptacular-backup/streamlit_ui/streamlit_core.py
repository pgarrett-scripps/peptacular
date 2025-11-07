from typing import Sequence

import streamlit as st

from ..constants import IonType
from ..fragment.types import Fragment
from ..plotting import (
    apply_spectral_filters,
    plot_annotated_spectra_from_dataframe,
    prepare_spectrum_dataframe_from_matches,
)
from ..proforma.annotation import ProFormaAnnotation
from ..score import MatchMode, Scorer, ToleranceType
from ..sequence.util import get_annotation_input


def annotated_spectrum_plot(
    peptide: str | ProFormaAnnotation,
    mz_spectra: Sequence[float],
    intensity_spectra: Sequence[float],
    min_charge: int = 1,
    max_charge: int = 5,
    default_charge: tuple[int, int] = (1, 2),
):
    input_container = st.container()
    chart_container = st.container()

    with input_container:
        c1, c2 = st.columns([3, 2])
        with c1:
            fragment_types = st.segmented_control(
                "Fragment Ions",
                selection_mode="multi",
                options=[
                    IonType.A,
                    IonType.B,
                    IonType.C,
                    IonType.X,
                    IonType.Y,
                    IonType.Z,
                ],
                default=[IonType.Y, IonType.B],
                help="Select fragment ion types to display",
                key="fragment_types",
                width="stretch",
            )
        with c2:
            neutral_losses = st.segmented_control(
                "Neutral Losses",
                options=["NH3", "H2O", "CO2"],
                selection_mode="multi",
                help="Select neutral losses to consider",
                key="neutral_losses",
                width="stretch",
            )

        c1, c2 = st.columns(2)
        with c1:
            charge_range = st.slider(
                "Charge Range",
                min_value=min_charge,
                max_value=max_charge,
                value=default_charge,
                help="Select the charge range for the peptide",
                key="charge_range",
            )
        with c2:
            isotope_range = st.slider(
                "Isotope Range",
                min_value=0,
                max_value=5,
                value=(0, 2),
                help="Select the isotope range for the peptide",
                key="isotope_range",
            )

        c1, c2, c3 = st.columns([1, 1, 1])
        tolerance_type: ToleranceType = c1.selectbox(
            "Tolerance Type", options=[ToleranceType.PPM, ToleranceType.TH], index=0
        )
        if tolerance_type == ToleranceType.PPM:
            tolerance_value = c2.number_input(
                "Tolerance Value", min_value=0, value=50, step=10
            )
        else:
            tolerance_value = c2.number_input(
                "Tolerance Value", min_value=0.0, value=0.1, step=0.1
            )

        match_mode = c3.selectbox(
            "Match Mode", options=[MatchMode.CLOSEST, MatchMode.LARGEST, MatchMode.ALL]
        )

        c1, c2 = st.columns(2)
        with c1:
            mono_isotopic = st.checkbox("Monoisotopic", value=True)
        with c2:
            immonium_ion = st.checkbox("Immonium Ion", value=False)

        annotation = get_annotation_input(peptide)

        charges = [c for c in range(charge_range[0], charge_range[1] + 1)]
        isotopes = [i for i in range(isotope_range[0], isotope_range[1] + 1)]
        if immonium_ion:
            fragment_types += [IonType.IMMONIUM]
        frags: list[Fragment] = list(
            annotation.fragment(
                ion_types=fragment_types,
                charges=charges,
                monoisotopic=mono_isotopic,
                water_loss="H2O" in neutral_losses,
                ammonia_loss="NH3" in neutral_losses,
                max_losses=1,
                isotopes=isotopes,
                return_type="fragment",
            )
        )

        # Apply spectral filtering
        mz_filtered, intensity_filtered = apply_spectral_filters(
            list(mz_spectra),
            list(intensity_spectra),
            intensity_threshold=(0.01, None),
            mz_range=(None, None),
        )

        # Create scorer to get fragment matches
        scorer = Scorer(
            experimental_spectra=(mz_filtered, intensity_filtered),
            fragments=frags,
            tolerance_type=tolerance_type,
            tolerance=tolerance_value,
            match_mode=match_mode,
            filter_fragments_with_iso_gap=True,
            filter_fragments_without_mono=True,
            remove_duplicate_matches=True,
            filter_losses_without_base=True,
        )

        # Generate DataFrame from fragment matches
        df = prepare_spectrum_dataframe_from_matches(
            mz_filtered, intensity_filtered, scorer.fragment_matches
        )

        # Get sequence string for display
        sequence_str = (
            annotation.sequence if hasattr(annotation, "sequence") else str(annotation)
        )

        # Plot using DataFrame-based function
        fig = plot_annotated_spectra_from_dataframe(
            df,
            peptide_sequence=sequence_str,
            show_error_plot=True,
            show_coverage_plot=True,
        )

    with chart_container:
        st.plotly_chart(fig, use_container_width=True)  # type: ignore
