from collections.abc import Sequence
from typing import NamedTuple

import ms2pip  # type: ignore
import numpy as np
from psm_utils import PSM, PSMList  # type: ignore

from ..proforma.annotation import ProFormaAnnotation
from ..sequence.util import get_annotation_input


class PredictedSpectrum(NamedTuple):
    theoretical_mz: dict[str, dict[str, np.ndarray]]
    predicted_intensity: dict[str, dict[str, np.ndarray]]


def predict_msms_spectra(
    peptides: Sequence[str | ProFormaAnnotation],
    model: str | None = "HCD",
    processes: int | None = None,
) -> dict[str, PredictedSpectrum]:
    """
    Predict MS/MS spectra for a list of peptides.

    Parameters:
    -----------
    peptides : list[str]
        A list of peptide sequences.

    Returns:
    --------
    list[float]
        A list of predicted MS/MS spectra intensities.
    """
    if isinstance(peptides, str):
        annotations = [get_annotation_input(peptides, copy=False)]
    elif isinstance(peptides, Sequence):  # type: ignore
        annotations = [get_annotation_input(p, copy=False) for p in peptides]
    else:
        raise TypeError("Invalid input type. Expected str or Sequence.")

    psm_list: list[PSM] = []
    for i, annotations in enumerate(annotations):
        psm_list.append(PSM(peptidoform=annotations.serialize(), spectrum_id=i))

    print(f"Predicting spectra for {len(psm_list)} unique peptides...")
    predictions = ms2pip.predict_batch(
        psms=PSMList(psm_list=psm_list), model=model, processes=processes
    )

    psm_lookup: dict[str, PredictedSpectrum] = {}
    for res in predictions:
        if not res.theoretical_mz or not res.predicted_intensity:
            continue

        intensities: dict[str, dict[str, np.ndarray]] = {}
        for ion_type in res.predicted_intensity:
            intensities[ion_type] = (2 ** res.predicted_intensity[ion_type]) - 0.001  #  type: ignore

        sequence: str = res.psm.peptidoform.proforma  # type: ignore

        psm_lookup[sequence] = PredictedSpectrum(
            theoretical_mz=res.theoretical_mz,  # type: ignore
            predicted_intensity=intensities,
        )

    return psm_lookup
