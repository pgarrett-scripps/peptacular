from collections.abc import Sequence
from typing import Any

from ..proforma.annotation import ProFormaAnnotation
from ..sequence.util import get_annotation_input

def predict_msms_spectra(
    peptides: Sequence[str | ProFormaAnnotation],
    model: str | None = "HCD",
    processes: int | None = None) -> Any:
    """
    Predict MS/MS spectra for a list of peptides.
    """

    import ms2pip  # type: ignore
    from psm_utils import PSM, PSMList  # type: ignore

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

    return predictions

        