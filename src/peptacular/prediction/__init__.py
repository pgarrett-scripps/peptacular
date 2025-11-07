"""
Streamlit UI components for peptacular.
Only available if streamlit is installed.
"""

from typing import Sequence

from ..sequence.basic import serialize

from ..proforma.annotation import ProFormaAnnotation
from .psm_utils_core import predict_msms_spectra
from .deeplc_core import predict_rt_deeplc
from .booster import RtPredictor, PeptideData, train_test_split

def predict_rt(
    peptides: Sequence[str | ProFormaAnnotation],
    calibration: tuple[Sequence[str | ProFormaAnnotation], Sequence[float]] | None = None,
    clip: tuple[float, float] | None = None,
    model: str | RtPredictor = '/home/patrick-garrett/Repos/peptacular/default.xgb',
    verbose: bool = True
) -> tuple[Sequence[float], object]:
    
    rt_predictor: RtPredictor = RtPredictor()

    if isinstance(model, str):
        rt_predictor.load_model(model)
    else:
        rt_predictor = model

    if any(isinstance(p, ProFormaAnnotation) for p in peptides):
        serialized_peptides = serialize(peptides)
    else:
        serialized_peptides = peptides  # type: ignore

    calib_data: PeptideData | None = None
    if calibration is not None:
        calib_peptides, calib_rts = calibration
        if any(isinstance(p, ProFormaAnnotation) for p in calib_peptides):
            serialized_calib_peptides = serialize(calib_peptides)
        else:
            serialized_calib_peptides = calib_peptides  # type: ignore
            
        calib_data = PeptideData(
            sequences=serialized_calib_peptides,
            retention_times=calib_rts
        )
        rt_predictor.align_to_observed(calib_data)

    preds = rt_predictor.predict(serialized_peptides, apply_alignment=True, clip=clip)

    if verbose and calib_data is not None:
        print(rt_predictor.score_summary(calib_data))
        
    return preds, rt_predictor