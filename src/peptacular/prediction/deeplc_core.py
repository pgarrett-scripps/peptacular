from collections.abc import Sequence
from typing import Any

from ..proforma.annotation import ProFormaAnnotation
from ..sequence import to_ms2_pip



def predict_rt_deeplc(
    peptides: Sequence[str | ProFormaAnnotation],
    calibration: tuple[Sequence[str | ProFormaAnnotation], Sequence[float]] | None = None,
    model: Any | None = None
) -> tuple[Sequence[float], object]:
    from deeplc import DeepLC
    import pandas as pd

    def find_bad_seqs(df: pd.DataFrame) -> pd.DataFrame:
        bad_aa_chars = {'X', 'B', 'J', 'U', 'Z', 'O'}
        def has_bad_aa(peptide: str) -> bool:
            return any(aa in bad_aa_chars for aa in peptide)
        df['is_bad'] = df['seq'].apply(has_bad_aa)        
        return df

    # Convert to MS2PIP format with original indices
    ms2_pip_results: list[tuple[str, str]] = to_ms2_pip(peptides)
    df = pd.DataFrame(ms2_pip_results, columns=['seq', 'modifications'])
    df['original_index'] = range(len(df))
    
    # Mark bad sequences
    df = find_bad_seqs(df)
    
    # Separate valid and invalid
    valid_df = df[~df['is_bad']].copy()
    
    # If all sequences are bad, return all -1.0
    if len(valid_df) == 0:
        return [-1.0] * len(peptides), model if model is not None else DeepLC()
    
    # Get unique (peptide, modifications) tuples for prediction
    # Keep track of first occurrence index for each unique tuple
    valid_df['tuple_key'] = valid_df['seq'] + '|' + valid_df['modifications']
    unique_df = valid_df.drop_duplicates(subset='tuple_key', keep='first').copy()
    unique_df = unique_df.drop(columns=['tuple_key', 'original_index', 'is_bad']).reset_index(drop=True)
    
    # Initialize or use existing model
    dlc = model if model is not None else DeepLC(batch_num_tf=124)
    
    # Calibrate if calibration data provided
    if calibration is not None and len(calibration[0]) > 0:
        calib_sequences, calib_rts = calibration
        calib_ms2_pip_results: list[tuple[str, str]] = to_ms2_pip(calib_sequences)
        calib_df = pd.DataFrame(calib_ms2_pip_results, columns=['seq', 'modifications'])
        calib_df['tr'] = list(calib_rts)
        
        # Only calibrate with valid sequences
        calib_df = find_bad_seqs(calib_df)
        calib_df = calib_df[~calib_df['is_bad']].drop(columns=['is_bad']).reset_index(drop=True)
        
        if len(calib_df) > 0:
            dlc.calibrate_preds(seq_df=calib_df)
    
    # Make predictions on unique valid sequences
    preds = dlc.make_preds(seq_df=unique_df, calibrate=(calibration is not None and len(calibration[0]) > 0))
    
    # Create mapping from (peptide, modifications) tuple to prediction
    unique_df['pred_rt'] = preds
    pred_map = dict(zip(
        unique_df['seq'] + '|' + unique_df['modifications'],
        unique_df['pred_rt']
    ))
    
    # Map predictions back to original order
    ordered_preds: list[float] = []
    for _, row in df.iterrows():
        if row['is_bad']:
            ordered_preds.append(-1.0)
        else:
            tuple_key = row['seq'] + '|' + row['modifications']
            ordered_preds.append(pred_map[tuple_key])
    
    return ordered_preds, dlc