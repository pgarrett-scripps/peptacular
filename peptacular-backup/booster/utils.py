from .adapters import PeptideData

from ...sequence import strip_mods


def train_test_split(
    data: PeptideData,
    test_size: float = 0.2,
    random_state: int | None = None,
    stratify_by_sequence: bool = True
) -> tuple[PeptideData, PeptideData]:
    """
    Split data into train and test sets.
    
    Args:
        data: PeptideData to split
        test_size: Fraction of data for test set (0.0 to 1.0)
        random_state: Random seed for reproducibility
        stratify_by_sequence: If True, ensures stripped sequences don't appear in both train and test
        
    Returns:
        Tuple of (train_data, test_data) as PeptideData objects
    """
    from sklearn.model_selection import train_test_split as sklearn_split
    import numpy as np
    
    if stratify_by_sequence:
        # Group by stripped sequence to avoid data leakage
        
        stripped = strip_mods(data.sequences)
        unique_stripped = list(set(stripped))
        
        # Split unique sequences
        train_stripped, test_stripped = sklearn_split(
            unique_stripped,
            test_size=test_size,
            random_state=random_state
        )
        
        train_stripped_set = set(train_stripped)
        test_stripped_set = set(test_stripped)
        
        # Assign sequences to train/test based on their stripped form
        train_indices = [i for i, s in enumerate(stripped) if s in train_stripped_set]
        test_indices = [i for i, s in enumerate(stripped) if s in test_stripped_set]
        
        train_data = PeptideData(
            sequences=[data.sequences[i] for i in train_indices],
            retention_times=[data.retention_times[i] for i in train_indices]
        )
        
        test_data = PeptideData(
            sequences=[data.sequences[i] for i in test_indices],
            retention_times=[data.retention_times[i] for i in test_indices]
        )
    else:
        # Simple random split
        indices = np.arange(len(data))
        train_idx, test_idx = sklearn_split(
            indices,
            test_size=test_size,
            random_state=random_state
        )
        
        train_data = PeptideData(
            sequences=[data.sequences[i] for i in train_idx],
            retention_times=[data.retention_times[i] for i in train_idx]
        )
        
        test_data = PeptideData(
            sequences=[data.sequences[i] for i in test_idx],
            retention_times=[data.retention_times[i] for i in test_idx]
        )
    
    return train_data, test_data