from abc import ABC, abstractmethod
from typing import Literal, Sequence
from .adapters import PeptideData


class Aligner(ABC):
    """Abstract base class for retention time aligners."""
    
    @abstractmethod
    def align(self, datasets: Sequence[PeptideData]) -> PeptideData:
        """
        Align multiple datasets and return aggregated result.
        
        Args:
            datasets: Sequence of PeptideData to align
            
        Returns:
            Single PeptideData with aligned and aggregated retention times
        """
        pass


class LinearAligner(Aligner):
    """Aligns retention times across multiple datasets using linear regression."""
    
    def __init__(
        self,
        aggregation: Literal['median', 'mean', 'min', 'max'] = 'median',
        min_shared_peptides: int = 50,
        verbose: bool = False,
    ):
        """
        Initialize linear aligner.
        
        Args:
            aggregation: How to aggregate RTs for peptides in multiple datasets
            min_shared_peptides: Minimum shared peptides required for alignment
            verbose: Whether to print alignment information
        """
        self.aggregation = aggregation
        self.min_shared_peptides = min_shared_peptides
        self.verbose = verbose
        
        if aggregation not in ['median', 'mean', 'min', 'max']:
            raise ValueError(f"aggregation must be 'median', 'mean', 'min', or 'max', got '{aggregation}'")
    
    def align(self, datasets: Sequence[PeptideData]) -> PeptideData:
        """
        Align multiple datasets to the first dataset and aggregate.
        
        Args:
            datasets: Sequence of PeptideData to align
            
        Returns:
            Single PeptideData with aligned and aggregated retention times
        """
        if not datasets:
            raise ValueError("datasets cannot be empty")
        
        if len(datasets) == 1:
            return datasets[0]
        
        # Use first dataset as reference
        reference = datasets[0]
        reference_dict = dict(zip(reference.sequences, reference.retention_times))
        
        # Store aligned RTs: {sequence: [rt1, rt2, ...]}
        aligned_rts: dict[str, list[float]] = {
            seq: [rt] for seq, rt in zip(reference.sequences, reference.retention_times)
        }
        
        # Align each subsequent dataset to reference
        for i, dataset in enumerate(datasets[1:], start=1):
            if self.verbose:
                print(f"\nAligning dataset {i} to reference...")
            
            aligned = self._align_to_reference(reference_dict, dataset, i)
            
            # Add aligned RTs to collection
            for seq, rt in zip(dataset.sequences, aligned):
                if seq not in aligned_rts:
                    aligned_rts[seq] = []
                aligned_rts[seq].append(rt)
        
        # Aggregate RTs for each peptide
        aggregated_sequences: list[str] = []
        aggregated_rts: list[float] = []

        for seq, rts in aligned_rts.items():
            aggregated_sequences.append(seq)
            aggregated_rts.append(self._aggregate_values(rts))
        
        if self.verbose:
            print(f"\nFinal dataset: {len(aggregated_sequences)} unique peptides")
        
        return PeptideData(
            sequences=aggregated_sequences,
            retention_times=aggregated_rts
        )
    
    def _aggregate_values(self, values: list[float]) -> float:
        """Aggregate multiple values using configured method."""
        import numpy as np

        if self.aggregation == 'median':
            return float(np.median(values))
        elif self.aggregation == 'mean':
            return float(np.mean(values))
        elif self.aggregation == 'min':
            return float(np.min(values))
        elif self.aggregation == 'max':
            return float(np.max(values))
        else:
            raise ValueError(f"Unknown aggregation: {self.aggregation}")
    
    def _align_to_reference(
        self,
        reference_dict: dict[str, float],
        dataset: PeptideData,
        dataset_idx: int
    ) -> list[float]:
        import numpy as np

        """
        Align a dataset to reference using shared peptides.
        
        Args:
            reference_dict: {sequence: rt} for reference dataset
            dataset: Dataset to align
            dataset_idx: Index of dataset (for logging)
            
        Returns:
            List of aligned retention times
        """
        # Find shared peptides
        shared_sequences = [
            seq for seq in dataset.sequences 
            if seq in reference_dict
        ]
        
        if len(shared_sequences) < self.min_shared_peptides:
            raise ValueError(
                f"Dataset {dataset_idx} has only {len(shared_sequences)} shared peptides "
                f"with reference (minimum required: {self.min_shared_peptides})"
            )
        
        # Get RTs for shared peptides
        dataset_dict = dict(zip(dataset.sequences, dataset.retention_times))
        reference_rts = np.array([reference_dict[seq] for seq in shared_sequences])
        dataset_rts = np.array([dataset_dict[seq] for seq in shared_sequences])
        
        # Fit linear regression: reference = slope * dataset + intercept
        slope, intercept = np.polyfit(dataset_rts, reference_rts, 1)
        
        if self.verbose:
            print(f"  Shared peptides: {len(shared_sequences)}")
            print(f"  Alignment: slope={slope:.4f}, intercept={intercept:.4f}")
        
        # Apply transformation to all RTs in dataset
        aligned_rts = [
            float(slope * rt + intercept)
            for rt in dataset.retention_times
        ]
        
        return aligned_rts
    

class StatAligner(Aligner):
    """Simple aligner that just concatenates datasets without transformation."""
    
    def __init__(self, aggregation: Literal['median', 'mean', 'min', 'max'] = 'median'):
        self.aggregation = aggregation
    
    def align(self, datasets: Sequence[PeptideData]) -> PeptideData:
        """Concatenate datasets and aggregate duplicates."""
        import numpy as np

        all_rts: dict[str, list[float]] = {}
        
        for dataset in datasets:
            for seq, rt in zip(dataset.sequences, dataset.retention_times):
                if seq not in all_rts:
                    all_rts[seq] = []
                all_rts[seq].append(rt)
        
        sequences = list(all_rts.keys())
        retention_times = [
            float(np.median(rts)) if self.aggregation == 'median' else float(np.mean(rts))
            for rts in all_rts.values()
        ]
        
        return PeptideData(sequences=sequences, retention_times=retention_times)

