from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Literal, NamedTuple, Sequence

from ...sequence.mod_builder import from_ms2_pip

if TYPE_CHECKING:
    import pandas as pd


class PeptideData(NamedTuple):
    """Standardized peptide data format."""
    sequences: Sequence[str]  # ProForma format
    retention_times: Sequence[float]
    
    def __len__(self) -> int:
        return len(self.sequences)
    
    def __post_init__(self):
        if len(self.sequences) != len(self.retention_times):
            raise ValueError("sequences and retention_times must have the same length")


class FormatAdapter(ABC):
    """Base class for data format adapters."""
    
    @abstractmethod
    def load(self, input_source: "str | Path | pd.DataFrame") -> PeptideData:
        """Load data and return standardized PeptideData."""
        pass

class SageAdapter(FormatAdapter):
    """Adapter for SAGE search engine output."""
    
    def __init__(
        self,
        qvalue_threshold: float = 0.01,
        qvalue_column: str = "peptide_q",
        rt_column: str = "rt",
        peptide_column: str = "peptide",
        filename_column: str = "filename",
        score_column: str = "hyperscore",
        aggregation: Literal['median', 'mean', 'min', 'max'] = "median",
    ):
        """
        Initialize SAGE adapter.
        
        Args:
            qvalue_threshold: Filter PSMs by q-value
            qvalue_column: Column name for q-value
            rt_column: Column name for retention time
            peptide_column: Column name for peptide sequence
            score_column: Column name for score (for selecting best hit)
            aggregation: How to aggregate RTs per peptide ('median', 'mean', 'min', 'max')
        """
        self.qvalue_threshold = qvalue_threshold
        self.qvalue_column = qvalue_column
        self.rt_column = rt_column
        self.peptide_column = peptide_column
        self.filename_column = filename_column
        self.score_column = score_column
        self.aggregation = aggregation
        
        if aggregation not in ['median', 'mean', 'min', 'max']:
            raise ValueError(f"aggregation must be 'median', 'mean', 'min', or 'max', got '{aggregation}'")
    
    def load(self, input_source: "str | Path | pd.DataFrame") -> PeptideData:
        """
        Load SAGE data.
        
        Args:
            input_source: Path to SAGE TSV file or DataFrame
            
        Returns:
            PeptideData with sequences and retention times
        """
        # Load data
        import pandas as pd
        if isinstance(input_source, (str, Path)):
            df = pd.read_csv(input_source, sep="\t")
        else:
            df = input_source.copy()
        
        # Filter by q-value
        df = df[df[self.qvalue_column] <= self.qvalue_threshold]
        
        # Keep best hit per peptide
        df = df.loc[df.groupby([self.peptide_column, self.filename_column])[self.score_column].idxmax()]
        
        # Aggregate RTs per peptide
        df = df.groupby(self.peptide_column, as_index=False).agg({
            self.rt_column: self.aggregation
        })
        
        return PeptideData(
            sequences=df[self.peptide_column].tolist(),
            retention_times=df[self.rt_column].tolist()
        )


class MS2PIPAdapter(FormatAdapter):
    """Adapter for MS2PIP format."""
    
    def __init__(self, static_mods: dict[str, float] | None = None):
        """
        Initialize MS2PIP adapter.
        
        Args:
            static_mods: Dict of static mods, e.g. {'C': 57.02146}
        """
        self.static_mods = static_mods or {'C': 57.02146}

    def load(self, input_source: "str | Path | pd.DataFrame", max_values: int | None = None) -> PeptideData:
        """
        Load MS2PIP data.
        
        Args:
            input_source: Path to CSV or DataFrame
            
        Returns:
            PeptideData with sequences and retention times
        """
        import pandas as pd

        # Load data
        if isinstance(input_source, (str, Path)):
            df = pd.read_csv(input_source)
        else:
            df = input_source.copy()

        if max_values is not None:
            df = df.head(max_values)

        # Clean
        df['modifications'] = df['modifications'].fillna("")
        df['seq'] = df['seq'].astype(str)
        df['modifications'] = df['modifications'].astype(str)
        df['tr'] = df['tr'].astype(float)
        
        # Convert to ProForma
        ms2pip_tuples = list(zip(df['seq'], df['modifications']))
        sequences = from_ms2_pip(ms2pip_tuples, static_mods=self.static_mods)
        
        return PeptideData(
            sequences=sequences,
            retention_times=df['tr'].tolist()
        )