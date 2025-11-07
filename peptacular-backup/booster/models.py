import os
from enum import StrEnum
from pathlib import Path
from typing import Any, Sequence, overload

import numpy as np
from xgboost import XGBRegressor
import pandas as pd

from ...sequence.mod_builder import from_ms2_pip

from ...proforma import ProFormaAnnotation
from .feature_vector import generate_feature_vector
from ...sequence import strip_mods


def load_sage_files(sage_input: list[str | Path | pd.DataFrame] | str | Path | pd.DataFrame) -> pd.DataFrame:
    """
    Load and merge multiple SAGE result files.
    
    If a path is a directory, recursively searches for 'results.sage.tsv' files.
    If a path is a file, loads it directly.
    
    Args:
        file_paths: List of file paths or directory paths
        
    Returns:
        Merged DataFrame with all PSMs
    """


    if isinstance(sage_input, (str, Path, pd.DataFrame)):
        sage_input = [sage_input]
    
    
    dfs: list[pd.DataFrame] = []
    all_files: list[Path] = []
    
    # Expand directories to find results.sage.tsv files
    for file_path in sage_input:


        if isinstance(file_path, pd.DataFrame):
            dfs.append(file_path)
            continue
        elif isinstance(file_path, str) or isinstance(file_path, Path):
            file_path = Path(file_path)
        else:
            print(f"ERROR: Invalid input type: {type(file_path)}")
            continue

        if not file_path.exists():
            print(f"ERROR: Path not found: {file_path}")
            continue
        
        if file_path.is_dir():
            # Search recursively for results.sage.tsv files
            sage_files = list(file_path.rglob("results.sage.tsv"))
            if sage_files:
                print(f"Found {len(sage_files)} SAGE result files in {file_path}/")
                all_files.extend(sage_files)
            else:
                print(f"WARNING: No 'results.sage.tsv' files found in {file_path}/")
        else:
            # It's a file, add it directly
            all_files.append(file_path)
    
    print(f"\nLoading {len(all_files)} files...")
    
    # Load each file
    for file_path in all_files:
        try:
            print(f"  Loading {file_path.parent.name}/{file_path.name}...")
            df = pd.read_csv(file_path, sep="\t")
            
            # Add filename column if not present
            # Use parent directory name to distinguish files with same name
            if 'filename' not in df.columns:
                # Use parent directory name as identifier (more informative than "results.sage")
                df['filename'] = file_path.parent.name
            
            dfs.append(df)
            
        except Exception as e:
            print(f"  ERROR loading {file_path}: {e}")
            continue
    
    if not dfs:
        raise ValueError("No valid files loaded")
    
    merged_df = pd.concat(dfs, ignore_index=True)
    print(f"\nLoaded {len(merged_df)} total PSMs from {len(dfs)} files")
    
    return merged_df


class RtModelType(StrEnum):
    """Enum for predefined retention time models."""
    DEFAULT = "default"


# Mapping of model types to their file paths
RT_MODEL_PATHS = {
    RtModelType.DEFAULT: "models/rt_model_default.json",
}


class RtPredictor:
    """Retention time predictor using XGBoost."""
    
    def __init__(
        self,
        model: str | RtModelType | XGBRegressor = RtModelType.DEFAULT,
        verbose: bool = False,
    ):
        """
        Initialize predictor.
        
        Args:
            model: Either a RtModelType enum (loads from predefined path),
                   a string path to a model file, or a trained XGBRegressor instance
            verbose: Whether to print verbose output
        """
        self.verbose = verbose
        self.model = self._initialize_model(model)
        self.alignment_params: dict[str, float] | None = None  # For linear alignment
    
    def align_to_dataset(
        self,
        sequences: Sequence[str | ProFormaAnnotation],
        observed_rts: Sequence[float],
        use_quantile_filtering: bool = True,
        lower_quantile: float = 0.05,
        upper_quantile: float = 0.95,
    ) -> dict[str, float]:
        """
        Calculate linear alignment parameters to map predictions to observed RTs.
        
        This fits a linear function: aligned_rt = slope * predicted_rt + intercept
        
        Args:
            sequences: Peptide sequences to use for alignment
            observed_rts: Observed retention times for these sequences
            use_quantile_filtering: Whether to filter outliers before fitting
            lower_quantile: Lower quantile for filtering (default 0.05)
            upper_quantile: Upper quantile for filtering (default 0.95)
            
        Returns:
            Dictionary with alignment parameters: slope, intercept, min_rt, max_rt
        """
        # Get predictions
        predicted_rts = np.array(self.predict(sequences, apply_alignment=False))
        
        observed_rts = np.array(observed_rts)
        min_rt = np.min(observed_rts)
        max_rt = np.max(observed_rts)
        
        # Optional: filter outliers based on prediction error
        # DO THIS BEFORE NORMALIZATION
        if use_quantile_filtering:
            # Normalize for comparison
            observed_rts_normalized = (observed_rts - min_rt) / (max_rt - min_rt)
            errors = np.abs(predicted_rts - observed_rts_normalized)
            lower_bound = np.quantile(errors, lower_quantile)
            upper_bound = np.quantile(errors, upper_quantile)
            mask = (errors >= lower_bound) & (errors <= upper_bound)
            predicted_rts = predicted_rts[mask]
            observed_rts = observed_rts[mask]  # Keep original scale
            
            if self.verbose:
                print(f"Filtered to {len(predicted_rts)} points ({mask.sum()/len(mask)*100:.1f}%)")
        
        # Fit linear regression: observed = slope * predicted + intercept
        # Use ORIGINAL observed_rts (not normalized)
        slope, intercept = np.polyfit(predicted_rts, observed_rts, 1)
        
        # Store alignment parameters
        self.alignment_params = {
            'slope': float(slope),
            'intercept': float(intercept),
            'min_observed_rt': float(np.min(observed_rts)),
            'max_observed_rt': float(np.max(observed_rts)),
        }
        
        if self.verbose:
            print(f"Alignment parameters: slope={slope:.4f}, intercept={intercept:.4f}")
            print(f"RT range: [{self.alignment_params['min_observed_rt']:.2f}, "
                f"{self.alignment_params['max_observed_rt']:.2f}]")
        
        return self.alignment_params
    
    def clear_alignment(self) -> None:
        """Remove alignment parameters."""
        self.alignment_params = None
    
    def _apply_alignment(self, predictions: np.ndarray) -> np.ndarray:
        """Apply linear alignment to predictions."""
        if self.alignment_params is None:
            return predictions
        
        aligned = (predictions * self.alignment_params['slope'] + 
                   self.alignment_params['intercept'])
        
        return aligned
    
    @overload
    def predict(
        self, 
        sequences: str | ProFormaAnnotation, 
        gradient: float | None = None,
        apply_alignment: bool = True,
    ) -> float: ...
    
    @overload
    def predict(
        self, 
        sequences: Sequence[str | ProFormaAnnotation], 
        gradient: float | None = None,
        apply_alignment: bool = True,
    ) -> list[float]: ...
    
    def predict(
        self,
        sequences: Sequence[str | ProFormaAnnotation] | str | ProFormaAnnotation,
        gradient: float | None = None,
        apply_alignment: bool = True,
    ) -> list[float] | float:
        """
        Predict retention times for given sequences.
        
        Args:
            sequences: Peptide sequence(s) to predict
            gradient: Optional gradient multiplier (applied before alignment)
            apply_alignment: Whether to apply linear alignment (if available)
            
        Returns:
            Predicted retention time(s)
        """
        is_single = isinstance(sequences, (str, ProFormaAnnotation))
        sequences_list = [sequences] if is_single else list(sequences)
        
        X = self._sequences_to_features(sequences_list)
        predictions = self.model.predict(X)  # 0-1 range from training
        
        # Apply gradient if specified (legacy support)
        if gradient is not None:
            predictions = predictions * gradient
        
        # Apply alignment if available and requested
        if apply_alignment and self.alignment_params is not None:
            predictions = self._apply_alignment(predictions)
        
        res = float(predictions[0]) if is_single else predictions.tolist()
        return res
    
    def align_to_sage(
        self,
        sage_input: str | pd.DataFrame | list[str | Path | pd.DataFrame],
        qvalue_threshold: float = 0.01,
        qvalue_column: str = 'peptide_q',
        use_quantile_filtering: bool = True,
        lower_quantile: float = 0.05,
        upper_quantile: float = 0.95,
    ) -> dict[str, float]:
        """
        Align predictions to SAGE data.
        
        Args:
            sage_input: Path(s) to SAGE file(s) or DataFrame
            qvalue_threshold: Q-value threshold for filtering
            qvalue_column: Column name for q-value filtering
            use_quantile_filtering: Whether to filter outliers
            lower_quantile: Lower quantile for filtering
            upper_quantile: Upper quantile for filtering
            
        Returns:
            Alignment parameters
        """
        if isinstance(sage_input, (str, Path, pd.DataFrame)):
            df = load_sage_files([sage_input])
        else:
            df = load_sage_files(sage_input)
        
        # Filter by q-value
        df = df[df[qvalue_column] <= qvalue_threshold]
        
        # Get best hit per peptide per file
        df = df.loc[df.groupby(['peptide', 'filename'])['hyperscore'].idxmax()]
        
        # Get median RT per peptide
        df = df.groupby('peptide').agg({'rt': 'median'}).reset_index()
        
        sequences = df['peptide'].tolist()
        observed_rts = df['rt'].tolist()
        
        return self.align_to_dataset(
            sequences=sequences,
            observed_rts=observed_rts,
            use_quantile_filtering=use_quantile_filtering,
            lower_quantile=lower_quantile,
            upper_quantile=upper_quantile,
        )
    
    def align_to_ms2_pip(
        self,
        ms2pip_input: pd.DataFrame | str,
        use_quantile_filtering: bool = True,
        lower_quantile: float = 0.05,
        upper_quantile: float = 0.95,
    ) -> dict[str, float]:
        """
        Align predictions to MS2PIP data.
        
        Args:
            ms2pip_input: Path to MS2PIP CSV or DataFrame
            use_quantile_filtering: Whether to filter outliers
            lower_quantile: Lower quantile for filtering
            upper_quantile: Upper quantile for filtering
            
        Returns:
            Alignment parameters
        """
        if isinstance(ms2pip_input, str):
            df = pd.read_csv(ms2pip_input)
        else:
            df = ms2pip_input
        
        # Validate columns
        required_columns = {'seq', 'modifications', 'tr'}
        if not required_columns.issubset(set(df.columns)):
            raise ValueError(f"Input DataFrame must contain columns: {required_columns}")
        
        # Clean data
        df['modifications'] = df['modifications'].replace({np.nan: ""})
        df['seq'] = df['seq'].astype(str)
        df['modifications'] = df['modifications'].astype(str)
        df['tr'] = df['tr'].astype(float)
        
        # Convert to ProForma
        ms2_pip_tuples: list[tuple[str, str]] = list(
            zip(df['seq'].tolist(), df['modifications'].tolist())
        )
        sequences: list[str] = from_ms2_pip(ms2_pip_tuples)
        observed_rts = df['tr'].tolist()
        
        return self.align_to_dataset(
            sequences=sequences,
            observed_rts=observed_rts,
            use_quantile_filtering=use_quantile_filtering,
            lower_quantile=lower_quantile,
            upper_quantile=upper_quantile,
        )
        
    
    def _initialize_model(self, model: str | RtModelType | XGBRegressor) -> XGBRegressor:
        """Load or set the XGBoost model."""
        if isinstance(model, XGBRegressor):
            return model
        elif isinstance(model, RtModelType):
            return self.load_predictor(RT_MODEL_PATHS[model])
        elif isinstance(model, str):
            return self.load_predictor(model)
        else:
            raise ValueError(f"Invalid model type: {type(model)}")
    
    @staticmethod
    def train_new_predictor(
        sequences: Sequence[str | ProFormaAnnotation],
        retention_times: Sequence[float],
        max_retention_time: float | None = None,
        **xgb_params: dict[str, Any],
    ) -> "RtPredictor":
        """
        Train a new XGBoost model.
        
        Args:
            sequences: Peptide sequences or ProForma annotations
            retention_times: Target retention times
            **xgb_params: Additional parameters for XGBRegressor
            
        Returns:
            New RtPredictor with trained model
        """
        X = RtPredictor._sequences_to_features(sequences)
        y = np.array(retention_times)

        if max_retention_time is not None:
            y = y / max_retention_time
        else:
            max_retention_time = np.max(y)
            y = y / max_retention_time


        model = XGBRegressor(**xgb_params)
        model.fit(X, y)
        
        return RtPredictor(model=model)
    
    @staticmethod
    def load_predictor(predictor_path: str | Path) -> XGBRegressor:
        """Load XGBoost model from file."""
        if not os.path.exists(predictor_path):
            raise FileNotFoundError(f"Model file not found: {predictor_path}")
        
        model = XGBRegressor()
        model.load_model(str(predictor_path))
        return model
    
    def save_predictor(self, predictor_path: str | Path) -> None:
        """Save XGBoost model to file."""
        predictor_path = Path(predictor_path)
        predictor_path.parent.mkdir(parents=True, exist_ok=True)
        self.model.save_model(str(predictor_path))
    
    @staticmethod
    def _sequences_to_features(
        sequences: Sequence[str | ProFormaAnnotation]
    ) -> np.ndarray:
        """Convert sequences to feature vectors."""
        feature_vectors = [f.vector for f in generate_feature_vector(sequences)]
        return np.array(feature_vectors)
        

    def finetune(
        self,
        sequences: Sequence[str | ProFormaAnnotation],
        retention_times: Sequence[float],
        max_retention_time: float | None = None,
        **xgb_fit_params: dict[str, Any],
    ) -> None:
        """
        Fine-tune the existing model with new data.
        
        Args:
            sequences: Peptide sequences or ProForma annotations
            retention_times: Target retention times
            **xgb_fit_params: Additional parameters for XGBRegressor.fit()
        """
        X = self._sequences_to_features(sequences)
        y = np.array(retention_times)
        if max_retention_time is not None:
            y = y / max_retention_time
        else:
            max_retention_time = np.max(y)
            y = y / max_retention_time

        self.model.fit(X, y, **xgb_fit_params)


    @staticmethod
    def train_sage(
        sage_input: str | pd.DataFrame,
        train_test_split: float = 0.8,
        qvalue_threshold: float = 0.01,
        qvalue_column: str = 'peptide_q',
        **xgb_params: dict[str, Any]
    ) -> "RtPredictor":
        """
        Train the model using SAGE input data with RT alignment across files.
        
        Args:
            sage_input: Path to SAGE input TSV file or a pandas DataFrame
            max_retention_time: Optional max RT for normalization
            train_test_split: Fraction of data to use for training
            qvalue_threshold: Q-value threshold for filtering
            qvalue_column: Column name for q-value filtering
            **xgb_params: Additional parameters for XGBRegressor
            
        Returns:
            Trained RtPredictor instance
        """
        if isinstance(sage_input, str):
            df = load_sage_files([sage_input])
        else:
            df = sage_input

        # get min/max rt for each file
        min_rt_per_file = df.groupby('filename')['rt'].min()
        max_rt_per_file = df.groupby('filename')['rt'].max()

        # scale rt to 0-1 per file
        def scale_rt(row):
            min_rt = min_rt_per_file[row['filename']]
            max_rt = max_rt_per_file[row['filename']]
            if max_rt > min_rt:
                return (row['rt'] - min_rt) / (max_rt - min_rt)
            else:
                return 0.5  # If all RTs are the same, assign middle value

        df['rt'] = df.apply(scale_rt, axis=1)

        # Filter by spectrum q-value
        df = df[df[qvalue_column] <= qvalue_threshold]

        # for each peptide, filename take the rt with the best hyperscore
        df = df.loc[df.groupby(['peptide', 'filename'])['hyperscore'].idxmax()]

        # for each peptide calculate the median rt
        df['rt'] = df.groupby('peptide')['rt'].transform('median')
        df['stripped_sequence'] = strip_mods(df['peptide'].tolist())

        # train test split on stripped sequence
        unique_sequences = df['stripped_sequence'].unique()
        np.random.shuffle(unique_sequences)
        split_index = int(len(unique_sequences) * train_test_split)
        train_sequences = set(unique_sequences[:split_index])
        test_sequences = set(unique_sequences[split_index:])

        train_df = df[df['stripped_sequence'].isin(train_sequences)]
        test_df = df[df['stripped_sequence'].isin(test_sequences)]

        print(f"Training on {len(train_df)} PSMs, testing on {len(test_df)} PSMs")

        # Train model

        xgboost_model = XGBRegressor(**xgb_params)

        rt_model = RtPredictor(xgboost_model)
        X_train = RtPredictor._sequences_to_features(train_df['peptide'].tolist())
        y_train = train_df['rt'].values

        X_test = RtPredictor._sequences_to_features(test_df['peptide'].tolist())
        y_test = test_df['rt'].values

        rt_model.model = XGBRegressor(**xgb_params)
        rt_model.model.fit(X_train, y_train, eval_set=[(X_test, y_test)], verbose=True)

        print(df.head())

        return rt_model
    
    def finetune_sage(
        self,
        sage_input: str | pd.DataFrame,
        qvalue_threshold: float = 0.01,
        qvalue_column: str = 'peptide_q',
        **xgb_fit_params: dict[str, Any],
    ) -> None:
        """
        Fine-tune the existing model with SAGE input data.
        
        Args:
            sage_input: Path to SAGE input TSV file or a pandas DataFrame
            qvalue_threshold: Q-value threshold for filtering
            qvalue_column: Column name for q-value filtering
            **xgb_fit_params: Additional parameters for XGBRegressor.fit()
        """
        if isinstance(sage_input, str):
            df = load_sage_files([sage_input])
        else:
            df = sage_input

        # Filter by spectrum q-value
        df = df[df[qvalue_column] <= qvalue_threshold]

        # for each peptide, filename take the rt with the best hyperscore
        df = df.loc[df.groupby(['peptide', 'filename'])['hyperscore'].idxmax()]

        # for each peptide calculate the median rt
        df['rt'] = df.groupby('peptide')['rt'].transform('median')

        # min_max scaling of rt
        min_rt = df['rt'].min()
        max_rt = df['rt'].max()
        df['rt'] = (df['rt'] - min_rt) / (max_rt - min_rt)

        sequences = df['peptide'].tolist()
        retention_times = df['rt'].tolist()

        self.finetune(
            sequences,
            retention_times,
            **xgb_fit_params,
        )
    
    @staticmethod
    def train_ms2_pip(ms2pip_input: pd.DataFrame | str, **xgb_params: dict[str, Any]) -> "RtPredictor":
        """
        Train the model using MS2PIP input data.
        
        Args:
            df: DataFrame with 'seq', 'modifications' and 'tr' columns
            **xgb_params: Additional parameters for XGBRegressor
            
        Returns:
            Trained RtPredictor instance
        """

        if isinstance(ms2pip_input, str):
            df = pd.read_csv(ms2pip_input)
        else:
            df = ms2pip_input
        # ensure seq, modifications, tr columns exist
        required_columns = {'seq', 'modifications', 'tr'}
        if not required_columns.issubset(set(df.columns)):
            raise ValueError(f"Input DataFrame must contain columns: {required_columns}")
        
        df['modifications'] = df['modifications'].replace({np.nan: ""})

        # ensure seq and mod cols are str
        df['seq'] = df['seq'].astype(str)
        df['modifications'] = df['modifications'].astype(str)
        df['tr'] = df['tr'].astype(float)

        # min_max scaling of tr
        min_tr = df['tr'].min()
        max_tr = df['tr'].max()
        df['tr'] = (df['tr'] - min_tr) / (max_tr - min_tr)

        # convert any nan mods to ""


        ms2_pip_tuples: list[tuple[str, str]] = list(zip(df['seq'].tolist(), df['modifications'].tolist()))

        rts = df['tr'].tolist()

        sequences: list[str] = from_ms2_pip(ms2_pip_tuples, static_mods={'C': 57.02146})
        X = RtPredictor._sequences_to_features(sequences)
        y = np.array(rts)
        model = XGBRegressor(**xgb_params)
        model.fit(X, y)
        return RtPredictor(model=model)
    

    def fine_tune_ms2_pip(
        self,
        ms2pip_input: pd.DataFrame | str,
        **xgb_fit_params: dict[str, Any],
    ) -> None:
        """
        Fine-tune the existing model with MS2PIP input data.
        
        Args:
            df: DataFrame with 'seq', 'modifications' and 'tr' columns
            **xgb_fit_params: Additional parameters for XGBRegressor.fit()
        """

        if isinstance(ms2pip_input, str):
            df = pd.read_csv(ms2pip_input)
        else:
            df = ms2pip_input
        # ensure seq, modifications, tr columns exist
        required_columns = {'seq', 'modifications', 'tr'}
        if not required_columns.issubset(set(df.columns)):
            raise ValueError(f"Input DataFrame must contain columns: {required_columns}")
        
        df['modifications'] = df['modifications'].replace({np.nan: ""})

        # ensure seq and mod cols are str
        df['seq'] = df['seq'].astype(str)
        df['modifications'] = df['modifications'].astype(str)
        df['tr'] = df['tr'].astype(float)

        ms2_pip_tuples: list[tuple[str, str]] = list(zip(df['seq'].tolist(), df['modifications'].tolist()))

        rts = df['tr'].tolist()

        sequences: list[str] = from_ms2_pip(ms2_pip_tuples, static_mods={'C': 57.02146})
        X = RtPredictor._sequences_to_features(sequences)
        y = np.array(rts)

        self.model.fit(X, y, **xgb_fit_params)

