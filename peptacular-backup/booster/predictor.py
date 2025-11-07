from typing import TYPE_CHECKING, Literal, Sequence

if TYPE_CHECKING:
    import numpy as np

from ...proforma.annotation import ProFormaAnnotation
from .adapters import PeptideData
from .feature_vector import generate_feature_vector
from .utils import train_test_split

from xgboost import XGBRegressor
from typing import NamedTuple

class Score(NamedTuple):
    """Scoring metrics for retention time predictions."""
    r2: float
    mae: float
    rmse: float
    mse: float
    median_ae: float
    n_samples: int

def sequences_to_features(
    sequences: Sequence[str | ProFormaAnnotation]
) -> "np.ndarray":
    """Convert sequences to feature vectors."""
    import numpy as np
    feature_vectors = [f.vector for f in generate_feature_vector(sequences)]
    return np.array(feature_vectors)



class RtPredictor:
    """Retention time predictor using XGBoost."""
    
    def __init__(self, **xgb_params):
        """
        Initialize RT predictor with XGBoost.
        
        Args:
            **xgb_params: Parameters for XGBRegressor
        """

        # if early_stopping is not mentioned, add it
        if 'early_stopping_rounds' not in xgb_params:
            xgb_params['early_stopping_rounds'] = 5

        self.estimator = XGBRegressor(**xgb_params)
        self.alignment_slope_ = None
        self.alignment_intercept_ = None
    
    def _fit(self, data: PeptideData, **fit_params):
        """
        Fit the model.
        
        Args:
            data: PeptideData with sequences and retention times
            **fit_params: Additional parameters for XGBRegressor.fit()
        """
        import numpy as np
        
        X = sequences_to_features(data.sequences)
        y = np.array(data.retention_times)
        
        self.estimator.fit(X, y, **fit_params)
        
        return self
        
    def fit(self, data: PeptideData, test_ratio: float = 0.15, early_stopping: int | None = 5, **fit_params) -> "RtPredictor":
        import numpy as np
        
        train_data, test_data = train_test_split(
            data,
            test_size=test_ratio,
            stratify_by_sequence=True
        )

        # Set up early stopping with eval set
        if early_stopping is not None:
            fit_params = fit_params.copy()
            fit_params['eval_set'] = [
                (sequences_to_features(test_data.sequences), 
                np.array(test_data.retention_times))
            ]
            fit_params.setdefault('verbose', False)

        self._fit(train_data, **fit_params)
        return self
    
    def finetune(self, data: PeptideData, learning_rate: float = 0.01, test_ratio: float = 0.15,
                    n_estimators_add: int = 50,
                early_stopping: int | None = 5, **fit_params) -> "RtPredictor":
        """
        Fine-tune the model on new data by adding additional trees with low learning rate.
        This preserves the existing model and makes gentle adjustments.
        
        Args:
            data: PeptideData with sequences and retention times
            learning_rate: Low learning rate for gentle updates (default: 0.01)
            n_estimators_add: Number of additional boosting rounds to add (default: 50)
            **fit_params: Additional parameters for XGBRegressor.fit()
        """
        import numpy as np
        
        # Store original parameters and update for fine-tuning
        original_params = self.estimator.get_params()
        finetune_params = {
            'learning_rate': learning_rate,
            'n_estimators': n_estimators_add,  # Only the NEW trees
            # Copy other params except n_estimators
            **{k: v for k, v in self.estimator.get_params().items() 
            if k != 'n_estimators' and k != 'learning_rate'}
        }
        # Create new estimator that will continue from existing model
        new_estimator = XGBRegressor(**finetune_params)

        train_data, test_data = train_test_split(
            data,
            test_size=test_ratio,
            stratify_by_sequence=True
        )

        X = sequences_to_features(train_data.sequences)
        y = np.array(train_data.retention_times)
        # Set up early stopping with eval set
        if early_stopping is not None:
            fit_params = fit_params.copy()
            fit_params['eval_set'] = [
                (sequences_to_features(test_data.sequences), 
                np.array(test_data.retention_times))
            ]
            fit_params.setdefault('verbose', False)

        # Continue training from existing model (adds trees, doesn't replace)
        fit_params['xgb_model'] = self.estimator.get_booster()
        new_estimator.fit(X, y, **fit_params)
        
        self.estimator = new_estimator
        
        return self

    def predict(self, sequences: list[str] | PeptideData, apply_alignment: bool = True, clip: tuple[float, float] | None = None) -> list[float]:
        """
        Predict retention times.
        
        Args:
            sequences: List of sequences (str/ProForma) or PeptideData
            apply_alignment: Whether to apply alignment
            gradient: Optional gradient multiplier (only works without alignment)
            
        Returns:
            List of predicted retention times
        """
        import numpy as np
        
        if isinstance(sequences, PeptideData):
            sequences = sequences.sequences
        
        X = sequences_to_features(sequences)
        predictions = self.estimator.predict(X)
        
        # Apply alignment if available
        if apply_alignment and self.alignment_slope_ is not None:
            print(f"Applying linear alignment to predictions: slope={self.alignment_slope_}, intercept={self.alignment_intercept_}")
            predictions = predictions * self.alignment_slope_ + self.alignment_intercept_
            

        # Apply clipping if specified
        if clip is not None:
            predictions = np.clip(predictions, clip[0], clip[1])
        
        return predictions.tolist()
    
    def align_to_observed(
        self, 
        data: PeptideData,
        remove_outliers: bool = True,
        outlier_method: str = "iqr",
        outlier_threshold: float = 1.5,
        density_bins: int | None = 20,
        samples_per_bin: int | None = 50
    ):
        """
        Fit linear alignment to observed data.
        
        Args:
            data: PeptideData with sequences and observed retention times
            remove_outliers: Whether to remove outliers before fitting
            outlier_method: 'iqr', 'std', or 'percentile'
            outlier_threshold: Threshold for outlier detection
            density_bins: Number of bins to divide RT range into for even sampling (None = no binning)
            samples_per_bin: Max samples per bin (None = no limit, uses all data in each bin)
        
        Returns:
            Self for method chaining
        """
        import numpy as np
        
        predicted = np.array(self.predict(data.sequences, apply_alignment=False))
        observed = np.array(data.retention_times)
        
        # Remove outliers if requested
        if remove_outliers:
            mask = self._get_outlier_mask(predicted, observed, outlier_method, outlier_threshold)
            predicted = predicted[mask]
            observed = observed[mask]
            
            if len(predicted) < 2:
                raise ValueError(
                    f"Too few points remaining after outlier removal ({len(predicted)}). "
                    "Try a less aggressive threshold or disable outlier removal."
                )
        
        # Apply density-based sampling if requested
        if density_bins is not None and density_bins > 0:
            # Bin by observed RT to ensure even coverage across the RT range
            bin_edges = np.linspace(observed.min(), observed.max(), density_bins + 1)
            bin_indices = np.digitize(observed, bin_edges) - 1
            bin_indices = np.clip(bin_indices, 0, density_bins - 1)  # Handle edge cases
            
            selected_indices = []
            for bin_idx in range(density_bins):
                bin_mask = bin_indices == bin_idx
                bin_points = np.where(bin_mask)[0]
                
                if len(bin_points) > 0:
                    if samples_per_bin is not None and len(bin_points) > samples_per_bin:
                        # Randomly sample from this bin
                        sampled = np.random.choice(bin_points, samples_per_bin, replace=False)
                        selected_indices.extend(sampled)
                    else:
                        # Use all points in this bin
                        selected_indices.extend(bin_points)
            
            selected_indices = np.array(selected_indices)
            predicted = predicted[selected_indices]
            observed = observed[selected_indices]
            
            print(f"Density sampling: selected {len(selected_indices)} points from {density_bins} bins")
        
        if len(predicted) < 2:
            raise ValueError(
                f"Too few points remaining after density sampling ({len(predicted)}). "
                "Try fewer bins or more samples per bin."
            )
        
        # Fit linear alignment
        self.alignment_slope_, self.alignment_intercept_ = np.polyfit(predicted, observed, 1)
        
        return self
    
    def _get_outlier_mask(self, predicted, observed, method, threshold):
        """Get boolean mask for non-outlier points."""
        import numpy as np
        
        residuals = np.abs(predicted - observed)
        
        if method == "iqr":
            q1 = np.percentile(residuals, 25)
            q3 = np.percentile(residuals, 75)
            iqr = q3 - q1
            threshold_val = q3 + threshold * iqr
            return residuals <= threshold_val
        
        elif method == "std":
            mean_residual = np.mean(residuals)
            std_residual = np.std(residuals)
            threshold_val = mean_residual + threshold * std_residual
            return residuals <= threshold_val
        
        elif method == "percentile":
            threshold_val = np.percentile(residuals, threshold)
            return residuals <= threshold_val
        
        else:
            raise ValueError(f"Unknown outlier_method '{method}'")
    
    def clear_alignment(self):
        """Remove alignment parameters."""
        self.alignment_slope_ = None
        self.alignment_intercept_ = None
        return self
    
    def save_model(self, path: str):
        """Save XGBoost model to file."""
        self.estimator.save_model(path)
    
    def load_model(self, path: str):
        """Load XGBoost model from file."""
        self.estimator.load_model(path)

    def score_summary(self, data: PeptideData) -> Score:
        """
        Calculate multiple scoring metrics.
        
        Args:
            data: PeptideData with sequences and observed retention times
            
        Returns:
            Score namedtuple with all available metrics
        """
        import numpy as np
        from sklearn.metrics import (
            r2_score,
            mean_absolute_error,
            mean_squared_error,
            median_absolute_error
        )
        
        predictions = np.array(self.predict(data.sequences))
        observed = np.array(data.retention_times)

        # raise error if nan in observed or predictions
        if np.any(np.isnan(observed)):
            raise ValueError(f"Observed retention times contain NaN values: {len(observed[np.isnan(observed)])}")

        if np.any(np.isnan(predictions)):
            raise ValueError(f"Predicted retention times contain NaN values: {len(predictions[np.isnan(predictions)])}")

        return Score(
            r2=r2_score(observed, predictions),
            mae=mean_absolute_error(observed, predictions),
            rmse=np.sqrt(mean_squared_error(observed, predictions)),
            mse=mean_squared_error(observed, predictions),
            median_ae=median_absolute_error(observed, predictions),
            n_samples=len(observed)
        )
    
    def plot_predictions(
        self, 
        data: PeptideData,
        title: str = "Predicted vs Observed Retention Times",
        figsize: tuple[int, int] = (8, 8),
        show_metrics: bool = True,
        ax=None
    ) -> tuple["plt.Figure", "plt.Axes"]:
        """
        Create a scatter plot of predicted vs observed retention times.
        
        Args:
            data: PeptideData with sequences and observed retention times
            title: Plot title
            figsize: Figure size (width, height) in inches
            show_metrics: Whether to display metrics on the plot
            ax: Optional matplotlib axis to plot on
            
        Returns:
            matplotlib Figure and Axis objects
        """
        import numpy as np
        import matplotlib.pyplot as plt
        
        predictions = np.array(self.predict(data.sequences))
        observed = np.array(data.retention_times)
        
        # Create figure if ax not provided
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
        
        # Scatter plot
        ax.scatter(observed, predictions, alpha=0.5, s=20, edgecolors='none')
        
        # Perfect prediction line
        min_val = min(observed.min(), predictions.min())
        max_val = max(observed.max(), predictions.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2)
        
        # Labels and title
        ax.set_xlabel('Observed Retention Time', fontsize=12)
        ax.set_ylabel('Predicted Retention Time', fontsize=12)
        ax.set_title(title, fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add metrics text box
        if show_metrics:
            scores = self.score_summary(data)
            metrics_text = (
                f'R² = {scores.r2:.3f}\n'
                f'MAE = {scores.mae:.2f}\n'
                f'RMSE = {scores.rmse:.2f}\n'
                f'n = {scores.n_samples}'
            )
            ax.text(0.05, 0.95, metrics_text,
                    transform=ax.transAxes,
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                    fontsize=10)
        
        plt.tight_layout()
        return fig, ax
    

    # load and save methods (use opickle)
    def save(self, path: str):
        """Save the RtPredictor instance to a file."""
        import pickle
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, path: str) -> "RtPredictor":
        """Load a RtPredictor instance from a file."""
        import pickle
        with open(path, 'rb') as f:
            instance = pickle.load(f)
        return instance
    

    @classmethod
    def tune_hyperparameters(
        self,
        data: PeptideData,
        n_trials: int = 50,
        test_ratio: float = 0.15,
        cv_folds: int = 3,
        random_state: int | None = 42,
        verbose: bool = True
    ) -> tuple["RtPredictor", dict]:
        """
        Tune hyperparameters using Optuna (Bayesian optimization).
        Much faster and smarter than grid search.
        
        Args:
            data: PeptideData for hyperparameter tuning
            n_trials: Number of trials to run
            test_ratio: Fraction of data to hold out for final evaluation
            cv_folds: Number of cross-validation folds
            random_state: Random seed for reproducibility
            verbose: Whether to print progress
            
        Returns:
            Tuple of (best_predictor, results_dict)
        """
        import optuna
        from sklearn.model_selection import cross_val_score
        import numpy as np
        
        # Split data into train and test
        train_data, test_data = train_test_split(
            data,
            test_size=test_ratio,
            random_state=random_state,
            stratify_by_sequence=True
        )
        
        # Prepare features
        X = sequences_to_features(train_data.sequences)
        y_train = np.array(train_data.retention_times)
        
        def objective(trial):
            """Objective function for Optuna optimization."""
            params = {
                'n_estimators': trial.suggest_int('n_estimators', 100, 1000),
                'max_depth': trial.suggest_int('max_depth', 3, 12),
                'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
                'subsample': trial.suggest_float('subsample', 0.6, 1.0),
                'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
                'min_child_weight': trial.suggest_int('min_child_weight', 1, 10),
                'gamma': trial.suggest_float('gamma', 0, 0.5),
                'reg_alpha': trial.suggest_float('reg_alpha', 0, 1.0),
                'reg_lambda': trial.suggest_float('reg_lambda', 0.5, 5.0),
                'random_state': random_state
            }
            
            model = XGBRegressor(**params)
            
            # Use cross-validation to evaluate
            scores = cross_val_score(
                model, X, y_train,
                cv=cv_folds,
                scoring='neg_mean_absolute_error',
                n_jobs=-1
            )
            
            return -scores.mean()  # Return MAE (lower is better)
        
        # Create study and optimize
        if verbose:
            optuna.logging.set_verbosity(optuna.logging.INFO)
        else:
            optuna.logging.set_verbosity(optuna.logging.WARNING)
        
        study = optuna.create_study(
            direction='minimize',
            sampler=optuna.samplers.TPESampler(seed=random_state)
        )
        
        study.optimize(objective, n_trials=n_trials, show_progress_bar=verbose, n_jobs=10)
        
        # Get best parameters
        best_params = study.best_params

        # Create predictor with best parameters
        best_predictor = RtPredictor(**best_params)
        
        # Train on full training set
        best_predictor.fit(train_data, test_ratio=0.15, early_stopping=5)
        
        # Evaluate on held-out test set
        test_score = best_predictor.score_summary(test_data)
        
        # Prepare results
        results = {
            'best_params': best_params,
            'best_cv_score': study.best_value,
            'test_score': test_score,
            'study': study,
            'n_trials': len(study.trials),
            'optimization_history': [trial.value for trial in study.trials]
        }
        
        if verbose:
            print(f"\nBest parameters: {best_params}")
            print(f"Best CV MAE: {study.best_value:.4f}")
            print(f"\nTest set performance:")
            print(f"  R² = {test_score.r2:.4f}")
            print(f"  MAE = {test_score.mae:.4f}")
            print(f"  RMSE = {test_score.rmse:.4f}")
        
        return best_predictor, results


    def _get_feature_names(self) -> list[str]:
        """Get human-readable feature names."""
        from .feature_vector import PeptideFeatures
        return PeptideFeatures.get_feature_names()

    def explain_predictions(
        self,
        data: PeptideData,
        background_size: int = 100,
        plot_type: Literal['summary', 'bar', 'waterfall', 'force', 'beeswarm'] = "summary",
        max_display: int = 20,
        show_plot: bool = True,
        feature_groups: dict[str, list[str]] | None = None
    ):
        """
        Generate SHAP explanations for predictions with proper feature names.
        
        Args:
            data: PeptideData to explain
            background_size: Number of background samples for SHAP explainer
            plot_type: Type of plot - 'summary', 'bar', 'waterfall', 'force', 'beeswarm'
            max_display: Maximum number of features to display
            show_plot: Whether to display the plot
            feature_groups: Optional dict mapping group names to feature name patterns
            
        Returns:
            Dictionary with SHAP values and explainer
        """
        import shap
        import numpy as np
        import matplotlib.pyplot as plt
        
        # Prepare features
        X = sequences_to_features(data.sequences)
        feature_names = self._get_feature_names()
        
        # Create background dataset for SHAP (sample if too large)
        if len(X) > background_size:
            background_indices = np.random.choice(len(X), background_size, replace=False)
            X_background = X[background_indices]
        else:
            X_background = X
        
        # Create SHAP explainer
        explainer = shap.TreeExplainer(self.estimator, X_background)
        shap_values = explainer(X, check_additivity=False)
        
        # Add feature names to SHAP values
        shap_values.feature_names = feature_names
        
        if show_plot:
            if plot_type == "summary":
                shap.summary_plot(shap_values, X, feature_names=feature_names, 
                                max_display=max_display, show=False)
            elif plot_type == "bar":
                shap.plots.bar(shap_values, max_display=max_display, show=False)
            elif plot_type == "beeswarm":
                shap.plots.beeswarm(shap_values, max_display=max_display, show=False)
            elif plot_type == "waterfall":
                # Show waterfall for first prediction
                shap.plots.waterfall(shap_values[0], max_display=max_display, show=False)
            elif plot_type == "force":
                # Force plot for first prediction
                shap.force_plot(
                    explainer.expected_value,
                    shap_values.values[0],
                    X[0],
                    feature_names=feature_names,
                    matplotlib=True,
                    show=False
                )
            
            plt.tight_layout()
            plt.show()
        
        return {
            'shap_values': shap_values,
            'explainer': explainer,
            'expected_value': explainer.expected_value,
            'feature_names': feature_names
        }
    
    def plot_alignment(self, data: PeptideData, figsize=(10, 6)):
        """
        Visualize the linear alignment.
        
        Args:
            data: PeptideData to plot alignment for
            figsize: Figure size
        """
        import numpy as np
        import matplotlib.pyplot as plt
        
        predicted = np.array(self.predict(data.sequences, apply_alignment=False))
        observed = np.array(data.retention_times)
        aligned = np.array(self.predict(data.sequences, apply_alignment=True))
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Before alignment
        ax1.scatter(predicted, observed, alpha=0.5, s=20)
        ax1.plot([predicted.min(), predicted.max()], 
                 [predicted.min(), predicted.max()], 'r--', label='Perfect')
        ax1.set_xlabel('Predicted RT')
        ax1.set_ylabel('Observed RT')
        ax1.set_title('Before Alignment')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # After alignment
        ax2.scatter(aligned, observed, alpha=0.5, s=20)
        ax2.plot([aligned.min(), aligned.max()], 
                 [aligned.min(), aligned.max()], 'r--', label='Perfect')
        ax2.set_xlabel('Aligned RT')
        ax2.set_ylabel('Observed RT')
        ax2.set_title('After Linear Alignment')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Show alignment line
        if self.alignment_slope_ is not None:
            pred_sorted = np.sort(predicted)
            aligned_sorted = pred_sorted * self.alignment_slope_ + self.alignment_intercept_
            ax1.plot(pred_sorted, aligned_sorted, 'g-', linewidth=2, label='Alignment line')
            ax1.legend()
        
        plt.tight_layout()
        plt.show()