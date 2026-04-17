import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

from .compositional import split_score_regimes
from .process_discovery import ProcessDiscoveryEngine, run_process_discovery
from .process_regression import build_process_regressor, run_process_contribution_model


def bootstrap_process_stability(comp_data: pd.DataFrame, n_components: int = 2, n_bootstraps: int = 100):
    """
    Evaluates stability of process discovery using bootstrapping.
    """
    loadings_history = []

    for _ in range(n_bootstraps):
        boot_idx = np.random.choice(comp_data.index, size=len(comp_data), replace=True)
        boot_data = comp_data.loc[boot_idx]

        try:
            res = run_process_discovery(boot_data, n_components=n_components)
            loadings_history.append(np.abs(res["loadings_clr"].values))
        except Exception:
            continue

    if not loadings_history:
        return None

    stacked = np.stack(loadings_history)
    mean_loadings = np.mean(stacked, axis=0)
    std_loadings = np.std(stacked, axis=0)

    return {
        "mean_loadings": pd.DataFrame(
            mean_loadings,
            columns=comp_data.columns,
            index=[f"Process_{i+1}" for i in range(n_components)],
        ),
        "std_loadings": pd.DataFrame(
            std_loadings,
            columns=comp_data.columns,
            index=[f"Process_{i+1}" for i in range(n_components)],
        ),
    }


def _apply_training_positive_floor(train_df: pd.DataFrame, test_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    train_out = train_df.copy()
    test_out = test_df.copy()
    for col in train_out.columns:
        pos = train_out.loc[train_out[col] > 0, col]
        floor = pos.min() / np.sqrt(2) if not pos.empty else 1e-6
        train_out.loc[train_out[col] <= 0, col] = floor
        test_out.loc[test_out[col] <= 0, col] = floor
    return train_out, test_out


def cross_validate_regression(
    comp_data: pd.DataFrame,
    target: pd.Series,
    n_splits: int = 5,
    use_ridge: bool = False,
    n_components: int | None = None,
    robust: bool = True,
):
    """
    Cross-validates the endpoint model with fold-wise process rediscovery.
    Discovery is fit on each training fold only, then applied to the held-out fold.
    """
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = []
    for train_idx, test_idx in kf.split(comp_data):
        train_comp = comp_data.iloc[train_idx].copy()
        test_comp = comp_data.iloc[test_idx].copy()
        y_train, y_test = target.iloc[train_idx], target.iloc[test_idx]

        if target.name in train_comp.columns:
            train_comp = train_comp.drop(columns=[target.name])
            test_comp = test_comp.drop(columns=[target.name])

        train_comp, test_comp = _apply_training_positive_floor(train_comp, test_comp)

        discovery = ProcessDiscoveryEngine(n_components=n_components, robust=robust)
        train_scores = discovery.fit_transform(train_comp)
        test_scores = discovery.transform(test_comp)

        X_train = split_score_regimes(train_scores)
        X_test = split_score_regimes(test_scores)
        inner_cv = min(5, len(y_train)) if use_ridge and len(y_train) >= 3 else None
        model = build_process_regressor(use_ridge=use_ridge, cv_folds=inner_cv)
        model.fit(X_train, y_train)
        scores.append(float(model.score(X_test, y_test)))

    return {
        "cv_explained_variance_mean": float(np.mean(scores)),
        "cv_explained_variance_std": float(np.std(scores)),
    }


def report_model_diagnostics(regression_results: dict):
    """
    Summarizes the statistical accuracy of the process contribution models.
    """
    diagnostics = []
    for ep, res in regression_results.items():
        diagnostics.append(
            {
                "Endpoint": ep,
                "R2_Explained": res["explained_variance"],
                "Model_Type": res["model_type"],
                "Regularization_Alpha": res.get("regularization_alpha", 0.0),
                "Average_Fitted_Contribution_Sum": res["relative_process_contributions"].drop("Unexplained").sum(),
                "Explained_Variability_Sum": res["process_variability_importance"].sum(),
            }
        )
    return pd.DataFrame(diagnostics)


def check_endpoint_leakage(comp_vars: list, target_var: str) -> bool:
    """
    Checks if the target endpoint is accidentally part of the compositional discovery block.
    """
    return target_var in comp_vars
