import math

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, RidgeCV

from .compositional import split_score_regimes


def build_process_regressor(use_ridge: bool = False, cv_folds: int | None = None):
    if use_ridge:
        ridge_kwargs = {"alphas": np.logspace(-3, 4, 16)}
        if cv_folds is not None and cv_folds >= 2:
            ridge_kwargs["cv"] = cv_folds
        return RidgeCV(**ridge_kwargs)
    return LinearRegression(positive=True)


def infer_process_groups(split_scores: pd.DataFrame) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for col in split_scores.columns:
        if "(" in col:
            group = col.split("(")[0]
        else:
            group = col
        groups.setdefault(group, []).append(col)
    return groups


def score_process_subset(
    split_scores: pd.DataFrame,
    target: pd.Series,
    process_groups: dict[str, list[str]],
    included_groups: frozenset[str],
    use_ridge: bool = False,
    cv_folds: int | None = None,
) -> float:
    selected_cols = [col for group in included_groups for col in process_groups[group]]
    if not selected_cols:
        return 0.0
    model = build_process_regressor(use_ridge=use_ridge, cv_folds=cv_folds)
    model.fit(split_scores[selected_cols], target)
    return max(float(model.score(split_scores[selected_cols], target)), 0.0)


def grouped_shapley_importance(
    split_scores: pd.DataFrame,
    target: pd.Series,
    use_ridge: bool = False,
    cv_folds: int | None = None,
) -> tuple[pd.Series, pd.Series]:
    process_groups = infer_process_groups(split_scores)
    group_names = list(process_groups.keys())
    n_groups = len(group_names)
    if n_groups == 0:
        empty = pd.Series(dtype=float)
        return empty, empty

    coalition_scores: dict[frozenset[str], float] = {frozenset(): 0.0}
    for mask in range(1, 1 << n_groups):
        subset = frozenset(group_names[idx] for idx in range(n_groups) if mask & (1 << idx))
        coalition_scores[subset] = score_process_subset(
            split_scores,
            target,
            process_groups,
            subset,
            use_ridge=use_ridge,
            cv_folds=cv_folds,
        )

    shapley_values: dict[str, float] = {}
    factorial_den = math.factorial(n_groups)
    for group in group_names:
        contribution = 0.0
        others = [g for g in group_names if g != group]
        for mask in range(1 << len(others)):
            subset = frozenset(others[idx] for idx in range(len(others)) if mask & (1 << idx))
            subset_with_group = subset | {group}
            weight = (
                math.factorial(len(subset))
                * math.factorial(n_groups - len(subset) - 1)
                / factorial_den
            )
            contribution += weight * (coalition_scores[subset_with_group] - coalition_scores[subset])
        shapley_values[group] = contribution

    shapley_series = pd.Series(shapley_values).sort_index()
    full_r2 = float(coalition_scores[frozenset(group_names)])
    positive_shapley = shapley_series.clip(lower=0.0)
    positive_sum = float(positive_shapley.sum())
    if positive_sum > 0 and full_r2 > 0:
        explained_pct = (positive_shapley / positive_sum) * 100.0
        total_pct = (positive_shapley / positive_sum) * full_r2 * 100.0
    else:
        explained_pct = pd.Series(0.0, index=shapley_series.index)
        total_pct = pd.Series(0.0, index=shapley_series.index)
    return explained_pct, total_pct


def run_process_contribution_model(
    process_scores: pd.DataFrame,
    target: pd.Series,
    use_ridge: bool = False,
) -> dict:
    """
    Fits constrained regression on split process scores and returns both:
    1. Average fitted contribution shares for compatibility with sample-level outputs.
    2. Grouped Shapley importance shares for process-level variability ranking.
    """
    split_scores = split_score_regimes(process_scores)
    cv_folds = min(5, len(target)) if use_ridge and len(target) >= 3 else None
    model = build_process_regressor(use_ridge=use_ridge, cv_folds=cv_folds)
    model.fit(split_scores, target)

    pole_contributions = np.zeros(split_scores.shape)
    for j in range(split_scores.shape[1]):
        pole_contributions[:, j] = model.coef_[j] * split_scores.iloc[:, j].values

    unexplained = model.intercept_ + (target.values - model.predict(split_scores))
    mean_contributions = np.maximum(0, pole_contributions.mean(axis=0))
    mean_unexplained = np.maximum(0, unexplained.mean())
    total_predicted = mean_contributions.sum() + mean_unexplained
    if total_predicted > 0:
        percent_contributions = (mean_contributions / total_predicted) * 100.0
        percent_unexplained = (mean_unexplained / total_predicted) * 100.0
    else:
        percent_contributions = np.zeros_like(mean_contributions)
        percent_unexplained = 0.0

    contrib_df = pd.Series(percent_contributions, index=split_scores.columns)
    contrib_df["Unexplained"] = percent_unexplained

    sample_level_df = pd.DataFrame(
        pole_contributions,
        columns=split_scores.columns,
        index=process_scores.index,
    )
    sample_level_df["Unexplained"] = unexplained

    explained_importance, total_importance = grouped_shapley_importance(
        split_scores,
        target,
        use_ridge=use_ridge,
        cv_folds=cv_folds,
    )

    return {
        "endpoint": target.name,
        "model_type": "Ridge" if use_ridge else "OLS_Positive",
        "regularization_alpha": float(getattr(model, "alpha_", 0.0)),
        "coefficients": pd.Series(model.coef_, index=split_scores.columns),
        "intercept": model.intercept_,
        "relative_process_contributions": contrib_df,
        "sample_level_contributions": sample_level_df,
        "process_variability_importance": explained_importance,
        "process_total_variability_share": total_importance,
        "explained_variance": float(model.score(split_scores, target)),
        "process_groups": infer_process_groups(split_scores),
    }
