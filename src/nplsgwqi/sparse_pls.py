from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

from .compositional import clr_transform


def _replace_nonpositive(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        pos = out.loc[out[col] > 0, col]
        floor = pos.min() / np.sqrt(2) if not pos.empty else 1e-6
        out.loc[out[col] <= 0, col] = floor
    return out


def prepare_clr_predictors(comp_data: pd.DataFrame, endpoint_name: str | None = None, include_scale: bool = True) -> pd.DataFrame:
    predictors = comp_data.copy()
    if endpoint_name in predictors.columns:
        predictors = predictors.drop(columns=[endpoint_name])
    predictors = _replace_nonpositive(predictors)
    
    clr_df = clr_transform(predictors)
    if include_scale:
        # Log of geometric mean captures absolute concentration scale
        scale_val = np.mean(np.log(predictors), axis=1)
        clr_df["Total_Mineralization_Scale"] = scale_val
        
    return clr_df


def default_keep_grid(n_features: int) -> list[int]:
    grid = {1, 2, 3, 4, 5, max(2, n_features // 2), n_features}
    return sorted(k for k in grid if 1 <= k <= n_features)


def default_component_grid(n_features: int, n_samples: int) -> list[int]:
    max_comp = max(1, min(4, n_features, n_samples - 1))
    return list(range(1, max_comp + 1))


def compute_vip(pls: PLSRegression, x_scaled: np.ndarray) -> np.ndarray:
    t_scores = np.asarray(pls.x_scores_)
    w_weights = np.asarray(pls.x_weights_)
    q_loadings = np.asarray(pls.y_loadings_)

    if q_loadings.ndim == 1:
        q_loadings = q_loadings.reshape(1, -1)

    n_predictors, n_components = w_weights.shape
    explained = np.diag(t_scores.T @ t_scores @ q_loadings.T @ q_loadings)
    total_explained = float(explained.sum())
    if total_explained <= 0:
        return np.zeros(n_predictors)

    vip = np.zeros(n_predictors)
    for pred_idx in range(n_predictors):
        weight = 0.0
        for comp_idx in range(n_components):
            denom = float(np.sum(w_weights[:, comp_idx] ** 2))
            if denom <= 0:
                continue
            weight += explained[comp_idx] * (w_weights[pred_idx, comp_idx] ** 2) / denom
        vip[pred_idx] = np.sqrt(n_predictors * weight / total_explained)
    return vip


@dataclass
class SparsePLSFit:
    selected_features: list[str]
    n_components: int
    keep_k: int
    x_scaler: StandardScaler
    y_scaler: StandardScaler
    model: PLSRegression
    coefficients: pd.Series
    vip_scores: pd.Series
    driver_contributions: pd.DataFrame
    fitted_values: pd.Series
    in_sample_r2: float

    def predict(self, x_df: pd.DataFrame) -> np.ndarray:
        x_sel = x_df.loc[:, self.selected_features]
        x_scaled = self.x_scaler.transform(x_sel)
        y_scaled = self.model.predict(x_scaled).reshape(-1, 1)
        return self.y_scaler.inverse_transform(y_scaled).ravel()


def fit_sparse_pls_model(
    x_df: pd.DataFrame,
    y: pd.Series,
    n_components: int,
    keep_k: int,
) -> SparsePLSFit:
    keep_k = max(1, min(keep_k, x_df.shape[1]))
    n_components = max(1, min(n_components, keep_k, x_df.shape[1], len(x_df) - 1))

    x_scaler_full = StandardScaler()
    x_scaled_full = x_scaler_full.fit_transform(x_df)
    y_scaler = StandardScaler()
    y_scaled = y_scaler.fit_transform(y.to_numpy().reshape(-1, 1)).ravel()

    dense_pls = PLSRegression(n_components=n_components, scale=False)
    dense_pls.fit(x_scaled_full, y_scaled)
    dense_vip = compute_vip(dense_pls, x_scaled_full)
    top_idx = np.argsort(dense_vip)[::-1][:keep_k]
    selected_features = list(x_df.columns[top_idx])
    selected_features.sort(key=lambda col: (-dense_vip[x_df.columns.get_loc(col)], col))

    x_selected = x_df.loc[:, selected_features]
    x_scaler = StandardScaler()
    x_scaled = x_scaler.fit_transform(x_selected)
    final_components = max(1, min(n_components, x_scaled.shape[1], len(x_selected) - 1))
    sparse_pls = PLSRegression(n_components=final_components, scale=False)
    sparse_pls.fit(x_scaled, y_scaled)

    y_pred_scaled = sparse_pls.predict(x_scaled).reshape(-1, 1)
    y_pred = y_scaler.inverse_transform(y_pred_scaled).ravel()
    coef_selected = np.asarray(sparse_pls.coef_).reshape(-1)
    vip_selected = compute_vip(sparse_pls, x_scaled)

    coef_full = pd.Series(0.0, index=x_df.columns, dtype=float)
    coef_full.loc[selected_features] = coef_selected
    vip_full = pd.Series(0.0, index=x_df.columns, dtype=float)
    vip_full.loc[selected_features] = vip_selected

    y_scale = float(y_scaler.scale_[0]) if float(y_scaler.scale_[0]) != 0 else 1.0
    contrib = pd.DataFrame(0.0, index=x_df.index, columns=x_df.columns, dtype=float)
    contrib.loc[:, selected_features] = x_scaled * coef_selected * y_scale
    contrib["Fitted"] = y_pred
    contrib["Observed"] = y.to_numpy()
    contrib["Residual"] = y.to_numpy() - y_pred

    return SparsePLSFit(
        selected_features=selected_features,
        n_components=final_components,
        keep_k=keep_k,
        x_scaler=x_scaler,
        y_scaler=y_scaler,
        model=sparse_pls,
        coefficients=coef_full,
        vip_scores=vip_full,
        driver_contributions=contrib,
        fitted_values=pd.Series(y_pred, index=y.index, name=f"{y.name}_Fitted"),
        in_sample_r2=float(r2_score(y, y_pred)),
    )


def _safe_kfold(n_samples: int, requested_splits: int, random_state: int) -> KFold:
    n_splits = max(2, min(requested_splits, n_samples))
    return KFold(n_splits=n_splits, shuffle=True, random_state=random_state)


def tune_sparse_pls(
    x_df: pd.DataFrame,
    y: pd.Series,
    component_grid: list[int] | None = None,
    keep_grid: list[int] | None = None,
    n_splits: int = 3,
    random_state: int = 42,
) -> dict:
    component_grid = component_grid or default_component_grid(x_df.shape[1], len(x_df))
    keep_grid = keep_grid or default_keep_grid(x_df.shape[1])
    inner_cv = _safe_kfold(len(x_df), n_splits, random_state)

    rows = []
    best_params = None
    best_score = -np.inf
    best_penalty = None

    for n_components in component_grid:
        for keep_k in keep_grid:
            if keep_k < n_components:
                continue
            fold_scores = []
            for fold_idx, (train_idx, test_idx) in enumerate(inner_cv.split(x_df), start=1):
                x_train, x_test = x_df.iloc[train_idx], x_df.iloc[test_idx]
                y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
                try:
                    fit = fit_sparse_pls_model(x_train, y_train, n_components=n_components, keep_k=keep_k)
                    pred = fit.predict(x_test)
                    score = float(r2_score(y_test, pred))
                except Exception:
                    score = -np.inf
                fold_scores.append(score)
                rows.append(
                    {
                        "n_components": n_components,
                        "keep_k": keep_k,
                        "fold": fold_idx,
                        "r2": score,
                    }
                )
            mean_score = float(np.mean(fold_scores))
            penalty = (keep_k, n_components)
            if mean_score > best_score or (np.isclose(mean_score, best_score) and penalty < best_penalty):
                best_score = mean_score
                best_penalty = penalty
                best_params = {"n_components": n_components, "keep_k": keep_k}

    return {
        "best_params": best_params,
        "grid_results": pd.DataFrame(rows),
        "best_inner_r2_mean": float(best_score),
    }


def nested_cv_sparse_pls(
    x_df: pd.DataFrame,
    y: pd.Series,
    outer_splits: int = 5,
    inner_splits: int = 3,
    component_grid: list[int] | None = None,
    keep_grid: list[int] | None = None,
    random_state: int = 42,
) -> dict:
    component_grid = component_grid or default_component_grid(x_df.shape[1], len(x_df))
    keep_grid = keep_grid or default_keep_grid(x_df.shape[1])
    outer_cv = _safe_kfold(len(x_df), outer_splits, random_state)

    outer_scores = []
    selected_counts = pd.Series(0, index=x_df.columns, dtype=float)
    fold_rows = []
    chosen_params = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(x_df), start=1):
        x_train, x_test = x_df.iloc[train_idx], x_df.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        tuned = tune_sparse_pls(
            x_train,
            y_train,
            component_grid=component_grid,
            keep_grid=keep_grid,
            n_splits=inner_splits,
            random_state=random_state + fold_idx,
        )
        params = tuned["best_params"]
        fit = fit_sparse_pls_model(
            x_train,
            y_train,
            n_components=params["n_components"],
            keep_k=params["keep_k"],
        )
        pred = fit.predict(x_test)
        score = float(r2_score(y_test, pred))
        outer_scores.append(score)
        selected_counts.loc[fit.selected_features] += 1
        chosen_params.append(params)
        fold_rows.append(
            {
                "fold": fold_idx,
                "r2": score,
                "n_components": params["n_components"],
                "keep_k": params["keep_k"],
            }
        )

    return {
        "outer_r2_mean": float(np.mean(outer_scores)),
        "outer_r2_std": float(np.std(outer_scores)),
        "fold_results": pd.DataFrame(fold_rows),
        "outer_selection_frequency": (selected_counts / len(fold_rows)) * 100.0,
        "chosen_params": chosen_params,
    }


def fixed_cv_sparse_pls(
    x_df: pd.DataFrame,
    y: pd.Series,
    n_components: int,
    keep_k: int,
    n_splits: int = 3,
    random_state: int = 42,
) -> dict:
    cv = _safe_kfold(len(x_df), n_splits, random_state)
    scores = []
    selected_counts = pd.Series(0, index=x_df.columns, dtype=float)

    for train_idx, test_idx in cv.split(x_df):
        x_train, x_test = x_df.iloc[train_idx], x_df.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        fit = fit_sparse_pls_model(x_train, y_train, n_components=n_components, keep_k=keep_k)
        pred = fit.predict(x_test)
        scores.append(float(r2_score(y_test, pred)))
        selected_counts.loc[fit.selected_features] += 1

    return {
        "cv_r2_mean": float(np.mean(scores)),
        "cv_r2_std": float(np.std(scores)),
        "selection_frequency": (selected_counts / len(scores)) * 100.0,
    }


def bootstrap_stability_selection(
    x_df: pd.DataFrame,
    y: pd.Series,
    n_components: int,
    keep_k: int,
    n_bootstraps: int = 100,
    sample_fraction: float = 0.8,
    random_state: int = 42,
) -> pd.Series:
    rng = np.random.default_rng(random_state)
    counts = pd.Series(0, index=x_df.columns, dtype=float)
    sample_size = max(5, int(np.ceil(sample_fraction * len(x_df))))

    for _ in range(n_bootstraps):
        sample_idx = rng.choice(len(x_df), size=sample_size, replace=True)
        fit = fit_sparse_pls_model(
            x_df.iloc[sample_idx],
            y.iloc[sample_idx],
            n_components=n_components,
            keep_k=keep_k,
        )
        counts.loc[fit.selected_features] += 1

    return (counts / n_bootstraps) * 100.0


def permutation_test_sparse_pls(
    x_df: pd.DataFrame,
    y: pd.Series,
    observed_score: float,
    n_components: int,
    keep_k: int,
    outer_splits: int = 3,
    n_permutations: int = 50,
    random_state: int = 42,
) -> tuple[float, pd.Series]:
    rng = np.random.default_rng(random_state)
    perm_scores = []

    for perm_idx in range(n_permutations):
        permuted = pd.Series(rng.permutation(y.to_numpy()), index=y.index, name=y.name)
        cv_res = fixed_cv_sparse_pls(
            x_df,
            permuted,
            n_components=n_components,
            keep_k=keep_k,
            n_splits=outer_splits,
            random_state=random_state + perm_idx + 1,
        )
        perm_scores.append(cv_res["cv_r2_mean"])

    perm_series = pd.Series(perm_scores, name="Permutation_R2")
    p_value = float((1 + np.sum(perm_series >= observed_score)) / (n_permutations + 1))
    return p_value, perm_series


def run_sparse_pls_endpoint_model(
    comp_data: pd.DataFrame,
    target: pd.Series,
    endpoint_name: str,
    *,
    include_scale: bool = True,
    use_augmented: bool = False,
    outer_splits: int = 5,
    inner_splits: int = 3,
    n_bootstraps: int = 100,
    n_permutations: int = 50,
    random_state: int = 42,
) -> dict:
    if use_augmented:
        x_df = comp_data.copy()
        block_label = "Augmented Field+Compositional"
    else:
        x_df = prepare_clr_predictors(
            comp_data, 
            endpoint_name=endpoint_name if endpoint_name in comp_data.columns else None,
            include_scale=include_scale
        )
        block_label = "Scale-Aware CLR chemistry"
        if endpoint_name in comp_data.columns:
            block_label = f"{endpoint_name}-safe {block_label}"
            
    component_grid = default_component_grid(x_df.shape[1], len(x_df))
    keep_grid = default_keep_grid(x_df.shape[1])

    tuned = tune_sparse_pls(
        x_df,
        target,
        component_grid=component_grid,
        keep_grid=keep_grid,
        n_splits=inner_splits,
        random_state=random_state,
    )
    best_params = tuned["best_params"]
    final_fit = fit_sparse_pls_model(
        x_df,
        target,
        n_components=best_params["n_components"],
        keep_k=best_params["keep_k"],
    )
    nested = nested_cv_sparse_pls(
        x_df,
        target,
        outer_splits=outer_splits,
        inner_splits=inner_splits,
        component_grid=component_grid,
        keep_grid=keep_grid,
        random_state=random_state,
    )
    stability = bootstrap_stability_selection(
        x_df,
        target,
        n_components=best_params["n_components"],
        keep_k=best_params["keep_k"],
        n_bootstraps=n_bootstraps,
        random_state=random_state,
    )
    permutation_p, permutation_scores = permutation_test_sparse_pls(
        x_df,
        target,
        observed_score=nested["outer_r2_mean"],
        n_components=best_params["n_components"],
        keep_k=best_params["keep_k"],
        outer_splits=outer_splits,
        n_permutations=n_permutations,
        random_state=random_state,
    )

    driver_raw = stability.clip(lower=0.0) * final_fit.vip_scores.clip(lower=0.0)
    if float(driver_raw.sum()) > 0:
        driver_score = (driver_raw / driver_raw.sum()) * 100.0
    else:
        driver_score = pd.Series(0.0, index=x_df.columns)

    return {
        "endpoint": endpoint_name,
        "transform": "CLR+Scale",
        "predictor_columns": list(x_df.columns),
        "predictor_block_label": block_label,
        "best_n_components": int(best_params["n_components"]),
        "best_keep_k": int(best_params["keep_k"]),
        "in_sample_r2": final_fit.in_sample_r2,
        "nested_cv_r2_mean": nested["outer_r2_mean"],
        "nested_cv_r2_std": nested["outer_r2_std"],
        "selection_frequency": stability,
        "outer_selection_frequency": nested["outer_selection_frequency"],
        "vip_scores": final_fit.vip_scores,
        "coefficients": final_fit.coefficients,
        "driver_score": driver_score.sort_values(ascending=False),
        "selected_features": final_fit.selected_features,
        "sample_level_contributions": final_fit.driver_contributions,
        "fitted_values": final_fit.fitted_values,
        "permutation_p_value": permutation_p,
        "permutation_scores": permutation_scores,
        "n_permutations": int(n_permutations),
        "grid_results": tuned["grid_results"],
    }
