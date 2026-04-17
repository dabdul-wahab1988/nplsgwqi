from __future__ import annotations

import math
import pickle
from pathlib import Path
import sys

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = ROOT / "manuscript"
ARTIFACTS_DIR = MANUSCRIPT_DIR / "artifacts"
MODEL_STORE_DIR = MANUSCRIPT_DIR / "model_store"

SRC_DIR = ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


def _read_data() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "data.csv")
    if "SampleID" in df.columns:
        df = df.set_index("SampleID", drop=False)
    return df


def _load_npls_gwqi_target() -> pd.Series:
    wqi = pd.read_csv(ARTIFACTS_DIR / "TableS1_Full_WQI.csv")
    wqi = wqi.set_index("SampleID")
    return wqi["NPLS_GWQI"].astype(float)


def _standardize(X: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    Xv = X.to_numpy(dtype=float)
    mean = Xv.mean(axis=0)
    std = Xv.std(axis=0, ddof=0)
    std = np.where(std == 0, 1.0, std)
    return (Xv - mean) / std, mean, std


def _standardize_y(y: pd.Series) -> tuple[np.ndarray, float, float]:
    yv = y.to_numpy(dtype=float)
    mean = float(np.mean(yv))
    std = float(np.std(yv, ddof=0))
    std = std if std != 0 else 1.0
    return (yv - mean) / std, mean, std


def _inverse_y(y_scaled: np.ndarray, mean: float, std: float) -> np.ndarray:
    return y_scaled * std + mean


def _compute_loglik(y_scaled: np.ndarray, mu: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    # y_scaled: (n_obs,)
    # mu: (n_draws, n_obs)
    # sigma: (n_draws,)
    sigma = np.maximum(sigma, 1e-12)
    resid = (y_scaled[None, :] - mu) / sigma[:, None]
    return -0.5 * np.log(2 * math.pi) - np.log(sigma)[:, None] - 0.5 * resid**2


def _stack_posterior(idata: az.InferenceData, var: str) -> xr.DataArray:
    da = idata.posterior[var]
    return da.stack(sample=("chain", "draw"))


def _build_ppc_figure(
    y_obs: np.ndarray,
    y_rep: np.ndarray,
    *,
    title: str,
    x_label: str,
    out_path: Path,
) -> None:
    # y_rep shape: (n_draws, n_obs)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.5, 5.2))

    # Posterior predictive overlays (thin, low alpha)
    bins = 18
    for i in np.linspace(0, y_rep.shape[0] - 1, num=min(80, y_rep.shape[0]), dtype=int):
        ax.hist(
            y_rep[i],
            bins=bins,
            density=True,
            histtype="step",
            linewidth=0.8,
            alpha=0.15,
            color="#1f77b4",
        )

    # Observed distribution
    ax.hist(
        y_obs,
        bins=bins,
        density=True,
        histtype="step",
        linewidth=2.2,
        color="black",
        label="Observed",
    )

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Density")
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    data = _read_data()

    from nplsgwqi.preprocessing import Preprocessor

    prep = Preprocessor(
        main_compositional_vars=["Ca", "Mg", "Na", "K", "Cl", "SO4", "HCO3", "NO3", "F"],
        non_compositional_vars=["pH", "EC", "TDS"],
    )

    specs = [
        {
            "endpoint": "NPLS_GWQI",
            "pkl": MODEL_STORE_DIR / "NPLS_GWQI_model.pkl",
            "y": _load_npls_gwqi_target().reindex(data.index),
            "x": prep.get_augmented_predictor_block(data.reset_index(drop=True), endpoint_name="NPLS_GWQI"),
            "xlabel": "NPLS-GWQI",
        },
        {
            "endpoint": "F",
            "pkl": MODEL_STORE_DIR / "F_model.pkl",
            "y": data["F"].astype(float),
            "x": prep.get_augmented_predictor_block(data.reset_index(drop=True), endpoint_name="F"),
            "xlabel": "F (mg/L)",
        },
    ]

    diag_rows: list[dict[str, object]] = []

    for spec in specs:
        obj = pickle.load(open(spec["pkl"], "rb"))
        idata: az.InferenceData = obj["trace"]

        x_df: pd.DataFrame = spec["x"].copy()
        expected_cols = list(obj.get("predictor_columns", []))
        if expected_cols:
            missing = [c for c in expected_cols if c not in x_df.columns]
            extra = [c for c in x_df.columns if c not in expected_cols]
            if missing:
                raise ValueError(f"{spec['endpoint']}: augmented block missing expected columns: {missing}")
            # Align exactly to the fitted model's feature order.
            x_df = x_df.loc[:, expected_cols]
        y: pd.Series = spec["y"].astype(float).copy()

        X_scaled, _, _ = _standardize(x_df)
        y_scaled, y_mean, y_std = _standardize_y(y)

        beta = _stack_posterior(idata, "beta").transpose("sample", ...).values  # (n_draws, n_features)
        alpha = _stack_posterior(idata, "alpha").values  # (n_draws,)
        sigma = _stack_posterior(idata, "sigma").values  # (n_draws,)

        linear = X_scaled @ beta.T  # (n_obs, n_draws)
        mu = (linear + alpha[None, :]).T  # (n_draws, n_obs)

        rng = np.random.default_rng(20260417)
        y_rep_scaled = rng.normal(loc=mu, scale=np.maximum(sigma, 1e-12)[:, None])
        y_rep = _inverse_y(y_rep_scaled, y_mean, y_std)

        # PPC figure (fixes blank-output issue in prior run)
        ppc_path = ARTIFACTS_DIR / f"FigureS10_PPC_{spec['endpoint']}.png"
        _build_ppc_figure(
            y_obs=y.to_numpy(dtype=float),
            y_rep=y_rep,
            title=f"PPC (distribution overlay): {spec['endpoint']}",
            x_label=spec["xlabel"],
            out_path=ppc_path,
        )

        # Diagnostics from existing posterior
        summary = az.summary(idata, var_names=["alpha", "sigma", "tau", "beta"], round_to=None)
        rhat_max = float(summary["r_hat"].max())
        ess_bulk_min = float(summary["ess_bulk"].min())
        ess_tail_min = float(summary["ess_tail"].min())

        # Divergences / treedepth
        div = int(idata.sample_stats.get("diverging", xr.DataArray([0])).values.sum())
        td = idata.sample_stats.get("tree_depth", None)
        max_td = int(td.values.max()) if td is not None else None
        reached = idata.sample_stats.get("reached_max_treedepth", None)
        n_reached = int(reached.values.sum()) if reached is not None else None

        # Add pointwise log-likelihood (Normal in y_scaled) to compute LOO/WAIC reproducibly
        ll = _compute_loglik(y_scaled=y_scaled, mu=mu, sigma=sigma)  # (n_draws, n_obs)
        # reshape to chain/draw dims for ArviZ
        chain = idata.posterior.dims.get("chain", 1)
        draw = idata.posterior.dims.get("draw", ll.shape[0] // chain)
        ll = ll.reshape((chain, draw, ll.shape[1]))
        ll_da = xr.DataArray(
            ll,
            dims=("chain", "draw", "obs_id"),
            coords={"chain": idata.posterior["chain"], "draw": idata.posterior["draw"], "obs_id": np.arange(ll.shape[2])},
            name="Y_obs",
        )
        idata_ll = idata.copy()
        idata_ll.add_groups({"log_likelihood": {"Y_obs": ll_da}})

        loo = az.loo(idata_ll, pointwise=True)
        waic = az.waic(idata_ll, pointwise=True)
        pareto_k = loo.pareto_k.values

        diag_rows.append(
            {
                "endpoint": spec["endpoint"],
                "chains": int(chain),
                "draws_per_chain": int(draw),
                "rhat_max": rhat_max,
                "ess_bulk_min": ess_bulk_min,
                "ess_tail_min": ess_tail_min,
                "divergences": div,
                "max_tree_depth": max_td,
                "reached_max_treedepth": n_reached,
                "loo_elpd": float(loo.elpd_loo),
                "loo_p": float(loo.p_loo),
                "pareto_k_max": float(np.nanmax(pareto_k)),
                "pareto_k_gt_0_7": int(np.sum(pareto_k > 0.7)),
                "pareto_k_gt_1_0": int(np.sum(pareto_k > 1.0)),
                "waic_elpd": float(waic.elpd_waic),
                "waic_p": float(waic.p_waic),
            }
        )

    diag_df = pd.DataFrame(diag_rows)
    out_csv = ARTIFACTS_DIR / "Table19_Bayes_Diagnostics_DriverModels.csv"
    diag_df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
