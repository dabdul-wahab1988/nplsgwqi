from __future__ import annotations

import math
import pickle
from pathlib import Path
import subprocess
import sys

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from sklearn.preprocessing import StandardScaler

try:
    from scipy.special import logsumexp
except Exception:  # pragma: no cover
    logsumexp = None

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


def _logmeanexp(x: np.ndarray) -> float:
    if logsumexp is None:
        m = float(np.max(x))
        return m + float(np.log(np.mean(np.exp(x - m))))
    return float(logsumexp(x) - np.log(x.size))


def _logpdf_normal(y: float, mu: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    sigma = np.maximum(sigma, 1e-12)
    resid = (y - mu) / sigma
    return -0.5 * np.log(2 * math.pi) - np.log(sigma) - 0.5 * resid**2


def _logpdf_studentt(y: float, mu: np.ndarray, sigma: np.ndarray, nu: float) -> np.ndarray:
    sigma = np.maximum(sigma, 1e-12)
    z = (y - mu) / sigma
    try:
        from scipy.stats import t as student_t  # type: ignore

        return student_t.logpdf(z, df=float(nu)) - np.log(sigma)
    except Exception:
        # Fallback closed-form; relies on scipy.special.gammaln.
        from scipy.special import gammaln  # type: ignore

        nu = float(nu)
        return (
            gammaln((nu + 1) / 2)
            - gammaln(nu / 2)
            - 0.5 * np.log(nu * math.pi)
            - np.log(sigma)
            - ((nu + 1) / 2) * np.log1p((z**2) / nu)
        )


def _exact_lpd_leftout(
    *,
    endpoint: str,
    x_df: pd.DataFrame,
    y: pd.Series,
    leftout_i: int,
    fit_config: dict[str, object],
    random_seed: int,
) -> float:
    from nplsgwqi.bayesian_regression import run_bayesian_endpoint_model

    # Drop by *position*, not by label (SampleID can repeat in this dataset).
    x_train = x_df.drop(index=x_df.index[leftout_i]).reset_index(drop=True)
    y_arr = y.to_numpy(dtype=float)
    y_train = pd.Series(np.delete(y_arr, leftout_i))

    # Exact-LOO (reloo-style) refits can be expensive; we use a lighter refit that is still
    # adequate for pointwise predictive density replacement on small n (n=34).
    reloo_draws = int(min(int(fit_config.get("draws", 1500)), 400))
    reloo_tune = int(min(int(fit_config.get("tune", 1500)), 400))

    res = run_bayesian_endpoint_model(
        x_train,
        y_train,
        endpoint,
        include_scale=bool(fit_config.get("include_scale", True)),
        use_augmented=bool(fit_config.get("use_augmented", True)),
        draws=reloo_draws,
        tune=reloo_tune,
        chains=1,
        cores=1,
        random_state=int(random_seed),
        log_likelihood=bool(fit_config.get("log_likelihood", True)),
        posterior_predictive=False,
        tau_beta=float(fit_config.get("tau_beta", 0.5)),
        lam_beta=float(fit_config.get("lam_beta", 1.0)),
        target_accept=float(fit_config.get("target_accept", 0.99)),
        max_treedepth=int(fit_config["max_treedepth"]) if fit_config.get("max_treedepth", None) is not None else None,
        likelihood=str(fit_config.get("likelihood", "studentt")),
        studentt_nu=float(fit_config.get("studentt_nu", 4.0)),
        regularized_horseshoe=bool(fit_config.get("regularized_horseshoe", True)),
        slab_scale=float(fit_config.get("slab_scale", 2.0)),
        shrinkage_family=str(fit_config.get("shrinkage_family", "normal")),
    )
    idata: az.InferenceData = res["trace"]

    x_scaler = StandardScaler()
    x_scaler.fit(x_train.to_numpy(dtype=float))
    x_left = x_scaler.transform(x_df.iloc[[leftout_i]].to_numpy(dtype=float)).reshape(-1)

    y_scaler = StandardScaler()
    y_scaler.fit(y_train.to_numpy(dtype=float).reshape(-1, 1))
    y_left = float(y_scaler.transform(np.array([[float(y_arr[leftout_i])]])).ravel()[0])

    beta = _stack_posterior(idata, "beta").transpose("sample", ...).values  # (n_draws, n_features)
    alpha = _stack_posterior(idata, "alpha").values  # (n_draws,)
    sigma = _stack_posterior(idata, "sigma").values  # (n_draws,)
    mu = alpha + beta.dot(x_left)

    like = str(fit_config.get("likelihood", "studentt")).lower()
    if like == "normal":
        logpdf = _logpdf_normal(y_left, mu, sigma)
    else:
        logpdf = _logpdf_studentt(y_left, mu, sigma, nu=float(fit_config.get("studentt_nu", 4.0)))
    return _logmeanexp(logpdf)


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

    # Use a clean, journal-friendly PPC: predictive density band + observed line.
    n_bins = 24
    lo, hi = np.nanmin(np.r_[y_obs, y_rep.ravel()]), np.nanmax(np.r_[y_obs, y_rep.ravel()])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = 0.0, 1.0
    pad = 0.03 * (hi - lo)
    edges = np.linspace(lo - pad, hi + pad, num=n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Predictive histogram density per draw
    dens = []
    for i in range(y_rep.shape[0]):
        h, _ = np.histogram(y_rep[i], bins=edges, density=True)
        dens.append(h)
    dens = np.asarray(dens, dtype=float)  # (n_draws, n_bins)
    q05, q50, q95 = np.nanpercentile(dens, [5, 50, 95], axis=0)

    h_obs, _ = np.histogram(y_obs, bins=edges, density=True)

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    ax.fill_between(centers, q05, q95, color="#1f77b4", alpha=0.25, label="Posterior predictive (90%)")
    ax.plot(centers, q50, color="#1f77b4", linewidth=2.0, label="Posterior predictive (median)")
    ax.plot(centers, h_obs, color="black", linewidth=2.0, label="Observed")

    ax.set_title(title, fontsize=13)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.legend(frameon=False, fontsize=10, loc="best")
    fig.tight_layout()

    fig.savefig(out_path, dpi=600, bbox_inches="tight")
    if out_path.suffix.lower() == ".png":
        fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    # Optional subcommand: compute one exact-LOO pointwise LPD in a fresh process.
    # This avoids repeated PyMC model compilation in a single long-running process.
    if len(sys.argv) == 4 and sys.argv[1] == "--reloo-one":
        endpoint = str(sys.argv[2])
        leftout_i = int(sys.argv[3])

        data = _read_data()
        from nplsgwqi.preprocessing import Preprocessor

        prep = Preprocessor(
            main_compositional_vars=["Ca", "Mg", "Na", "K", "Cl", "SO4", "HCO3", "NO3", "F"],
            non_compositional_vars=["pH", "EC", "TDS"],
        )

        if endpoint == "NPLS_GWQI":
            y = _load_npls_gwqi_target().reindex(data.index)
            x = prep.get_augmented_predictor_block(data.reset_index(drop=True), endpoint_name="NPLS_GWQI")
            pkl = MODEL_STORE_DIR / "NPLS_GWQI_model.pkl"
        elif endpoint == "F":
            y = data["F"].astype(float)
            x = prep.get_augmented_predictor_block(data.reset_index(drop=True), endpoint_name="F")
            pkl = MODEL_STORE_DIR / "F_model.pkl"
        else:
            raise ValueError(f"Unknown endpoint: {endpoint}")

        obj = pickle.load(open(pkl, "rb"))
        fit_config = obj.get("fit_config", None)
        if not isinstance(fit_config, dict):
            raise ValueError(f"{endpoint}: fit_config missing from model_store; rerun generation_script.py.")

        x_df: pd.DataFrame = x.copy()
        expected_cols = list(obj.get("predictor_columns", []))
        if expected_cols:
            x_df = x_df.loc[:, expected_cols]

        seed = int(fit_config.get("random_state", 42)) + 10_000 + int(leftout_i)
        lpd = _exact_lpd_leftout(
            endpoint=endpoint,
            x_df=x_df,
            y=y.astype(float),
            leftout_i=leftout_i,
            fit_config=fit_config,
            random_seed=seed,
        )
        print(lpd, flush=True)
        # Work around occasional low-level crashes during interpreter cleanup after PyMC sampling
        # by exiting immediately once we've emitted the scalar result.
        import os

        os._exit(0)

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
            "nonneg": True,
        },
        {
            "endpoint": "F",
            "pkl": MODEL_STORE_DIR / "F_model.pkl",
            "y": data["F"].astype(float),
            "x": prep.get_augmented_predictor_block(data.reset_index(drop=True), endpoint_name="F"),
            "xlabel": "F (mg/L)",
            "nonneg": True,
        },
    ]

    diag_rows: list[dict[str, object]] = []

    for spec in specs:
        obj = pickle.load(open(spec["pkl"], "rb"))
        idata: az.InferenceData = obj["trace"]
        fit_config = obj.get("fit_config", None)
        chain = int(idata.posterior.sizes.get("chain", 1))
        draw = int(idata.posterior.sizes.get("draw", 0))

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
        like = str(fit_config.get("likelihood", "normal")).lower() if isinstance(fit_config, dict) else "normal"
        if like == "studentt":
            nu = float(fit_config.get("studentt_nu", 4.0)) if isinstance(fit_config, dict) else 4.0
            y_rep_scaled = rng.standard_t(df=nu, size=mu.shape) * np.maximum(sigma, 1e-12)[:, None] + mu
        else:
            y_rep_scaled = rng.normal(loc=mu, scale=np.maximum(sigma, 1e-12)[:, None])
        y_rep = _inverse_y(y_rep_scaled, y_mean, y_std)
        if bool(spec.get("nonneg", False)):
            y_rep = np.maximum(y_rep, 0.0)

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

        # Prefer stored log-likelihood if present (PyMC can persist it into InferenceData).
        if hasattr(idata, "log_likelihood") and ("Y_obs" in idata.log_likelihood.data_vars):
            idata_ll = idata
        else:
            # Reconstruct pointwise log-likelihood (Normal in y_scaled) to compute LOO/WAIC reproducibly.
            ll = _compute_loglik(y_scaled=y_scaled, mu=mu, sigma=sigma)  # (n_draws, n_obs)
            # reshape to chain/draw dims for ArviZ
            chain = int(idata.posterior.sizes.get("chain", 1))
            draw = int(idata.posterior.sizes.get("draw", ll.shape[0] // chain))
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

        loo_elpd_reloo = None
        # Exact refits are costly; prioritize only the "very bad" PSIS points (Pareto-k > 1.0).
        reloo_idx = np.where(pareto_k > 1.0)[0]
        reloo_n = int(reloo_idx.size)
        if reloo_n > 0:
            if not isinstance(fit_config, dict):
                raise ValueError(f"{spec['endpoint']}: fit_config missing from model_store; rerun generation_script.py.")

            # Replace PSIS-LOO for problematic points with exact LOO refits (reloo-style).
            loo_i = np.asarray(loo.loo_i.values, dtype=float).copy()
            script_path = str(Path(__file__).resolve())
            for leftout_i in reloo_idx:
                last_err = None
                for attempt in range(1, 4):
                    proc = subprocess.run(
                        [sys.executable, script_path, "--reloo-one", str(spec["endpoint"]), str(int(leftout_i))],
                        capture_output=True,
                        text=True,
                    )
                    stdout = (proc.stdout or "").strip().splitlines()
                    if stdout:
                        try:
                            loo_i[int(leftout_i)] = float(stdout[-1])
                            last_err = None
                            break
                        except ValueError as exc:
                            last_err = exc
                    else:
                        last_err = RuntimeError(f"no stdout; returncode={proc.returncode}")
                if last_err is not None:
                    raise RuntimeError(f"{spec['endpoint']} reloo-one({int(leftout_i)}): failed after 3 attempts") from last_err
            loo_elpd_reloo = float(np.sum(loo_i))

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
                "loo_elpd_reloo": loo_elpd_reloo,
                "loo_p": float(loo.p_loo),
                "pareto_k_max": float(np.nanmax(pareto_k)),
                "pareto_k_gt_0_7": int(np.sum(pareto_k > 0.7)),
                "pareto_k_gt_1_0": int(np.sum(pareto_k > 1.0)),
                "reloo_n": int(reloo_n),
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
