from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from phreeqpython import PhreeqPython

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

from nplsgwqi.bayesian_regression import run_bayesian_endpoint_model
from nplsgwqi.preprocessing import Preprocessor
from sklearn.decomposition import PCA

from nplsgwqi.compositional import ilr_transform
from nplsgwqi.process_discovery import helmert_basis
from nplsgwqi.irrigation import calculate_irrigation_indices, assess_irrigation_suitability


MANUSCRIPT_DIR = ROOT / "manuscript"
ARTIFACTS_DIR = MANUSCRIPT_DIR / "artifacts"


def _read_data() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "data.csv")
    if "SampleID" in df.columns:
        df = df.set_index("SampleID", drop=False)
    return df


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _write_plot(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def generate_cbe_tables() -> dict[str, Path]:
    df = _read_data()
    prep = Preprocessor(
        main_compositional_vars=["Ca", "Mg", "Na", "K", "Cl", "SO4", "HCO3", "NO3", "F"],
        non_compositional_vars=["pH", "EC", "TDS"],
    )
    data_meq = prep.to_meq(df)
    cbe = prep.calculate_charge_balance(data_meq).astype(float)
    abs_cbe = cbe.abs()

    per_sample = pd.DataFrame(
        {
            "SampleID": df["SampleID"].astype(str),
            "CBE_pct": cbe.values,
            "Abs_CBE_pct": abs_cbe.values,
        }
    ).sort_values("Abs_CBE_pct", ascending=False)

    summary = pd.DataFrame(
        [
            {
                "Metric": "n_samples",
                "Value": int(len(cbe)),
            },
            {
                "Metric": "mean_CBE_pct",
                "Value": float(cbe.mean()),
            },
            {
                "Metric": "median_CBE_pct",
                "Value": float(cbe.median()),
            },
            {
                "Metric": "min_CBE_pct",
                "Value": float(cbe.min()),
            },
            {
                "Metric": "max_CBE_pct",
                "Value": float(cbe.max()),
            },
            {
                "Metric": "pct_within_±5",
                "Value": float((abs_cbe <= 5).mean() * 100.0),
            },
            {
                "Metric": "pct_within_±10",
                "Value": float((abs_cbe <= 10).mean() * 100.0),
            },
        ]
    )

    p1 = ARTIFACTS_DIR / "Table16_CBE_PerSample.csv"
    p2 = ARTIFACTS_DIR / "Table17_CBE_Summary.csv"
    _write_csv(per_sample, p1)
    _write_csv(summary, p2)
    return {"per_sample": p1, "summary": p2}


def _phreeqc_solution_composition(row: pd.Series, temp_c: float) -> dict[str, str]:
    def as_float(x) -> float:
        try:
            return float(x)
        except Exception:
            return float("nan")

    # Use PHREEQC "as" syntax for species that are commonly reported as compounds.
    # Values are supplied in mg/L and interpreted under "-units mg/L".
    return {
        "-units": "mg/L",
        "-temp": f"{temp_c}",
        "pH": f"{as_float(row['pH']):.6g}",
        "Na": f"{as_float(row['Na']):.6g}",
        "K": f"{as_float(row['K']):.6g}",
        "Ca": f"{as_float(row['Ca']):.6g}",
        "Mg": f"{as_float(row['Mg']):.6g}",
        "Cl": f"{as_float(row['Cl']):.6g}",
        "S(6)": f"{as_float(row['SO4']):.6g} as SO4",
        "N(5)": f"{as_float(row['NO3']):.6g} as NO3",
        "F": f"{as_float(row['F']):.6g}",
        "Alkalinity": f"{as_float(row['HCO3']):.6g} as HCO3",
    }


def generate_saturation_indices() -> dict[str, Path]:
    df = _read_data()

    # NOTE: The bundled VIPHREEQC backend used by phreeqpython can fail to load some
    # databases. We use phreeqc.dat as the most compatible choice and report SI for
    # key carbonate/fluoride phases.
    pp = PhreeqPython(database="phreeqc.dat")
    phases = ["Calcite", "Dolomite", "Fluorite"]
    temps = [20.0, 25.0, 30.0]

    rows: list[dict[str, object]] = []
    for sample_id, row in df.iterrows():
        for temp in temps:
            sol = pp.add_solution(_phreeqc_solution_composition(row, temp))
            out: dict[str, object] = {"SampleID": str(sample_id), "Temp_C": float(temp)}
            for phase in phases:
                try:
                    out[f"SI_{phase}"] = float(sol.si(phase))
                except Exception:
                    out[f"SI_{phase}"] = float("nan")
            rows.append(out)
            sol.forget()

    si_df = pd.DataFrame(rows).sort_values(["SampleID", "Temp_C"])
    out_csv = ARTIFACTS_DIR / "Table18_Saturation_Indices.csv"
    _write_csv(si_df, out_csv)

    # Simple summary plot: fluorite SI vs fluoride at 25C + histogram
    si_25 = si_df.loc[si_df["Temp_C"] == 25.0].copy()
    si_25 = si_25.merge(df[["SampleID", "F"]].reset_index(drop=True), on="SampleID", how="left")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].scatter(si_25["F"], si_25["SI_Fluorite"], alpha=0.85)
    axes[0].axhline(0.0, color="black", linewidth=1)
    axes[0].set_xlabel("F (mg/L)")
    axes[0].set_ylabel("SI(Fluorite) at 25°C")
    axes[0].set_title("Fluoride vs Fluorite SI")

    vals = si_25["SI_Fluorite"].dropna().to_numpy()
    axes[1].hist(vals, bins=12, edgecolor="black", alpha=0.9)
    axes[1].axvline(0.0, color="black", linewidth=1)
    axes[1].set_xlabel("SI(Fluorite) at 25°C")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Distribution")

    fig_path = ARTIFACTS_DIR / "FigureS9_Saturation_Indices.png"
    _write_plot(fig, fig_path)
    return {"table": out_csv, "figure": fig_path}


@dataclass(frozen=True)
class DriverModelSpec:
    endpoint: str
    target_col: str
    use_target_from_wqi_table: bool = False


def _load_npls_gwqi_target() -> pd.Series:
    # Use the already-reported NPLS-GWQI values used elsewhere in manuscript artifacts.
    wqi = pd.read_csv(ARTIFACTS_DIR / "TableS1_Full_WQI.csv")
    wqi = wqi.set_index("SampleID")
    return wqi["NPLS_GWQI"].astype(float)


def _fit_driver_model(
    data: pd.DataFrame,
    spec: DriverModelSpec,
    *,
    draws: int,
    tune: int,
    chains: int,
    random_state: int,
    tau_beta: float,
    target_accept: float,
    max_treedepth: int | None = None,
) -> dict:
    prep = Preprocessor(
        main_compositional_vars=["Ca", "Mg", "Na", "K", "Cl", "SO4", "HCO3", "NO3", "F"],
        non_compositional_vars=["pH", "EC", "TDS"],
    )
    X_aug = prep.get_augmented_predictor_block(data.reset_index(drop=True), endpoint_name=spec.endpoint)
    X_aug.index = data.index

    if spec.use_target_from_wqi_table:
        y = _load_npls_gwqi_target().reindex(data.index)
    else:
        y = data[spec.target_col].astype(float)

    return run_bayesian_endpoint_model(
        X_aug,
        y,
        spec.endpoint,
        use_augmented=True,
        draws=draws,
        tune=tune,
        chains=chains,
        random_state=random_state,
        log_likelihood=True,
        posterior_predictive=True,
        tau_beta=tau_beta,
        lam_beta=1.0,
        target_accept=target_accept,
        max_treedepth=max_treedepth,
    )


def generate_driver_model_diagnostics() -> dict[str, Path]:
    data = _read_data()

    specs = [
        DriverModelSpec(endpoint="NPLS_GWQI", target_col="NPLS_GWQI", use_target_from_wqi_table=True),
        DriverModelSpec(endpoint="F", target_col="F", use_target_from_wqi_table=False),
    ]

    diag_rows: list[dict[str, object]] = []
    ppc_paths: list[Path] = []

    for spec in specs:
        draws = 1000
        tune = 1000
        chains = 4
        target_accept = 0.995
        max_treedepth = 15
        fit = _fit_driver_model(
            data,
            spec,
            draws=draws,
            tune=tune,
            chains=chains,
            random_state=42,
            tau_beta=1.0,
            target_accept=target_accept,
            max_treedepth=max_treedepth,
        )
        idata = fit["trace"]

        summary = az.summary(idata, var_names=["alpha", "sigma", "tau", "beta"], round_to=None)
        rhat_max = float(summary["r_hat"].max())
        ess_bulk_min = float(summary["ess_bulk"].min())
        ess_tail_min = float(summary["ess_tail"].min())

        loo = az.loo(idata, pointwise=True)
        waic = az.waic(idata, pointwise=True)
        pareto_k = loo.pareto_k
        pareto_k_max = float(np.nanmax(pareto_k.values))

        diag_rows.append(
            {
                "endpoint": spec.endpoint,
                "draws": draws,
                "tune": tune,
                "chains": chains,
                "target_accept": target_accept,
                "max_treedepth": max_treedepth,
                "rhat_max": rhat_max,
                "ess_bulk_min": ess_bulk_min,
                "ess_tail_min": ess_tail_min,
                "loo_elpd": float(loo.elpd_loo),
                "loo_p": float(loo.p_loo),
                "pareto_k_max": pareto_k_max,
                "pareto_k_gt_0_7": int(np.sum(pareto_k.values > 0.7)),
                "waic_elpd": float(waic.elpd_waic),
                "waic_p": float(waic.p_waic),
                "n_selected_pd_gt_95": int(len(fit["selected_features"])),
            }
        )

        # PPC plot
        fig = plt.figure(figsize=(7, 3.6))
        az.plot_ppc(idata, data_pairs={"Y_obs": "Y_obs"}, num_pp_samples=100)
        plt.title(f"PPC: {spec.endpoint}")
        out_fig = ARTIFACTS_DIR / f"FigureS10_PPC_{spec.endpoint}.png"
        _write_plot(fig, out_fig)
        ppc_paths.append(out_fig)

    diag_df = pd.DataFrame(diag_rows)
    diag_path = ARTIFACTS_DIR / "Table19_Bayes_Diagnostics_DriverModels.csv"
    _write_csv(diag_df, diag_path)

    # Prior sensitivity (global scale tau beta): lighter runs
    base = diag_df.set_index("endpoint")
    sens_rows: list[dict[str, object]] = []
    for spec in specs:
        baseline_fit = _fit_driver_model(
            data,
            spec,
            draws=400,
            tune=400,
            chains=2,
            random_state=42,
            tau_beta=1.0,
            target_accept=0.95,
            max_treedepth=12,
        )
        base_score = baseline_fit["driver_score"]
        for tau_beta in [0.5, 2.0]:
            alt_fit = _fit_driver_model(
                data,
                spec,
                draws=400,
                tune=400,
                chains=2,
                random_state=42,
                tau_beta=tau_beta,
                target_accept=0.95,
                max_treedepth=12,
            )
            alt_score = alt_fit["driver_score"]
            # Spearman rank correlation on shared features
            common = base_score.index.intersection(alt_score.index)
            base_rank = base_score.loc[common].rank(ascending=False)
            alt_rank = alt_score.loc[common].rank(ascending=False)
            rho = float(base_rank.corr(alt_rank, method="spearman"))
            sens_rows.append(
                {
                    "endpoint": spec.endpoint,
                    "tau_beta": float(tau_beta),
                    "rank_spearman_rho": rho,
                    "n_selected_pd_gt_95": int(len(alt_fit["selected_features"])),
                }
            )

    sens_df = pd.DataFrame(sens_rows)
    sens_path = ARTIFACTS_DIR / "Table20_Prior_Sensitivity_Horseshoe.csv"
    _write_csv(sens_df, sens_path)

    return {"diagnostics": diag_path, "ppc_figures": ppc_paths[0], "sensitivity": sens_path}


def generate_ilr_basis_sensitivity() -> dict[str, Path]:
    data = _read_data()
    comp_cols = ["Ca", "Mg", "Na", "K", "Cl", "SO4", "HCO3", "NO3", "F"]
    comp = data[comp_cols].copy()
    # Replace nonpositive with a small floor per column to keep logs stable.
    for col in comp.columns:
        pos = comp.loc[comp[col] > 0, col]
        floor = float(pos.min() / math.sqrt(2)) if len(pos) else 1e-6
        comp.loc[comp[col] <= 0, col] = floor

    n_components = 5

    V = helmert_basis(len(comp_cols))
    ilr_base = ilr_transform(comp, V=V)
    pca_base = PCA(n_components=n_components).fit(ilr_base)
    base_load = pd.DataFrame(
        pca_base.components_ @ V.T,
        columns=comp.columns,
        index=[f"Process_{i+1}" for i in range(n_components)],
    )

    rng = np.random.default_rng(12345)

    rows: list[dict[str, object]] = []
    for k in range(1, 3):
        A = rng.normal(size=(len(comp_cols) - 1, len(comp_cols) - 1))
        Q, _ = np.linalg.qr(A)
        V_alt = V @ Q
        ilr_alt = ilr_transform(comp, V=V_alt)
        pca_alt = PCA(n_components=n_components).fit(ilr_alt)
        alt_load = pd.DataFrame(
            pca_alt.components_ @ V_alt.T,
            columns=comp.columns,
            index=[f"Process_{i+1}" for i in range(n_components)],
        )

        # Compare absolute CLR loading patterns (component order may permute/sign-flip)
        base_abs = base_load.abs().to_numpy()
        alt_abs = alt_load.abs().to_numpy()
        corr = np.corrcoef(base_abs.reshape(n_components, -1), alt_abs.reshape(n_components, -1))[:n_components, n_components:]
        best = corr.max(axis=1)
        rows.append(
            {
                "basis_variant": f"rotation_{k}",
                "mean_abs_corr_best_match": float(np.mean(best)),
                "min_abs_corr_best_match": float(np.min(best)),
            }
        )

    out = pd.DataFrame(rows)
    out_path = ARTIFACTS_DIR / "Table21_ILR_Basis_Sensitivity.csv"
    _write_csv(out, out_path)
    return {"table": out_path}


def generate_wqi_benchmarking() -> dict[str, Path]:
    data = _read_data()

    standards = {
        "EC": 1500.0,
        "TDS": 1000.0,
        "Na": 200.0,
        "K": 12.0,
        "Mg": 150.0,
        "Ca": 200.0,
        "Cl": 250.0,
        "SO4": 250.0,
        "HCO3": 500.0,
        "NO3": 50.0,
        "F": 1.5,
        # pH handled separately
    }

    cols = ["pH"] + list(standards.keys())
    df = data[["SampleID"] + cols].copy()

    # Hazard ratios Y (lower better), with bipolar pH hazard
    Y = pd.DataFrame(index=df.index)
    Y[list(standards.keys())] = df[list(standards.keys())].div(pd.Series(standards), axis=1)
    ph = df["pH"].astype(float)
    acid_haz = (7.0 - ph) / 0.5
    alk_haz = (ph - 7.0) / 1.5
    Y["pH"] = np.maximum(acid_haz, alk_haz).clip(lower=0.0)

    features = Y.columns.tolist()
    p = len(features)
    w_equal = np.ones(p) / p

    # Conventional arithmetic WQI
    wqi_conv = (Y * 100.0).dot(w_equal)

    # Simple fuzzy-only compliance score (no indeterminacy/falsity): T = exp(-1.5 * Y)
    T = np.exp(-1.5 * Y)
    wqi_fuzzy = (1.0 - T.dot(w_equal)) * 100.0

    # NPLS-GWQI from existing artifact table
    npls = _load_npls_gwqi_target().reindex(df.index)

    def cls(x: pd.Series) -> pd.Series:
        return pd.cut(
            x,
            bins=[0, 50, 100, 200, 300, np.inf],
            labels=["Excellent", "Good", "Poor", "Very Poor", "Unsuitable"],
            include_lowest=True,
        ).astype(str)

    out = pd.DataFrame(
        {
            "SampleID": df["SampleID"].astype(str),
            "WQI_Conventional": wqi_conv.values,
            "WQI_Fuzzy": wqi_fuzzy.values,
            "NPLS_GWQI": npls.values,
            "Class_Conventional": cls(wqi_conv).values,
            "Class_Fuzzy": cls(wqi_fuzzy).values,
            "Class_NPLS": cls(npls).values,
        }
    )

    # Monte Carlo weight perturbation (Dirichlet around equal weights)
    rng = np.random.default_rng(20260417)
    n_mc = 5000
    alpha = np.ones(p) * 20.0
    W = rng.dirichlet(alpha, size=n_mc)
    scores = (Y.to_numpy() * 100.0) @ W.T  # (n_samples, n_mc)
    score_df = pd.DataFrame(scores, index=df.index)
    # class probability per sample
    class_bins = [0, 50, 100, 200, 300, np.inf]
    class_labels = ["Excellent", "Good", "Poor", "Very Poor", "Unsuitable"]
    probs_rows = []
    for idx, arr in score_df.iterrows():
        counts = pd.cut(arr.to_numpy(), bins=class_bins, labels=class_labels, include_lowest=True).value_counts()
        probs = {lab: float(counts.get(lab, 0) / n_mc * 100.0) for lab in class_labels}
        most = max(probs, key=probs.get)
        probs_rows.append(
            {
                "SampleID": data.loc[idx, "SampleID"],
                **probs,
                "Most_Probable_Class": most,
                "Confidence_Pct": float(probs[most]),
            }
        )
    probs_df = pd.DataFrame(probs_rows)

    out_path = ARTIFACTS_DIR / "Table22_WQI_Benchmarking.csv"
    probs_path = ARTIFACTS_DIR / "Table23_ConventionalWQI_WeightUncertainty.csv"
    _write_csv(out, out_path)
    _write_csv(probs_df, probs_path)

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.scatter(out["WQI_Conventional"], out["NPLS_GWQI"], alpha=0.8)
    ax.set_xlabel("Conventional WQI (equal weights)")
    ax.set_ylabel("NPLS-GWQI")
    ax.set_title("WQI benchmarking (sample-level)")
    fig_path = ARTIFACTS_DIR / "FigureS11_WQI_Benchmarking.png"
    _write_plot(fig, fig_path)

    return {"benchmark": out_path, "weight_uncertainty": probs_path, "figure": fig_path}


def generate_irrigation_disagreement() -> dict[str, Path]:
    data = _read_data()
    prep = Preprocessor(
        main_compositional_vars=["Ca", "Mg", "Na", "K", "Cl", "SO4", "HCO3", "NO3", "F"],
        non_compositional_vars=["pH", "EC", "TDS"],
    )
    data_meq = prep.to_meq(data)
    indices = calculate_irrigation_indices(data_meq)
    nis = assess_irrigation_suitability(indices)

    # Map classical classes onto the 4-class N-ISI scale for comparison.
    sar_map = {"Excellent": "Excellent", "Good": "Good", "Doubtful": "Marginal", "Unsuitable": "Unsuitable"}
    ssp_map = {"Excellent": "Excellent", "Good": "Good", "Permissible": "Good", "Doubtful": "Marginal", "Unsuitable": "Unsuitable"}
    rsc_map = {"Safe": "Excellent", "Marginal": "Marginal", "Unsuitable": "Unsuitable"}

    per_sample = pd.DataFrame(
        {
            "SampleID": data["SampleID"].astype(str),
            "N_ISI_Score": nis["N_ISI_Score"].astype(float).values,
            "N_ISI_Class": nis["ISI_Class"].astype(str).values,
            "SAR": indices["SAR"].astype(float).values,
            "SAR_Class": indices["SAR_Class"].astype(str).values,
            "SSP": indices["SSP"].astype(float).values,
            "SSP_Class": indices["SSP_Class"].astype(str).values,
            "RSC": indices["RSC"].astype(float).values,
            "RSC_Class": indices["RSC_Class"].astype(str).values,
        }
    )
    per_sample["SAR_Class_4"] = per_sample["SAR_Class"].map(sar_map).fillna("Unknown")
    per_sample["SSP_Class_4"] = per_sample["SSP_Class"].map(ssp_map).fillna("Unknown")
    per_sample["RSC_Class_4"] = per_sample["RSC_Class"].map(rsc_map).fillna("Unknown")

    # Summary disagreement counts
    def disagree(col: str) -> pd.Series:
        return (per_sample[col] != per_sample["N_ISI_Class"]).astype(int)

    summary = pd.DataFrame(
        [
            {"Criterion": "SAR", "Disagree_Count": int(disagree("SAR_Class_4").sum()), "Disagree_Pct": float(disagree("SAR_Class_4").mean() * 100.0)},
            {"Criterion": "SSP", "Disagree_Count": int(disagree("SSP_Class_4").sum()), "Disagree_Pct": float(disagree("SSP_Class_4").mean() * 100.0)},
            {"Criterion": "RSC", "Disagree_Count": int(disagree("RSC_Class_4").sum()), "Disagree_Pct": float(disagree("RSC_Class_4").mean() * 100.0)},
        ]
    )

    per_path = ARTIFACTS_DIR / "Table24_Irrigation_PerSample.csv"
    sum_path = ARTIFACTS_DIR / "Table25_Irrigation_Disagreement.csv"
    _write_csv(per_sample, per_path)
    _write_csv(summary, sum_path)
    return {"per_sample": per_path, "summary": sum_path}


def generate_risk_uncertainty_decomposition() -> dict[str, Path]:
    idata = az.from_netcdf(MANUSCRIPT_DIR / "model_store" / "hbmpra_trace.nc")
    # Focus on HI_skeletal_dental (fluoride) and HI_hemato (nitrate)
    # Posterior arrays shape: (chain, draw, group)
    hi_sd = idata.posterior["HI_skeletal_dental"].values
    hi_he = idata.posterior["HI_hemato"].values
    bw = idata.posterior["BW_g"].values
    ir = idata.posterior["IR_perkg_g"].values
    c_f = idata.posterior["C_F"].values
    c_no3 = idata.posterior["C_NO3"].values

    groups = ["Adults", "Children", "Teens"]

    def var_share(log_terms: dict[str, np.ndarray]) -> dict[str, float]:
        # Approximate share via variance of log term (ignoring covariance) for interpretability.
        vars_ = {k: float(np.nanvar(v)) for k, v in log_terms.items()}
        denom = sum(vars_.values()) or 1.0
        return {k: (vv / denom) * 100.0 for k, vv in vars_.items()}

    rows: list[dict[str, object]] = []
    thr_rows: list[dict[str, object]] = []
    thresholds = [0.5, 1.0, 2.0]
    for g_idx, g in enumerate(groups):
        hi_sd_g = hi_sd[:, :, g_idx].reshape(-1)
        hi_he_g = hi_he[:, :, g_idx].reshape(-1)

        # Decompose log HI into log C + log IR_perkg - (log BW cancels in IR_perkg definition)
        # Here we use IR_perkg directly; BW still appears in dermal terms in other contaminants, but F/NO3 ingestion-only.
        c_f_g = c_f.reshape(-1) if c_f.ndim == 2 else c_f[:, :, g_idx].reshape(-1)
        c_no3_g = c_no3.reshape(-1) if c_no3.ndim == 2 else c_no3[:, :, g_idx].reshape(-1)
        ir_g = ir[:, :, g_idx].reshape(-1)

        shares_sd = var_share({"log_C_F": np.log(np.maximum(c_f_g, 1e-12)), "log_IR_perkg": np.log(np.maximum(ir_g, 1e-12))})
        shares_he = var_share({"log_C_NO3": np.log(np.maximum(c_no3_g, 1e-12)), "log_IR_perkg": np.log(np.maximum(ir_g, 1e-12))})

        rows.append(
            {
                "Group": g,
                "HI_SkeletalDental_Mean": float(np.mean(hi_sd_g)),
                "P(HI_SD>1)": float(np.mean(hi_sd_g > 1.0)),
                "VarShare_SD_log_C_pct": float(shares_sd["log_C_F"]),
                "VarShare_SD_log_IR_perkg_pct": float(shares_sd["log_IR_perkg"]),
                "HI_Hemato_Mean": float(np.mean(hi_he_g)),
                "P(HI_Hemato>1)": float(np.mean(hi_he_g > 1.0)),
                "VarShare_Hemato_log_C_pct": float(shares_he["log_C_NO3"]),
                "VarShare_Hemato_log_IR_perkg_pct": float(shares_he["log_IR_perkg"]),
            }
        )

        thr_row: dict[str, object] = {"Group": g}
        for t in thresholds:
            thr_row[f"P(HI_SD>{t})"] = float(np.mean(hi_sd_g > t))
            thr_row[f"P(HI_Hemato>{t})"] = float(np.mean(hi_he_g > t))
        thr_rows.append(thr_row)

    out = pd.DataFrame(rows)
    out_path = ARTIFACTS_DIR / "Table26_Risk_Uncertainty_Decomposition.csv"
    _write_csv(out, out_path)
    thr_df = pd.DataFrame(thr_rows)
    thr_path = ARTIFACTS_DIR / "Table27_Risk_Threshold_Probabilities.csv"
    _write_csv(thr_df, thr_path)
    return {"table": out_path, "thresholds": thr_path}


def main() -> None:
    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

    print("Generating CBE tables…")
    generate_cbe_tables()

    print("Generating PHREEQC saturation indices…")
    generate_saturation_indices()

    print("Generating Bayesian driver model diagnostics + PPC + prior sensitivity…")
    generate_driver_model_diagnostics()

    print("Generating ILR basis sensitivity summary…")
    generate_ilr_basis_sensitivity()

    print("Generating WQI benchmarking + weight uncertainty…")
    generate_wqi_benchmarking()

    print("Generating irrigation disagreement tables…")
    generate_irrigation_disagreement()

    print("Generating risk uncertainty decomposition…")
    generate_risk_uncertainty_decomposition()

    print("OK: reviewer revision artifacts generated.")


if __name__ == "__main__":
    main()
