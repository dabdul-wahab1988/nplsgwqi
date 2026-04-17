from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = ROOT / "manuscript"
ARTIFACTS_DIR = MANUSCRIPT_DIR / "artifacts"

SRC_DIR = ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


def _load_data() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "data.csv")
    if "SampleID" in df.columns:
        df = df.set_index("SampleID", drop=False)
    return df


def _load_cbe_table() -> pd.DataFrame:
    path = ARTIFACTS_DIR / "Table16_CBE_PerSample.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run generation_script.py to create it.")
    return pd.read_csv(path)


def _truncnorm_nonneg(rng: np.random.Generator, mu: float, sigma: float, size: int) -> np.ndarray:
    sigma = float(max(sigma, 1e-9))
    try:
        from scipy.stats import truncnorm  # type: ignore

        a, b = (0.0 - mu) / sigma, np.inf
        return truncnorm.rvs(a, b, loc=mu, scale=sigma, size=size, random_state=rng)
    except Exception:
        # Fallback: clip (slightly underestimates tail mass vs truncnorm).
        return np.maximum(rng.normal(mu, sigma, size=size), 0.0)


def _mc_hbmpra_summary(df: pd.DataFrame, *, n_mc: int, seed: int) -> pd.DataFrame:
    """
    Fast Monte Carlo risk propagation that mirrors the HBMPRA priors used in `nplsgwqi.risk`
    without running NUTS. This is used ONLY for a sensitivity check to CBE filtering.
    """
    rng = np.random.default_rng(int(seed))

    groups = ["Adults", "Children", "Teens"]
    bw_mean = np.array([60.0, 15.0, 45.0], dtype=float)
    ir_mean = np.array([2.5, 1.0, 1.5], dtype=float)

    # Priors from nplsgwqi.risk.assess_health_risk
    CV_BW = 0.21
    sigma_log_bw = np.sqrt(np.log(1.0 + CV_BW**2))
    mu_log_bw = np.log(bw_mean) - 0.5 * np.log(1.0 + CV_BW**2)

    IR_perkg_med = ir_mean / bw_mean
    sigma_log_ir = 0.6
    mu_log_ir = np.log(IR_perkg_med)

    # Concentration priors (truncated normal, >= 0)
    mu_f, sd_f = float(df["F"].mean()), float(df["F"].std(ddof=0) + 1e-6)
    mu_n, sd_n = float(df["NO3"].mean()), float(df["NO3"].std(ddof=0) + 1e-6)
    c_f = _truncnorm_nonneg(rng, mu_f, sd_f, n_mc)
    c_no3 = _truncnorm_nonneg(rng, mu_n, sd_n, n_mc)

    # RfD (oral) from nplsgwqi.risk
    rfd_f = 0.06
    rfd_no3 = 1.6

    rows = []
    for gi, g in enumerate(groups):
        bw_g = np.exp(mu_log_bw[gi] + rng.normal(0.0, 1.0, size=n_mc) * sigma_log_bw)
        ir_perkg_g = np.exp(mu_log_ir[gi] + rng.normal(0.0, 1.0, size=n_mc) * sigma_log_ir)

        # Ingestion-only pathway; IR_perkg is already normalized by BW.
        hi_skeletal = (c_f * ir_perkg_g) / rfd_f
        hi_hemato = (c_no3 * ir_perkg_g) / rfd_no3

        for organ, arr in [("skeletal_dental", hi_skeletal), ("hemato", hi_hemato)]:
            rows.append(
                {
                    "Demographic_Group": g,
                    "Organ_System": organ,
                    "Mean_HI": float(np.mean(arr)),
                    "CI_2_5": float(np.percentile(arr, 2.5)),
                    "CI_97_5": float(np.percentile(arr, 97.5)),
                    "P_HI_greater_1": float(np.mean(arr > 1.0)),
                    "Method": "Monte Carlo (HBMPRA priors)",
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    data = _load_data()
    cbe = _load_cbe_table()
    if "Abs_CBE_pct" not in cbe.columns:
        raise ValueError("Table16_CBE_PerSample.csv missing Abs_CBE_pct column")

    keep_ids = set(cbe.loc[cbe["Abs_CBE_pct"] <= 10.0, "SampleID"].astype(str))
    data_keep = data.loc[data["SampleID"].astype(str).isin(keep_ids)].copy()

    n_mc = 20000
    seed = 20260417

    out_rows = []
    for scenario, df in [("All samples (n=34)", data), ("CBE abs<=10% (n=26)", data_keep)]:
        t = _mc_hbmpra_summary(df[["F", "NO3"]], n_mc=n_mc, seed=seed)

        t.insert(0, "Scenario", scenario)
        out_rows.append(t)

    out = pd.concat(out_rows, ignore_index=True)

    # Add simple deltas (CBE-filtered minus full) for Mean_HI and exceedance probability.
    wide = out.pivot_table(
        index=["Demographic_Group", "Organ_System"],
        columns="Scenario",
        values=["Mean_HI", "P_HI_greater_1"],
        aggfunc="first",
    )
    if ("Mean_HI", "All samples (n=34)") in wide.columns and ("Mean_HI", "CBE abs<=10% (n=26)") in wide.columns:
        wide[("Delta", "Mean_HI")] = wide[("Mean_HI", "CBE abs<=10% (n=26)")] - wide[("Mean_HI", "All samples (n=34)")]
        wide[("Delta", "P_HI_greater_1")] = (
            wide[("P_HI_greater_1", "CBE abs<=10% (n=26)")] - wide[("P_HI_greater_1", "All samples (n=34)")]
        )

    # Flatten columns for CSV readability.
    wide.columns = [f"{a}__{b}" for a, b in wide.columns.to_flat_index()]
    wide = wide.reset_index()

    out_path = ARTIFACTS_DIR / "Table28_Risk_CBE_Filter_Sensitivity.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
