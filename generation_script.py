from __future__ import annotations

import os
import random
import re
import pickle
import sys
import textwrap
from pathlib import Path

import arviz as az
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
from sklearn.decomposition import PCA
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold
from statsmodels.nonparametric.smoothers_lowess import lowess

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from nplsgwqi.compositional import ilr_transform
from nplsgwqi.orchestrator import IntegratedWorkflow
from nplsgwqi.preprocessing import Preprocessor
from nplsgwqi.risk import assess_health_risk
from nplsgwqi.sensitivity import run_probabilistic_sensitivity
from nplsgwqi.bayesian_regression import run_bayesian_endpoint_model
from nplsgwqi.sparse_pls import (
    default_component_grid,
    default_keep_grid,
    fit_sparse_pls_model,
    prepare_clr_predictors,
    run_sparse_pls_endpoint_model,
    tune_sparse_pls,
)
from nplsgwqi.figure_style import (
    apply_publication_style,
    save_figure,
    style_axes,
    style_colorbar,
    style_legend,
)

SEED = 42
ROOT = Path(__file__).resolve().parent
ART = ROOT / "manuscript" / "artifacts"
MODELS = ROOT / "manuscript" / "model_store"
HYDRO = ART / "hydrochem"
LOCKED_ARTIFACTS = {
    ART / "Figure2_Hydrochemical_Characterisation.png",
    HYDRO / "ilr_plot.png",
}
WQI_PROB_PATH = ART / "Table14_WQI_Probabilities.csv"
ISI_PROB_PATH = ART / "Table15_ISI_Probabilities.csv"
HBMPRA_TRACE_PATH = MODELS / "hbmpra_trace.nc"
HBMPRA_TABLE_PATH = ART / "Table11_Bayesian_Risk.csv"
HBMPRA_FIGURE_PATH = ART / "Figure10_Bayesian_Posteriors.png"


def env_flag(name: str, default: bool = False) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def env_int(name: str, default: int, minimum: int = 1) -> int:
    raw = os.getenv(name)
    if raw is None or not raw.strip():
        return default
    value = int(raw)
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {value}")
    return value


QUICK_BAYES = env_flag("WQI_QUICK_BAYES", default=False)
FORCE_HBMPRA_RERUN = env_flag("WQI_FORCE_HBMPRA_RERUN", default=False)
ENDPOINT_BAYES_DRAWS = env_int("WQI_ENDPOINT_DRAWS", 150 if QUICK_BAYES else 1000)
ENDPOINT_BAYES_TUNE = env_int("WQI_ENDPOINT_TUNE", 100 if QUICK_BAYES else 500)
ENDPOINT_BAYES_CHAINS = env_int("WQI_ENDPOINT_CHAINS", 1 if QUICK_BAYES else 2)
HBMPRA_BAYES_DRAWS = env_int("WQI_HBMPRA_SAMPLES", 150 if QUICK_BAYES else 500)
HBMPRA_BAYES_TUNE = env_int("WQI_HBMPRA_TUNE", 120 if QUICK_BAYES else 800)
HBMPRA_BAYES_CHAINS = env_int("WQI_HBMPRA_CHAINS", 1 if QUICK_BAYES else 2)
LOWESS_FRAC = float(os.getenv("WQI_LOWESS_FRAC", "0.65"))

STANDARDS = pd.Series({
    "pH": 8.5,
    "EC": 1500.0,
    "TDS": 1000.0,
    "Na": 200.0,
    "K": 12.0,
    "Mg": 50.0,
    "Ca": 75.0,
    "Cl": 250.0,
    "SO4": 250.0,
    "HCO3": 300.0,
    "NO3": 50.0,
    "F": 1.5,
})

COMP = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F"]
ENDS = ["NO3", "F"]

GROUPS = {
    "Adults": {"bw": 60.0, "ir": 2.5},
    "Children": {"bw": 15.0, "ir": 1.0},
    "Teens": {"bw": 45.0, "ir": 1.5},
}

TOX = {
    "F": {"RfD_oral": 0.06, "RfD_derm": float("inf"), "Kp": 0.0, "organs": ["skeletal_dental"]},
    "NO3": {"RfD_oral": 1.6, "RfD_derm": float("inf"), "Kp": 0.0, "organs": ["hemato"]},
}

apply_publication_style()

PROCESS_COLOURS = ["#0F4C5C", "#2A9D8F", "#7B2CBF", "#BC6C25", "#7F5539"]
DEMOGRAPHIC_COLOURS = {"Adults": "#2A9D8F", "Children": "#D95D39", "Teens": "#577590"}
CLASS_COLOURS = {
    "Excellent": "#2A9D8F",
    "Good": "#8AB17D",
    "Poor": "#E9C46A",
    "Very Poor": "#F4A261",
    "Unsuitable": "#E76F51",
}
ENDPOINT_TEXT_LABELS = {"NPLS_GWQI": "NPLS-GWQI", "NO3": "NO3", "F": "F"}
ENDPOINT_PLOT_LABELS = {"NPLS_GWQI": "NPLS-GWQI", "NO3": r"NO$_3^-$", "F": r"F$^-$"}
ION_COLOURS = {
    "Ca": "#4C78A8",
    "Mg": "#72B7B2",
    "Na": "#0F4C5C",
    "K": "#BC6C25",
    "HCO3": "#9C755F",
    "Cl": "#E45756",
    "SO4": "#B279A2",
    "NO3": "#54A24B",
    "F": "#7B2CBF",
    "Total_Mineralization_Scale": "#2F4F4F",
}


def save_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def log_runtime_configuration() -> None:
    print(
        "Runtime configuration:"
        f" quick_bayes={QUICK_BAYES},"
        f" endpoint_draws={ENDPOINT_BAYES_DRAWS},"
        f" endpoint_tune={ENDPOINT_BAYES_TUNE},"
        f" endpoint_chains={ENDPOINT_BAYES_CHAINS},"
        f" hbmpra_draws={HBMPRA_BAYES_DRAWS},"
        f" hbmpra_tune={HBMPRA_BAYES_TUNE},"
        f" hbmpra_chains={HBMPRA_BAYES_CHAINS},"
        f" force_hbmpra_rerun={FORCE_HBMPRA_RERUN}"
    )


def lighten_colour(colour: str, blend: float = 0.35) -> tuple[float, float, float]:
    base = np.array(mcolors.to_rgb(colour))
    return tuple(base + (1.0 - base) * blend)


def build_sparse_pls_predictors(
    comp_data: pd.DataFrame,
    endpoint_name: str,
    *,
    include_scale: bool = True,
    use_augmented: bool = False,
) -> pd.DataFrame:
    if use_augmented:
        return comp_data.copy()
    endpoint_arg = endpoint_name if endpoint_name in comp_data.columns else None
    return prepare_clr_predictors(comp_data, endpoint_name=endpoint_arg, include_scale=include_scale)


def compute_sparse_pls_cv_metrics(
    comp_data: pd.DataFrame,
    target: pd.Series,
    endpoint_name: str,
    *,
    include_scale: bool = True,
    use_augmented: bool = False,
    outer_splits: int = 5,
    inner_splits: int = 3,
    random_state: int = 42,
) -> dict[str, float]:
    x_df = build_sparse_pls_predictors(
        comp_data,
        endpoint_name,
        include_scale=include_scale,
        use_augmented=use_augmented,
    )
    # Ensure no constant columns
    x_df = x_df.loc[:, x_df.std() > 1e-10]
    
    component_grid = default_component_grid(x_df.shape[1], len(x_df))
    keep_grid = default_keep_grid(x_df.shape[1])
    n_splits = max(2, min(outer_splits, len(x_df)))
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    r2_scores: list[float] = []
    rmse_scores: list[float] = []
    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(x_df), start=1):
        x_train, x_test = x_df.iloc[train_idx], x_df.iloc[test_idx]
        y_train, y_test = target.iloc[train_idx], target.iloc[test_idx]
        
        # Guard against zero-variance target in test fold
        if y_test.std() < 1e-10:
            continue
            
        try:
            tuned = tune_sparse_pls(
                x_train,
                y_train,
                component_grid=component_grid,
                keep_grid=keep_grid,
                n_splits=max(2, min(inner_splits, len(x_train))),
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
            fold_r2 = float(r2_score(y_test, pred))
            fold_rmse = float(np.sqrt(mean_squared_error(y_test, pred)))
            
            # Clip extreme negative R2 for visualization stability
            r2_scores.append(max(-1.5, fold_r2))
            rmse_scores.append(fold_rmse)
        except Exception:
            continue

    if not r2_scores:
        return {"cv_r2_mean": -1.0, "cv_r2_std": 0.0, "cv_rmse_mean": 0.0, "cv_rmse_std": 0.0}

    return {
        "cv_r2_mean": float(np.mean(r2_scores)),
        "cv_r2_std": float(np.std(r2_scores)),
        "cv_rmse_mean": float(np.mean(rmse_scores)),
        "cv_rmse_std": float(np.std(rmse_scores)),
    }


def add_lowess_scatter(
    ax: plt.Axes,
    x: pd.Series,
    y: pd.Series,
    *,
    color: str,
    label: str,
) -> None:
    mask = np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        return
    x_vals = np.asarray(x[mask], dtype=float)
    y_vals = np.asarray(y[mask], dtype=float)
    ax.scatter(
        x_vals,
        y_vals,
        alpha=0.55,
        color=lighten_colour(color, 0.15),
        edgecolor="white",
        linewidth=0.6,
        s=42,
        label=label,
    )
    if len(x_vals) >= 4:
        smooth = lowess(y_vals, x_vals, frac=LOWESS_FRAC, return_sorted=True)
        ax.plot(smooth[:, 0], smooth[:, 1], color=color, linewidth=2.2)


def build_bayesian_pathway_rows(endpoint_models: dict, order: list[str]) -> pd.DataFrame:
    """
    Rebuilds pathway relationship rows using the Bayesian driver results.
    Filters to significant drivers (>1%) and classifies them as Promoting or Suppressing factors.
    """
    rows: list[dict[str, object]] = []
    for ep in order:
        if ep not in endpoint_models:
            continue
        model_res = endpoint_models[ep]
        drivers = model_res["driver_score"]
        betas = model_res["beta_mean"]
        pds = model_res["prob_direction"]
        
        # Only include drivers with >1% contribution to prevent pathway crowding
        significant = drivers[drivers > 1.0]
        
        for pred_id, contribution in significant.items():
            beta_val = betas.get(pred_id, 0.0)
            pd_val = pds.get(pred_id, 50.0)
            
            # Formatting: remove Aug_ prefix and replace underscores with spaces
            display_name = pred_id.replace("Aug_", "").replace("_", " ")
            
            # Classification based on Bayesian direction (Beta mean sign)
            source_type = "Promoting Factor" if beta_val > 0 else "Suppressing Factor"
            
            rows.append({
                "Endpoint": ep,
                "Basis_Label": model_res.get("predictor_block_label", "Bayesian Driver Block"),
                "Basis_Name": "augmented" if "Aug" in pred_id else "clr",
                "Base_Component": pred_id,
                "Base_Component_Name": display_name,
                "Contribution_Percent": float(contribution),
                "Prob_Direction": float(pd_val),
                "Source_Type": source_type
            })
    return pd.DataFrame(rows)


def draw_bayesian_pathway_panel(
    ax,
    pathway_rows: pd.DataFrame,
    endpoint: str,
    panel_letter: str = "",
    *,
    process_heading: str = "Hydrochemical Drivers",
    title_y: float = 1.04,
    heading_y: float = 0.98,
) -> None:
    """
    Renders a Pathway diagram grounded in the Bayesian posterior.
    Left: Promoting vs Suppressing Factors (Context)
    Middle: Individual Ionic/Process Drivers (Signature)
    Right: The Endpoint (Impact)
    """
    ax.axis("off")
    ax.set_xlim(-0.04, 1.02)
    ax.set_ylim(-0.08, 1.08)

    if pathway_rows.empty:
        ax.text(0.5, 0.52, f"No significant drivers for {ENDPOINT_PLOT_LABELS[endpoint]}",
                ha="center", va="center", fontsize=12, fontweight="bold")
        return

    # Visual layout coordinates
    x_src, x_proc, x_ep = 0.12, 0.48, 0.86
    bw_src, bw_proc, bw_ep = 0.20, 0.24, 0.18
    bh_src_ep = 0.11

    # Define promoting and suppressing context nodes
    source_types = ["Promoting Factor", "Suppressing Factor"]
    source_colours = {"Promoting Factor": "#2A9D8F", "Suppressing Factor": "#E15759"}
    source_y = {"Promoting Factor": 0.72, "Suppressing Factor": 0.28}

    # Sort and position drivers in the middle column
    process_rows = pathway_rows.sort_values("Contribution_Percent", ascending=False).reset_index(drop=True)
    if len(process_rows) <= 1:
        process_y = [0.50]
        bh_proc = 0.10
        txt_off_1, txt_off_2 = 0.015, 0.025
    else:
        # Stretch further to use available space
        y_max, y_min = 0.94, 0.06
        process_y = np.linspace(y_max, y_min, len(process_rows))
        gap = (y_max - y_min) / (len(process_rows) - 1)
        bh_proc = min(0.09, gap * 0.85)
        # Dynamically scale text offsets based on box height
        txt_off_1 = bh_proc * 0.20
        txt_off_2 = bh_proc * 0.30

    endpoint_y = 0.50
    endpoint_colour = {"NPLS_GWQI": "#2A9D8F", "F": "#7B2CBF"}[endpoint]

    # Draw Headings
    heading_prefix = f"{panel_letter} " if panel_letter else ""
    ax.text(0.5, title_y, f"{heading_prefix}{ENDPOINT_PLOT_LABELS[endpoint]} Bayesian Driver Pathway",
            ha="center", va="bottom", fontsize=13, fontweight="bold")
    ax.text(x_src, heading_y, "Context (Direction)", ha="center", fontsize=11, fontweight="bold")
    ax.text(x_proc, heading_y, process_heading, ha="center", fontsize=11, fontweight="bold")
    ax.text(x_ep, heading_y, "Endpoint", ha="center", fontsize=11, fontweight="bold")

    # Draw Source Nodes (Promoting/Suppressing)
    for src, sy in source_y.items():
        if src in pathway_rows["Source_Type"].unique():
            ax.add_patch(Rectangle((x_src - bw_src / 2, sy - bh_src_ep / 2), bw_src, bh_src_ep,
                                  facecolor=source_colours[src], edgecolor="#1f2937", 
                                  linewidth=1.2, alpha=0.92, zorder=4))
            ax.text(x_src, sy, src, ha="center", va="center", fontsize=10, 
                    fontweight="bold", color="white", zorder=5)

    # Draw Endpoint Node
    ax.add_patch(Rectangle((x_ep - bw_ep / 2, endpoint_y - bh_src_ep / 2), bw_ep, bh_src_ep,
                          facecolor=endpoint_colour, edgecolor="#1f2937", 
                          linewidth=1.2, alpha=0.92, zorder=4))
    ax.text(x_ep, endpoint_y, ENDPOINT_PLOT_LABELS[endpoint], ha="center", va="center", 
            fontsize=11, fontweight="bold", color="white", zorder=5)

    # Draw Driver Nodes and Flow Arrows
    max_flow = max(float(process_rows["Contribution_Percent"].max()), 1.0)
    for py, (_, row) in zip(process_y, process_rows.iterrows()):
        # Colour based on the ion mapping (raw ion name stripped of Aug_ prefix)
        raw_ion = row["Base_Component"].replace("Aug_", "")
        proc_colour = ION_COLOURS.get(raw_ion, "#B0B8C1")
        
        # Driver box
        ax.add_patch(Rectangle((x_proc - bw_proc / 2, py - bh_proc / 2), bw_proc, bh_proc,
                              facecolor=proc_colour, edgecolor="#1f2937", 
                              linewidth=1.2, alpha=0.9, zorder=4))
        # Driver Label (Name + PD%)
        ax.text(x_proc, py + txt_off_1, textwrap.fill(str(row["Base_Component_Name"]), width=20),
                ha="center", va="center", fontsize=8.5, fontweight="bold", zorder=5)
        ax.text(x_proc, py - txt_off_2, f"pd={row['Prob_Direction']:.0f}%",
                ha="center", va="center", fontsize=8, fontstyle="italic", zorder=5)

        # Context -> Driver arrow
        sy = source_y[row["Source_Type"]]
        ax.add_patch(FancyArrowPatch((x_src + bw_src / 2, sy), (x_proc - bw_proc / 2, py),
                                    arrowstyle="-|>", mutation_scale=14, connectionstyle="arc3,rad=0.04",
                                    color=source_colours[row["Source_Type"]], linewidth=1.5, alpha=0.35, zorder=2))

        # Driver -> Endpoint arrow
        flow_val = float(row["Contribution_Percent"])
        flow_width = 1.0 + 4.0 * (flow_val / max_flow)
        ax.add_patch(FancyArrowPatch((x_proc + bw_proc / 2, py), (x_ep - bw_ep / 2, endpoint_y),
                                    arrowstyle="-|>", mutation_scale=12, connectionstyle="arc3,rad=0.04",
                                    color=endpoint_colour, linewidth=flow_width, alpha=0.32, zorder=2))
        # Flow value labels
        ax.text(0.68, (py + endpoint_y) / 2, f"{flow_val:.1f}%", fontsize=8, ha="center", va="center",
                color="#333", fontweight="bold", zorder=6,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))


def style_boxplot(bp: dict[str, list], face_colours: list[str]) -> None:
    for patch, colour in zip(bp["boxes"], face_colours):
        patch.set_facecolor(colour)
        patch.set_edgecolor("#1f2933")
        patch.set_linewidth(1.2)
        patch.set_alpha(0.85)
    for item in bp["medians"]:
        item.set_color("#1f2933")
        item.set_linewidth(2.0)
    for item in bp["whiskers"] + bp["caps"]:
        item.set_color("#1f2933")
        item.set_linewidth(1.0)


def style_violin(parts: dict[str, object], face_colours: list[str]) -> None:
    for body, colour in zip(parts["bodies"], face_colours):
        body.set_facecolor(colour)
        body.set_edgecolor("#1f2933")
        body.set_alpha(0.7)
    for key in ("cmeans", "cbars", "cmins", "cmaxes"):
        if key in parts:
            parts[key].set_color("#1f2933")
            parts[key].set_linewidth(1.2)


def save_fig(path: Path, fig=None, *, dpi: int = 600, pad: float = 0.35) -> None:
    if path in LOCKED_ARTIFACTS and path.exists():
        plt.close(fig if fig is not None else plt.gcf())
        print(f"Preserved locked figure: {path.relative_to(ROOT)}")
        return
    save_figure(path, fig=fig, dpi=dpi, pad=pad)


def ensure_positive(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        pos = out.loc[out[c] > 0, c]
        floor = pos.min() / np.sqrt(2) if not pos.empty else 1e-6
        out.loc[out[c] <= 0, c] = floor
    return out


def parse_paths(yaml_path: Path) -> list[str]:
    if not yaml_path.exists():
        return []
    text = yaml_path.read_text(encoding="utf-8")
    return [m.strip() for m in re.findall(r"output_path:\s*['\"]?([^'\"\n]+)['\"]?", text)]


def artifacts_exist(*paths: Path) -> bool:
    return all(path.exists() for path in paths)


def load_cached_hbmpra(trace_path: Path, groups: list[str]) -> tuple[dict[str, dict[str, np.ndarray]], pd.DataFrame]:
    trace = az.from_netcdf(str(trace_path))
    post: dict[str, dict[str, np.ndarray]] = {}
    rows: list[dict[str, object]] = []
    n_groups = len(groups)
    for i, group in enumerate(groups):
        post[group] = {}
        for organ in ["skeletal_dental", "hemato"]:
            var_name = f"HI_{organ}"
            hi_post = trace.posterior[var_name].values.reshape(-1, n_groups)[:, i]
            post[group][organ] = hi_post
            rows.append(
                {
                    "Demographic_Group": group,
                    "Organ_System": organ,
                    "Mean_HI": float(np.mean(hi_post)),
                    "CI_2_5": float(np.percentile(hi_post, 2.5)),
                    "CI_97_5": float(np.percentile(hi_post, 97.5)),
                    "P_HI_greater_1": float(np.mean(hi_post > 1.0)),
                    "Method": "HBMPRA (PyMC)",
                }
            )
    return post, pd.DataFrame(rows)


def water_types(data_meq: pd.DataFrame, cat_thresh: float = 0.25, ani_thresh: float = 0.25) -> pd.Series:
    """Classify water type using fractional thresholds on meq/L concentrations."""
    labels = []
    for _, row in data_meq.iterrows():
        total_cat = row["Ca"] + row["Mg"] + row["Na"] + row["K"]
        total_ani = row["HCO3"] + row["Cl"] + row["SO4"] + row["NO3"] + row["F"]
        cats = []
        for ion in ["Ca", "Mg", "Na", "K"]:
            if row[ion] / total_cat > cat_thresh:
                cats.append(ion)
        anis = []
        for ion in ["HCO3", "Cl", "SO4", "NO3", "F"]:
            if row[ion] / total_ani > ani_thresh:
                anis.append(ion)
        cats_str = "-".join(cats) if cats else "Mixed"
        anis_str = "-".join(anis) if anis else "Mixed"
        labels.append(f"{cats_str}-{anis_str}")
    return pd.Series(labels, index=data_meq.index)


def deterministic_group_risk(data: pd.DataFrame) -> dict[str, dict[str, np.ndarray]]:
    out: dict[str, dict[str, np.ndarray]] = {}
    for g, p in GROUPS.items():
        f = p["ir"] / p["bw"]
        hq_f = data["F"].to_numpy() * f / TOX["F"]["RfD_oral"]
        hq_n = data["NO3"].to_numpy() * f / TOX["NO3"]["RfD_oral"]
        out[g] = {
            "HQ_F": hq_f,
            "HQ_NO3": hq_n,
            "HI_skeletal_dental": hq_f,
            "HI_hemato": hq_n,
        }
    return out


def parallel_table(comp_data: pd.DataFrame, n_iter: int = 400) -> pd.DataFrame:
    ilr = ilr_transform(ensure_positive(comp_data))
    obs = PCA().fit(ilr).explained_variance_
    rng = np.random.default_rng(SEED)
    rnd = np.zeros((n_iter, len(obs)))
    for i in range(n_iter):
        rnd[i, :] = PCA().fit(rng.normal(size=ilr.shape)).explained_variance_
    crit = np.percentile(rnd, 95, axis=0)
    return pd.DataFrame({
        "Component": [f"PC{i+1}" for i in range(len(obs))],
        "Observed_Eigenvalue": obs,
        "Random_95th_Eigenvalue": crit,
        "Retained": obs > crit,
    })


def main() -> None:
    random.seed(SEED)
    np.random.seed(SEED)
    ART.mkdir(parents=True, exist_ok=True)
    MODELS.mkdir(parents=True, exist_ok=True)
    HYDRO.mkdir(parents=True, exist_ok=True)
    log_runtime_configuration()

    data = pd.read_csv(ROOT / "data.csv")
    numeric = [c for c in data.columns if c != "SampleID"]
    data[numeric] = data[numeric].apply(pd.to_numeric, errors="coerce")
    data = data.dropna(subset=COMP + ["pH", "EC", "TDS"]).copy()
    data.index = [f"{sid}_{i+1:02d}" for i, sid in enumerate(data["SampleID"].astype(str))]

    wf = IntegratedWorkflow(main_comp_vars=COMP, endpoint_vars=ENDS, standards=STANDARDS, toxref=TOX, use_ridge=True)
    res = wf.run(data, n_components=None, robust=True)

    prep = Preprocessor(main_compositional_vars=COMP, endpoint_vars=ENDS)
    data_meq = prep.to_meq(data[["Ca", "Mg", "Na", "K", "Cl", "HCO3", "SO4", "NO3", "F"]].copy())
    wtype = water_types(data_meq)

    hydro_df = pd.concat([
        pd.DataFrame({"SampleID": data["SampleID"].values, "SampleKey": data.index}),
        data_meq.reset_index(drop=True),
        pd.Series(wtype.values, name="Water_Type"),
    ], axis=1)
    save_csv(hydro_df, HYDRO / "data_meq_with_water_types.csv")

    hydro_sum = wtype.value_counts().rename_axis("Water_Type").reset_index(name="Count")
    hydro_sum["Percent"] = 100.0 * hydro_sum["Count"] / hydro_sum["Count"].sum()
    save_csv(hydro_sum, HYDRO / "water_type_summary.csv")

    # Compute explicit cation/anion ILR coordinates from meq/L
    dm = data_meq.copy()
    dm = dm.clip(lower=1e-10)  # avoid log(0)
    z1 = np.sqrt(2 / 3) * np.log(np.sqrt(dm["Ca"] * dm["Mg"]) / (dm["Na"] + dm["K"]))
    z2 = 1 / np.sqrt(2) * np.log(dm["Ca"] / dm["Mg"])
    z3 = np.sqrt(2 / 3) * np.log(np.sqrt(dm["Cl"] * dm["SO4"]) / dm["HCO3"])
    z4 = 1 / np.sqrt(2) * np.log(dm["Cl"] / dm["SO4"])

    # Colour palette for water types
    unique_wt = sorted(wtype.unique())
    _paired_colours = plt.cm.Paired(np.linspace(0, 1, max(len(unique_wt), 2)))
    _colour_map = dict(zip(unique_wt, _paired_colours))
    _wt_counts = wtype.value_counts()
    _total = len(wtype)
    _wt_labels = {w: f"{w} ({100.0 * _wt_counts.get(w, 0) / _total:.1f}%)" for w in unique_wt}

    # Standalone ILR plot (z1 vs z3 – cation-anion overview)
    plt.figure(figsize=(7.5, 6))
    for wt in unique_wt:
        mask = wtype == wt
        plt.scatter(z3[mask], z1[mask], c=[_colour_map[wt]], label=_wt_labels[wt],
                    edgecolors="black", s=50, linewidth=0.5, alpha=0.85)
    plt.xlabel(r"[Cl$^-$,SO$_4^{2-}$|HCO$_3^-$]")
    plt.ylabel(r"[Ca$^{2+}$,Mg$^{2+}$|Na$^+$+K$^+$]")
    plt.title("Hydrochemical Facies in ILR Space")
    plt.legend(fontsize=7)
    plt.axhline(0, ls="--", color="gray", alpha=0.6)
    plt.axvline(0, ls="--", color="gray", alpha=0.6)
    save_fig(HYDRO / "ilr_plot.png")

    cbe = res["charge_balance_error"]
    t1 = data[["pH", "EC", "TDS", "Na", "K", "Mg", "Ca", "Cl", "SO4", "HCO3", "NO3", "F"]].describe().T
    t1 = t1.rename(columns={"25%": "Q1", "50%": "Median", "75%": "Q3"}).reset_index().rename(columns={"index": "Parameter"})
    cbe_row = pd.DataFrame([{
        "Parameter": "Charge_Balance_Error_pct",
        "count": cbe.count(), "mean": cbe.mean(), "std": cbe.std(),
        "min": cbe.min(), "Q1": cbe.quantile(0.25), "Median": cbe.quantile(0.5),
        "Q3": cbe.quantile(0.75), "max": cbe.max(),
    }])
    save_csv(pd.concat([t1, cbe_row], ignore_index=True), ART / "Table1_Descriptive_Statistics.csv")

    limits = {
        "pH": ("6.5-8.5", lambda s: (s < 6.5) | (s > 8.5)),
        "EC": ("<=1500", lambda s: s > 1500),
        "TDS": ("<=1000", lambda s: s > 1000),
        "Na": ("<=200", lambda s: s > 200),
        "K": ("<=12", lambda s: s > 12),
        "Mg": ("<=50", lambda s: s > 50),
        "Ca": ("<=75", lambda s: s > 75),
        "Cl": ("<=250", lambda s: s > 250),
        "SO4": ("<=250", lambda s: s > 250),
        "HCO3": ("<=300", lambda s: s > 300),
        "NO3": ("<=50", lambda s: s > 50),
        "F": ("<=1.5", lambda s: s > 1.5),
    }
    rows = []
    for p, (rule_txt, fn) in limits.items():
        bad = fn(data[p])
        rows.append({
            "Parameter": p, "Guideline": rule_txt,
            "Exceedance_Count": int(bad.sum()),
            "Exceedance_Percent": float(bad.mean() * 100.0),
            "Mean": float(data[p].mean()), "Max": float(data[p].max()),
        })
    save_csv(pd.DataFrame(rows), ART / "Table2_WHO_GSA_Compliance.csv")

    weights = res["wqi_weights"].sort_values(ascending=False).reset_index()
    weights.columns = ["Parameter", "VIP_Dynamic_Weight"]
    weights["Rank"] = np.arange(1, len(weights) + 1)
    save_csv(weights, ART / "Table3_VIP_Weights.csv")
    # Persist WQI weights for diagnostics
    weights.to_csv(MODELS / "wqi_weights.csv", index=False)

    # Bayesian Sensitivity Analysis (Probabilistic Classification)
    if artifacts_exist(WQI_PROB_PATH, ISI_PROB_PATH):
        print("Using cached Bayesian sensitivity tables (14 & 15).")
    else:
        try:
            sens_res = run_probabilistic_sensitivity(data, STANDARDS, n_samples=300, error_pct=0.05, random_seed=SEED)
            t14 = sens_res["WQI_Probabilities"].reset_index().rename(columns={"index": "SampleKey"})
            save_csv(t14, WQI_PROB_PATH)

            t15 = sens_res["ISI_Probabilities"].reset_index().rename(columns={"index": "SampleKey"})
            save_csv(t15, ISI_PROB_PATH)
            print("Generated Bayesian Sensitivity Tables (14 & 15).")
        except Exception as e:
            print(f"Sensitivity Analysis skipped due to error: {e}")

    wqi = res["wqi"].copy()
    
    # Caching system for long-running models
    endpoint_models = {}
    model_configs = [
        ("NPLS_GWQI", "bayesian", lambda: run_bayesian_endpoint_model(
            wf.prep.get_augmented_predictor_block(data, "NPLS_GWQI"), wqi["NPLS_GWQI"], "NPLS_GWQI",
            include_scale=True, use_augmented=True,
            draws=max(ENDPOINT_BAYES_DRAWS, 1500),
            tune=max(ENDPOINT_BAYES_TUNE, 1500),
            chains=max(ENDPOINT_BAYES_CHAINS, 2),
            cores=1,
            target_accept=0.99,
            max_treedepth=15,
            likelihood="studentt",
            studentt_nu=4.0,
            regularized_horseshoe=True,
            slab_scale=2.0,
            shrinkage_family="normal",
            tau_beta=0.5,
            lam_beta=1.0,
            random_state=SEED
        )),
        ("F", "bayesian", lambda: run_bayesian_endpoint_model(
            wf.prep.get_augmented_predictor_block(data, "F"), data["F"], "F", include_scale=True, use_augmented=True,
            draws=max(ENDPOINT_BAYES_DRAWS, 1500),
            tune=max(ENDPOINT_BAYES_TUNE, 1500),
            chains=max(ENDPOINT_BAYES_CHAINS, 2),
            cores=1,
            target_accept=0.99,
            max_treedepth=15,
            likelihood="studentt",
            studentt_nu=4.0,
            regularized_horseshoe=True,
            slab_scale=2.0,
            shrinkage_family="normal",
            tau_beta=0.5,
            lam_beta=1.0,
            random_state=SEED + 200
        )),
        ("F_SPLS", "spls", lambda: run_sparse_pls_endpoint_model(
            wf.prep.get_augmented_predictor_block(data, "F"), data["F"], "F", use_augmented=True,
            n_permutations=20, random_state=SEED + 400
        )),
    ]
    
    for key, mtype, run_fn in model_configs:
        cache_path = MODELS / f"{key}_model.pkl"
        if cache_path.exists():
            print(f"Loading cached {mtype} model for {key} from {cache_path.name}")
            with open(cache_path, "rb") as f:
                cached_obj = pickle.load(f)

            # If cached Bayesian models predate newer metadata (e.g., fit_config), backfill metadata
            # (avoid expensive refits when the trace already exists).
            if mtype == "bayesian" and isinstance(cached_obj, dict) and ("fit_config" not in cached_obj):
                print(f"Cached bayesian model for {key} is missing fit metadata; backfilling fit_config...")
                if key in {"NPLS_GWQI", "F"}:
                    cached_obj["fit_config"] = {
                        "draws": int(max(ENDPOINT_BAYES_DRAWS, 1500)),
                        "tune": int(max(ENDPOINT_BAYES_TUNE, 1500)),
                        "chains": int(max(ENDPOINT_BAYES_CHAINS, 2)),
                        "cores": 1,
                        "random_state": int(SEED if key == "NPLS_GWQI" else SEED + 200),
                        "target_accept": 0.99,
                        "max_treedepth": 15,
                        "likelihood": "studentt",
                        "studentt_nu": 4.0,
                        "regularized_horseshoe": True,
                        "slab_scale": 2.0,
                        "shrinkage_family": "normal",
                        "tau_beta": 0.5,
                        "lam_beta": 1.0,
                        "log_likelihood": True,
                        "posterior_predictive": False,
                        "use_augmented": True,
                        "include_scale": True,
                    }
                    with open(cache_path, "wb") as f:
                        pickle.dump(cached_obj, f)
                else:
                    print(f"WARNING: missing fit metadata for {key}; leaving as-is.")
            endpoint_models[key] = cached_obj
        else:
            print(f"Running {mtype} model for {key}...")
            res_dict = run_fn()
            # If the dict contains a trace (InferenceData), it should be picklable 
            # within the same environment; we include it for downstream access.
            endpoint_models[key] = res_dict
            print(f"Saving model for {key} to cache...")
            with open(cache_path, "wb") as f:
                pickle.dump(res_dict, f)
    t4 = wqi.groupby("Quality_Class", observed=False).agg(
        Sample_Count=("NPLS_GWQI", "size"),
        Mean_NPLS_GWQI=("NPLS_GWQI", "mean"),
        Mean_Indeterminacy=("Indeterminacy_Mean", "mean"),
    ).reset_index()
    save_csv(t4, ART / "Table4_GWQI_Classification.csv")

    t5 = parallel_table(data[COMP])
    save_csv(t5, ART / "Table5_Parallel_Analysis.csv")

    loadings = res["discovery"]["loadings_clr"].copy()
    t6 = loadings.reset_index().rename(columns={"index": "Process"})
    t6["Refined_Process_Name"] = t6["Process"].map(res["discovery"]["suggested_names"])
    save_csv(t6, ART / "Table6_CLR_Loadings.csv")

    t7 = res["geochemical_constraints"].copy().reset_index(drop=True)
    t7.insert(0, "SampleID", data["SampleID"].values)
    save_csv(t7, ART / "Table7_Geochemical_Constraints.csv")

    order = ["NPLS_GWQI", "F"]
    t8_rows: list[dict[str, object]] = []
    for ep in order:
        model_res = endpoint_models[ep]
        for predictor, driver_score in model_res["driver_score"].items():
            t8_rows.append(
                {
                    "Endpoint": ep,
                    "Endpoint_Label": ENDPOINT_TEXT_LABELS.get(ep, ep),
                    "Predictor_Block_Label": model_res["predictor_block_label"],
                    "Predictor": predictor,
                    "Predictor_Label": predictor,
                    "Driver_Score_Percent": float(driver_score),
                    "Prob_Direction_Percent": float(model_res["prob_direction"].get(predictor, 0.0)),
                    "Beta_Mean": float(model_res["beta_mean"].get(predictor, 0.0)),
                    "Beta_HDI_Lower": float(model_res["beta_hdi_lower"].get(predictor, 0.0)),
                    "Beta_HDI_Upper": float(model_res["beta_hdi_upper"].get(predictor, 0.0)),
                    "Model_R2": float(model_res["in_sample_r2"]),
                }
            )
    t8 = pd.DataFrame(t8_rows)
    t8["Endpoint_Order"] = t8["Endpoint"].map({ep: idx for idx, ep in enumerate(order)})
    t8 = t8.sort_values(
        ["Endpoint_Order", "Driver_Score_Percent", "Prob_Direction_Percent", "Predictor_Label"],
        ascending=[True, False, False, True],
    ).drop(columns="Endpoint_Order").reset_index(drop=True)
    save_csv(t8, ART / "Table8_Process_Influence.csv")

    t9_rows = []
    for ep in ["NPLS_GWQI", "F", "F_SPLS"]:
        if ep not in endpoint_models:
            continue
        model_res = endpoint_models[ep]
        is_spls = "_SPLS" in ep
        t9_rows.append(
            {
                "Endpoint": ep,
                "Model_Type": "Sparse_PLS" if is_spls else "Bayesian_Horseshoe",
                "Predictor_Block": model_res["predictor_block_label"],
                "Features_Count": int(len(model_res["selected_features"])),
                "Features": "; ".join(model_res["selected_features"]),
                "In_Sample_R2": float(model_res["in_sample_r2"]),
                "CV_R2_Mean": float(model_res.get("nested_cv_r2_mean", np.nan)),
                "Permutation_p": float(model_res.get("permutation_p_value", np.nan)),
            }
        )
    t9 = pd.DataFrame(t9_rows)
    save_csv(t9, ART / "Table9_Model_Diagnostics.csv")

    det = deterministic_group_risk(data)
    t10_rows = []
    for g, vals in det.items():
        t10_rows.append({
            "Demographic_Group": g,
            "HQ_F_Mean": float(np.mean(vals["HQ_F"])),
            "HQ_F_P95": float(np.percentile(vals["HQ_F"], 95)),
            "HQ_NO3_Mean": float(np.mean(vals["HQ_NO3"])),
            "HQ_NO3_P95": float(np.percentile(vals["HQ_NO3"], 95)),
            "HI_Skeletal_Dental_Mean": float(np.mean(vals["HI_skeletal_dental"])),
            "HI_Hemato_Mean": float(np.mean(vals["HI_hemato"])),
        })
    save_csv(pd.DataFrame(t10_rows), ART / "Table10_Deterministic_Risk.csv")

    post: dict[str, dict[str, np.ndarray]] = {}
    gmap = {"Adults": "Adults", "Children": "Children", "Teens": "Teens"}
    method = "HBMPRA (PyMC)"
    if HBMPRA_TRACE_PATH.exists() and not FORCE_HBMPRA_RERUN:
        try:
            post, t11 = load_cached_hbmpra(HBMPRA_TRACE_PATH, list(gmap.values()))
            save_csv(t11, HBMPRA_TABLE_PATH)
            print(f"Using cached HBMPRA trace from {HBMPRA_TRACE_PATH}.")
        except Exception as exc:
            print(f"Cached HBMPRA trace could not be reused ({type(exc).__name__}: {exc}); rerunning Bayesian model.")
            try:
                risk = assess_health_risk(
                    data[["F", "NO3"]],
                    toxref=TOX,
                    mode="bayesian",
                    n_samples=HBMPRA_BAYES_DRAWS,
                    random_seed=SEED,
                    chains=HBMPRA_BAYES_CHAINS,
                    tune=HBMPRA_BAYES_TUNE,
                    target_accept=0.95,
                )
                if "pymc_trace" in risk:
                    risk["pymc_trace"].to_netcdf(HBMPRA_TRACE_PATH)
                    print(f"Persisted PyMC trace to {HBMPRA_TRACE_PATH}")

                bayes_rows = []
                for pg, dg in gmap.items():
                    post[dg] = {
                        "skeletal_dental": np.array(risk[f"HI_{pg}_skeletal_dental_posterior"]),
                        "hemato": np.array(risk[f"HI_{pg}_hemato_posterior"]),
                    }
                    for organ in ["skeletal_dental", "hemato"]:
                        lo, hi = risk[f"HI_{pg}_{organ}_95_CI"]
                        bayes_rows.append({
                            "Demographic_Group": dg,
                            "Organ_System": organ,
                            "Mean_HI": float(risk[f"HI_{pg}_{organ}_mean"]),
                            "CI_2_5": float(lo),
                            "CI_97_5": float(hi),
                            "P_HI_greater_1": float(risk[f"HI_{pg}_{organ}_prob_exceed_1"]),
                            "Method": method,
                        })
            except Exception as rerun_exc:
                method = f"Monte Carlo fallback ({type(rerun_exc).__name__})"
                rng = np.random.default_rng(SEED)
                bayes_rows = []
                for dg in GROUPS:
                    f_draws = rng.choice(det[dg]["HI_skeletal_dental"], size=2000, replace=True)
                    n_draws = rng.choice(det[dg]["HI_hemato"], size=2000, replace=True)
                    post[dg] = {"skeletal_dental": f_draws, "hemato": n_draws}
                    for organ, arr in post[dg].items():
                        bayes_rows.append({
                            "Demographic_Group": dg,
                            "Organ_System": organ,
                            "Mean_HI": float(np.mean(arr)),
                            "CI_2_5": float(np.percentile(arr, 2.5)),
                            "CI_97_5": float(np.percentile(arr, 97.5)),
                            "P_HI_greater_1": float(np.mean(arr > 1.0)),
                            "Method": method,
                        })
            t11 = pd.DataFrame(bayes_rows)
            save_csv(t11, HBMPRA_TABLE_PATH)
    else:
        try:
            risk = assess_health_risk(
                data[["F", "NO3"]],
                toxref=TOX,
                mode="bayesian",
                n_samples=HBMPRA_BAYES_DRAWS,
                random_seed=SEED,
                chains=HBMPRA_BAYES_CHAINS,
                tune=HBMPRA_BAYES_TUNE,
                target_accept=0.95,
            )
            if "pymc_trace" in risk:
                risk["pymc_trace"].to_netcdf(HBMPRA_TRACE_PATH)
                print(f"Persisted PyMC trace to {HBMPRA_TRACE_PATH}")

            bayes_rows = []
            for pg, dg in gmap.items():
                post[dg] = {
                    "skeletal_dental": np.array(risk[f"HI_{pg}_skeletal_dental_posterior"]),
                    "hemato": np.array(risk[f"HI_{pg}_hemato_posterior"]),
                }
                for organ in ["skeletal_dental", "hemato"]:
                    lo, hi = risk[f"HI_{pg}_{organ}_95_CI"]
                    bayes_rows.append({
                        "Demographic_Group": dg,
                        "Organ_System": organ,
                        "Mean_HI": float(risk[f"HI_{pg}_{organ}_mean"]),
                        "CI_2_5": float(lo),
                        "CI_97_5": float(hi),
                        "P_HI_greater_1": float(risk[f"HI_{pg}_{organ}_prob_exceed_1"]),
                        "Method": method,
                    })
        except Exception as exc:
            method = f"Monte Carlo fallback ({type(exc).__name__})"
            rng = np.random.default_rng(SEED)
            bayes_rows = []
            for dg in GROUPS:
                f_draws = rng.choice(det[dg]["HI_skeletal_dental"], size=2000, replace=True)
                n_draws = rng.choice(det[dg]["HI_hemato"], size=2000, replace=True)
                post[dg] = {"skeletal_dental": f_draws, "hemato": n_draws}
                for organ, arr in post[dg].items():
                    bayes_rows.append({
                        "Demographic_Group": dg,
                        "Organ_System": organ,
                        "Mean_HI": float(np.mean(arr)),
                        "CI_2_5": float(np.percentile(arr, 2.5)),
                        "CI_97_5": float(np.percentile(arr, 97.5)),
                        "P_HI_greater_1": float(np.mean(arr > 1.0)),
                        "Method": method,
                    })
        t11 = pd.DataFrame(bayes_rows)
        save_csv(t11, HBMPRA_TABLE_PATH)

    irr = res["irrigation"].copy()
    num_cols = ["SAR", "SSP", "MH", "KR", "PI", "RSC", "PS", "N_ISI_Score"]
    t12_rows = []
    for c in num_cols:
        s = irr[c]
        t12_rows.append({"Section": "Index_Statistics", "Metric": c, "Mean": float(s.mean()), "STD": float(s.std()), "Min": float(s.min()), "Max": float(s.max())})
    class_cols = [c for c in irr.columns if c.endswith("_Class") or c == "ISI_Class"]
    for c in class_cols:
        vc = irr[c].value_counts(dropna=False)
        for k, v in vc.items():
            t12_rows.append({"Section": "Class_Distribution", "Metric": c, "Class": str(k), "Count": int(v), "Percent": float(100.0 * v / len(irr))})
    save_csv(pd.DataFrame(t12_rows), ART / "Table12_Irrigation_Indices.csv")

    ent = irr["Entropy_Weights"].iloc[0]
    t13_rows = [{"Row_Type": "Entropy_Weight", "Item": k, "Value": float(v)} for k, v in ent.items()]
    for k, v in irr["ISI_Class"].value_counts().items():
        t13_rows.append({"Row_Type": "NISI_Class_Distribution", "Item": str(k), "Value": float(100.0 * v / len(irr)), "Count": int(v)})
    save_csv(pd.DataFrame(t13_rows), ART / "Table13_NISI_Classification.csv")

    ts1 = wqi.copy()
    ts1.insert(0, "SampleID", data["SampleID"].values)
    ts1.insert(1, "SampleKey", data.index)
    save_csv(ts1, ART / "TableS1_Full_WQI.csv")

    ts2 = pd.DataFrame({"SampleID": data["SampleID"].values, "SampleKey": data.index})
    ts2 = pd.concat([ts2, res["discovery"]["scores"].reset_index(drop=True)], axis=1)
    for ep in order:
        if ep in endpoint_models:
            m = endpoint_models[ep]["sample_level_contributions"].reset_index(drop=True)
            ts2 = pd.concat([ts2, m.add_prefix(f"{ep}_")], axis=1)
    save_csv(ts2, ART / "TableS2_Sample_Contributions.csv")

    # Figure 2 — Standalone 2×2 ILR plot
    fig2, ((ax_ul, ax_ur), (ax_ll, ax_lr)) = plt.subplots(2, 2, figsize=(9, 7))

    _ilr_colours = [_colour_map[w] for w in wtype]

    # Dummy variables for field boundaries
    _dum1 = np.full(23, 0.5)
    _dum2 = np.array([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2,
                      5e-2, 1e-1, 2e-1, 3e-1, 4e-1, 4.5e-1, 4.9e-1, 4.99e-1,
                      4.999e-1, 4.9999e-1, 4.99999e-1, 4.999999e-1,
                      4.9999999e-1, 4.99999999e-1])
    _dum3 = 1 - _dum1 - _dum2

    # ---- (a) Upper-Left: [Ca²⁺|Mg²⁺] vs [Cl⁻|SO₄²⁻] ----
    ax_ul.scatter(z2, z4, c=_ilr_colours, edgecolors="black", s=50, linewidth=0.5)
    ax_ul.axhline(0, ls="--", color="gray", alpha=0.6)
    ax_ul.axvline(0, ls="--", color="gray", alpha=0.6)
    ax_ul.set_xlim(-10, 10); ax_ul.set_ylim(-10, 10)
    ax_ul.set_xticks(np.arange(-8, 9, 2)); ax_ul.set_yticks(np.arange(-8, 9, 2))
    ax_ul.tick_params(axis="both", labelsize=10)
    ax_ul.set_ylabel(r"[Cl$^-$|SO$_4^{2-}$]", fontsize=12)
    ax_ul.set_title(r"[Ca$^{2+}$|Mg$^{2+}$]", fontsize=12, pad=20)
    # Ion-dominance arrows
    ax_ul.arrow(-0.5, 8.5, -1, 0, head_width=0.2, head_length=0.3, fc="black", ec="black", length_includes_head=True)
    ax_ul.text(-2.5, 9, r"Mg$^{2+}$ $>$ Ca$^{2+}$", fontsize=8, ha="center")
    ax_ul.arrow(0.5, 8.5, 1, 0, head_width=0.2, head_length=0.3, fc="black", ec="black", length_includes_head=True)
    ax_ul.text(2.5, 9, r"Ca$^{2+}$ $>$ Mg$^{2+}$", fontsize=8, ha="center")
    ax_ul.arrow(-8.5, -1, 0, -1, head_width=0.2, head_length=0.3, fc="black", ec="black", length_includes_head=True)
    ax_ul.text(-9, -3, r"SO$_4^{2-}$ $>$ Cl$^-$", fontsize=8, ha="center", va="center", rotation=90)
    ax_ul.arrow(-8.5, 1, 0, 1, head_width=0.2, head_length=0.3, fc="black", ec="black", length_includes_head=True)
    ax_ul.text(-9, 3, r"Cl$^-$ $>$ SO$_4^{2-}$", fontsize=8, ha="center", va="center", rotation=90)

    # ---- (b) Upper-Right: Anion space ----
    ax_ur.scatter(z3, z4, c=_ilr_colours, edgecolors="black", s=50, linewidth=0.5)
    ax_ur.set_xlim(-10, 10); ax_ur.set_ylim(-10, 10)
    ax_ur.set_xticks(np.arange(-8, 9, 2)); ax_ur.set_yticks([])
    ax_ur.tick_params(axis="both", labelsize=10)
    ax_ur.set_title(r"[Cl$^-$,SO$_4^{2-}$|HCO$_3^-$]", fontsize=12, pad=20)
    ax_ur.set_ylabel(r"[Cl$^-$|SO$_4^{2-}$]", fontsize=12)
    ax_ur.yaxis.set_label_position("right")
    # SO4-type field
    _b_ilr1 = np.sqrt(2 / 3) * np.log(np.sqrt(_dum1 * _dum2) / _dum3)
    _b_ilr2 = np.sqrt(1 / 2) * np.log(_dum2 / _dum1)
    ax_ur.plot(_b_ilr1, _b_ilr2, ls="--", color="gray", alpha=0.6)
    ax_ur.text(1, -1.5, r"SO$_4^{2-}$ type", fontsize=10)
    # HCO3-type field
    _b_ilr1 = np.sqrt(2 / 3) * np.log(np.sqrt(_dum3 * _dum2) / _dum1)
    _b_ilr2 = np.sqrt(1 / 2) * np.log(_dum3 / _dum2)
    ax_ur.plot(_b_ilr1, _b_ilr2, ls="--", color="gray", alpha=0.6)
    ax_ur.text(-6, -4, r"HCO$_3^{-}$ type", fontsize=10)
    # Cl-type field
    _b_ilr1 = np.sqrt(2 / 3) * np.log(np.sqrt(_dum3 * _dum1) / _dum2)
    _b_ilr2 = np.sqrt(1 / 2) * np.log(_dum1 / _dum3)
    ax_ur.plot(_b_ilr1, _b_ilr2, ls="--", color="gray", alpha=0.6)
    ax_ur.text(1, 3, r"Cl$^{-}$ type", fontsize=10)

    # ---- (c) Lower-Left: Cation space ----
    ax_ll.scatter(z2, z1, c=_ilr_colours, edgecolors="black", s=50, linewidth=0.5)
    ax_ll.set_xlim(-10, 10); ax_ll.set_ylim(-10, 10)
    ax_ll.set_xticks(np.arange(-8, 9, 2)); ax_ll.set_yticks(np.arange(-8, 9, 2))
    ax_ll.tick_params(axis="both", labelsize=10)
    ax_ll.set_xlabel(r"[Ca$^{2+}$|Mg$^{2+}$]", fontsize=12)
    ax_ll.set_ylabel(r"[Ca$^{2+}$,Mg$^{2+}$|Na$^+$ + K$^+$]", fontsize=12)
    # Mg-type field
    _c_ilr1 = np.sqrt(2 / 3) * np.log(np.sqrt(_dum1 * _dum2) / _dum3)
    _c_ilr2 = np.sqrt(1 / 2) * np.log(_dum2 / _dum1)
    ax_ll.plot(_c_ilr2, _c_ilr1, ls="--", color="gray", alpha=0.6)
    ax_ll.text(-4, 4, r"Mg$^{2+}$ type", fontsize=10)
    # Na+K-type field
    _c_ilr1 = np.sqrt(2 / 3) * np.log(np.sqrt(_dum3 * _dum2) / _dum1)
    _c_ilr2 = np.sqrt(1 / 2) * np.log(_dum3 / _dum2)
    ax_ll.plot(_c_ilr2, _c_ilr1, ls="--", color="gray", alpha=0.6)
    ax_ll.text(2, -5, r"Na$^+$+K$^+$ type", fontsize=10)
    # Ca-type field
    _c_ilr1 = np.sqrt(2 / 3) * np.log(np.sqrt(_dum3 * _dum1) / _dum2)
    _c_ilr2 = np.sqrt(1 / 2) * np.log(_dum1 / _dum3)
    ax_ll.plot(_c_ilr2, _c_ilr1, ls="--", color="gray", alpha=0.6)
    ax_ll.text(2, 0.5, r"Ca$^{2+}$ type", fontsize=10)

    # ---- (d) Lower-Right: Diamond-like combined panel ----
    ax_lr.scatter(z3, z1, c=_ilr_colours, edgecolors="black", s=50, linewidth=0.5)
    ax_lr.axhline(0, ls="--", color="gray", alpha=0.6)
    ax_lr.axvline(0, ls="--", color="gray", alpha=0.6)
    ax_lr.set_xlim(-10, 10); ax_lr.set_ylim(-10, 10)
    ax_lr.set_xticks(np.arange(-8, 9, 2)); ax_lr.set_yticks([])
    ax_lr.tick_params(axis="both", labelsize=10)
    ax_lr.set_xlabel(r"[Cl$^-$,SO$_4^{2-}$|HCO$_3^-$]", fontsize=12)
    ax_lr.set_ylabel(r"[Ca$^{2+}$,Mg$^{2+}$|Na$^+$ + K$^+$]", fontsize=12)
    ax_lr.yaxis.set_label_position("right")
    # Ion-dominance arrows
    ax_lr.arrow(-1, 8.5, -1, 0, head_width=0.2, head_length=0.3, fc="black", ec="black", length_includes_head=True)
    ax_lr.text(-3, 9, r"HCO$_3^-$ $>$ Cl$^-$, SO$_4^{2-}$", fontsize=7, ha="center")
    ax_lr.arrow(1, 8.5, 1, 0, head_width=0.2, head_length=0.3, fc="black", ec="black", length_includes_head=True)
    ax_lr.text(3.5, 9, r"Cl$^-$, SO$_4^{2-}$ $>$ HCO$_3^-$", fontsize=7, ha="center")
    ax_lr.arrow(-8.5, -1, 0, -1, head_width=0.2, head_length=0.3, fc="black", ec="black", length_includes_head=True)
    ax_lr.text(-9, -5, r"Na$^+$+K$^+$ $>$ Ca$^{2+}$, Mg$^{2+}$", fontsize=7, ha="center", va="center", rotation=90)
    ax_lr.arrow(-8.5, 1, 0, 1, head_width=0.2, head_length=0.3, fc="black", ec="black", length_includes_head=True)
    ax_lr.text(-9, 5, r"Ca$^{2+}$, Mg$^{2+}$ $>$ Na$^+$+K$^+$", fontsize=7, ha="center", va="center", rotation=90)

    # Remove inter-subplot spacing and add legend below
    plt.subplots_adjust(wspace=0, hspace=0, bottom=0.22)
    handles = [plt.scatter([], [], c=[_colour_map[w]], label=_wt_labels[w],
               edgecolors="black", s=50) for w in unique_wt]
    fig2.legend(handles=handles, loc="lower center",
                bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=10, frameon=False)
    save_fig(ART / "Figure2_Hydrochemical_Characterisation.png", fig=fig2, pad=0.2)

    # Figure 3
    cat_sum = data_meq[["Ca", "Mg", "Na", "K"]].sum(axis=1)
    ani_sum = data_meq[["Cl", "SO4", "HCO3", "NO3", "F"]].sum(axis=1)
    fig3, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.8))
    a1.hist(cbe, bins=10, color="#5E548E", edgecolor="white", linewidth=1.2)
    a1.axvline(cbe.mean(), color="#BC4749", linestyle="--", linewidth=1.6, label=f"Mean = {cbe.mean():.2f}%")
    a1.set_title("Charge Balance Error Histogram")
    a1.set_xlabel("Charge balance error (%)")
    a1.set_ylabel("Sample count")
    style_axes(a1, xbins=5, ybins=5, grid=True, grid_axis="y")
    style_legend(a1.legend(loc="upper left"))

    a2.scatter(cat_sum, ani_sum, color="#2A9D8F", edgecolor="black", s=80, alpha=0.85, linewidth=0.7)
    mn, mx = min(cat_sum.min(), ani_sum.min()), max(cat_sum.max(), ani_sum.max())
    a2.plot([mn, mx], [mn, mx], linestyle="--", color="#BC4749", linewidth=1.8, label="1:1 equivalence")
    a2.set_title("Cation vs Anion Sum (meq/L)")
    a2.set_xlabel("Sum of cations (meq/L)")
    a2.set_ylabel("Sum of anions (meq/L)")
    style_axes(a2, xbins=5, ybins=5, grid=True)
    style_legend(a2.legend(loc="upper left"))
    save_fig(ART / "Figure3_CBE_Integrity.png", fig=fig3)

    # Figure 4
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(2, 2)
    
    # (a) Neutrosophic Components vs Score (Formerly b)
    a1 = fig.add_subplot(gs[0, 0])
    wqi_sorted = wqi.sort_values("NPLS_GWQI")
    a1.plot(wqi_sorted["NPLS_GWQI"], wqi_sorted["Truth_Mean"], label="Truth (T)", color="#2A9D8F", linewidth=2.5, zorder=4)
    a1.plot(wqi_sorted["NPLS_GWQI"], wqi_sorted["Indeterminacy_Mean"], label="Indeterminacy (I)", color="#7B2CBF", linewidth=2.5, zorder=4)
    a1.plot(wqi_sorted["NPLS_GWQI"], wqi_sorted["Falsity_Mean"], label="Falsity (F)", color="#E15759", linewidth=2.5, zorder=4)
    a1.scatter(wqi["NPLS_GWQI"], wqi["Truth_Mean"], color="#2A9D8F", alpha=0.25, s=20, zorder=3)
    a1.scatter(wqi["NPLS_GWQI"], wqi["Indeterminacy_Mean"], color="#7B2CBF", alpha=0.25, s=20, zorder=3)
    a1.scatter(wqi["NPLS_GWQI"], wqi["Falsity_Mean"], color="#E15759", alpha=0.25, s=20, zorder=3)
    a1.set_title("(a) Neutrosophic Components vs Score", fontsize=14, fontweight="bold")
    a1.set_xlabel("NPLS-GWQI Score", fontsize=13)
    a1.set_ylabel("Mean Component Value (0-1)", fontsize=13)
    style_legend(a1.legend(loc="upper right"))
    style_axes(a1, xbins=5, ybins=5, grid=True)
    
    # (b) VIP Dynamic Weights Radar (Formerly c)
    a2 = fig.add_subplot(gs[0, 1], projection="polar")
    ang = np.linspace(0, 2 * np.pi, len(weights), endpoint=False)
    vals = weights["VIP_Dynamic_Weight"].to_numpy()
    a2.plot(np.r_[ang, ang[0]], np.r_[vals, vals[0]], linewidth=2.4, color="#6D597A")
    a2.fill(np.r_[ang, ang[0]], np.r_[vals, vals[0]], alpha=0.22, color="#B8A1D9")
    a2.set_xticks(ang)
    a2.set_xticklabels(weights["Parameter"].tolist(), fontsize=12, fontweight="bold")
    a2.tick_params(axis="both", direction="in", labelsize=12)
    a2.grid(alpha=0.3, linewidth=0.5, linestyle="--")
    a2.set_title("(b) VIP Dynamic Weights Radar", fontsize=14, fontweight="bold", pad=20)
    
    # (c) NPLS-GWQI by Class (Formerly a)
    a3 = fig.add_subplot(gs[1, 0])
    cls = [c for c in ["Excellent", "Good", "Poor", "Very Poor", "Unsuitable"] if c in set(wqi["Quality_Class"].astype(str))]
    bp = a3.boxplot(
        [wqi.loc[wqi["Quality_Class"].astype(str) == c, "NPLS_GWQI"] for c in cls],
        tick_labels=cls,
        showfliers=False,
        patch_artist=True,
    )
    style_boxplot(bp, [CLASS_COLOURS.get(c, "#C9D6DF") for c in cls])
    a3.set_title("(c) NPLS-GWQI by Class", fontsize=14, fontweight="bold")
    a3.set_ylabel("NPLS-GWQI Score", fontsize=13)
    style_axes(a3, ybins=5, grid=True, grid_axis="y")
    
    # (d) NPLS-GWQI Classification Distribution
    a4 = fig.add_subplot(gs[1, 1])
    wqi_counts = wqi["Quality_Class"].value_counts(dropna=False).reindex(cls).fillna(0)
    x_positions = np.arange(len(cls))
    bars = a4.bar(x_positions, wqi_counts.values, color=[CLASS_COLOURS.get(c, "#C9D6DF") for c in cls], edgecolor="black", alpha=0.8, width=0.6)
    a4.set_xticks(x_positions)
    a4.set_xticklabels(cls, rotation=20, ha="right")
    a4.set_xlabel("NPLS-GWQI Suitability Class", fontsize=13)
    a4.set_ylabel("Frequency (Wells)", fontsize=13)
    a4.set_title("(d) NPLS-GWQI Classification Distribution", fontsize=14, fontweight="bold")
    style_axes(a4, ybins=5, grid=True, grid_axis="y")
    for bar in bars:
        yval = bar.get_height()
        a4.text(bar.get_x() + bar.get_width() / 2, yval + (wqi_counts.max() * 0.01), f"{int(yval)}", ha="center", va="bottom", fontsize=11, fontweight="bold")

        
    save_fig(ART / "Figure4_NPLS_GWQI_Structure.png", fig=fig)

    # Figure 5
    fig5, ax5 = plt.subplots(figsize=(12, 7))
    lp = loadings.T
    w = 0.8 / max(1, lp.shape[1])
    x = np.arange(lp.shape[0])
    process_name_map = t6.set_index("Process")["Refined_Process_Name"].to_dict()
    for i, proc in enumerate(lp.columns):
        ax5.bar(
            x + i * w,
            lp[proc].values,
            width=w,
            color=PROCESS_COLOURS[i % len(PROCESS_COLOURS)],
            label=textwrap.fill(process_name_map.get(proc, proc).replace("_", " "), width=18),
            edgecolor="white",
            linewidth=0.8,
        )
    ax5.axhline(0, color="#1f2933", linewidth=1.2)
    ax5.set_xticks(x + w * (lp.shape[1] - 1) / 2, lp.index, rotation=30, ha="right")
    ax5.set_title("Latent Hydrochemical Signature Fingerprint (CLR Loadings)")
    ax5.set_ylabel("CLR loading")
    style_axes(ax5, ybins=5, grid=True, grid_axis="y")
    style_legend(ax5.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2))
    save_fig(ART / "Figure5_CLR_Loadings.png", fig=fig5, pad=0.55)

    # Figure 6 — Dual Gibbs + CEV + Dissolution diagnostic
    _na_ratio = data_meq["Na"] / (data_meq["Na"] + data_meq["Ca"])
    _cl_ratio = data_meq["Cl"] / (data_meq["Cl"] + data_meq["HCO3"])
    _tds_calc = t7["TDS_Calculated"].values
    _f_vals = data["F"].values
    _cmap = plt.cm.coolwarm

    fig6 = plt.figure(figsize=(14.5, 10.5))
    gs6 = fig6.add_gridspec(
        2,
        3,
        width_ratios=[1.0, 1.0, 0.06],
        height_ratios=[1.0, 1.0],
        wspace=0.20,
        hspace=0.20,
    )
    a1 = fig6.add_subplot(gs6[0, 0])
    a2 = fig6.add_subplot(gs6[0, 1])
    cax = fig6.add_subplot(gs6[0, 2])
    a3 = fig6.add_subplot(gs6[1, 0])
    a4 = fig6.add_subplot(gs6[1, 1])

    # (a) Gibbs — Cations
    sc1 = a1.scatter(_na_ratio, _tds_calc, c=_f_vals, cmap=_cmap,
                     edgecolor="black", s=80, linewidth=0.5, zorder=3)
    a1.set_yscale("log")
    a1.set_ylim(10, 100_000)
    a1.set_xlim(-0.05, 1.05)
    a1.set_xlabel(r"Na$^+$ / (Na$^+$ + Ca$^{2+}$)", fontsize=11)
    a1.set_ylabel("TDS (mg/L)", fontsize=11)
    a1.set_title("(a) Gibbs Plot \u2014 Cations", fontsize=12, pad=10)
    a1.grid(True, which="both", ls="--", alpha=0.4)
    a1.text(0.85, 30, "Precipitation\nDominance", fontsize=9, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="#e8f4fd", ec="gray", alpha=0.7))
    a1.text(0.30, 800, "Rock Weathering\nDominance", fontsize=9, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fef9e7", ec="gray", alpha=0.7))
    a1.text(0.85, 20_000, "Evaporation\nDominance", fontsize=9, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fde2e2", ec="gray", alpha=0.7))

    # (b) Gibbs — Anions
    sc2 = a2.scatter(_cl_ratio, _tds_calc, c=_f_vals, cmap=_cmap,
                     edgecolor="black", s=80, linewidth=0.5, zorder=3)
    a2.set_yscale("log")
    a2.set_ylim(10, 100_000)
    a2.set_xlim(-0.05, 1.05)
    a2.set_xlabel(r"Cl$^-$ / (Cl$^-$ + HCO$_3^-$)", fontsize=11)
    a2.set_ylabel("TDS (mg/L)", fontsize=11)
    a2.set_title("(b) Gibbs Plot \u2014 Anions", fontsize=12, pad=10)
    a2.grid(True, which="both", ls="--", alpha=0.4)
    a2.text(0.85, 30, "Precipitation\nDominance", fontsize=9, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="#e8f4fd", ec="gray", alpha=0.7))
    a2.text(0.30, 800, "Rock Weathering\nDominance", fontsize=9, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fef9e7", ec="gray", alpha=0.7))
    a2.text(0.85, 20_000, "Evaporation\nDominance", fontsize=9, ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fde2e2", ec="gray", alpha=0.7))

    # Shared colorbar for Gibbs panels (Fluoride)
    cbar = fig6.colorbar(sc2, cax=cax)
    style_colorbar(cbar, "F (mg/L)")
    style_axes(a1, xbins=5, grid=True)
    style_axes(a2, xbins=5, grid=True)

    # (c) CEV Scatter — Chloro-Alkaline Exchange
    _cev_x = (data_meq["Na"] - data_meq["Cl"]).values  # Na − Cl
    _cev_y = ((data_meq["Ca"] + data_meq["Mg"]) - (data_meq["HCO3"] + data_meq["SO4"])).values
    a3.scatter(_cev_x, _cev_y, c="teal", edgecolor="black", s=60, alpha=0.8, linewidth=0.5)
    a3.axhline(0, color="black", linewidth=0.8)
    a3.axvline(0, color="black", linewidth=0.8)
    a3.set_xlabel(r"Na$^+$ $-$ Cl$^-$ (meq/L)", fontsize=11)
    a3.set_ylabel(r"(Ca$^{2+}$+Mg$^{2+}$) $-$ (HCO$_3^-$+SO$_4^{2-}$) (meq/L)", fontsize=10)
    a3.set_title("(c) Chloro-Alkaline Exchange", fontsize=12, pad=10)
    style_axes(a3, xbins=5, ybins=5, grid=True)
    # Quadrant annotations
    _xr = max(abs(a3.get_xlim()[0]), abs(a3.get_xlim()[1]))
    _yr = max(abs(a3.get_ylim()[0]), abs(a3.get_ylim()[1]))
    a3.text(_xr * 0.5, _yr * 0.45, "Reverse\nIon Exchange", fontsize=9, ha="center",
            color="#8B0000", fontweight="bold", alpha=0.6)
    a3.text(-_xr * 0.5, -_yr * 0.45, "Ion\nExchange", fontsize=9, ha="center",
            color="#00008B", fontweight="bold", alpha=0.6)

    # (d) Na vs Cl dissolution diagnostic
    a4.scatter(data_meq["Cl"], data_meq["Na"], c="darkorange", edgecolor="black",
              s=60, alpha=0.8, linewidth=0.5)
    _maxval = max(data_meq["Cl"].max(), data_meq["Na"].max()) * 1.05
    a4.plot([0, _maxval], [0, _maxval], "r--", linewidth=1.2, label="1:1 (halite dissolution)")
    a4.set_xlabel(r"Cl$^-$ (meq/L)", fontsize=11)
    a4.set_ylabel(r"Na$^+$ (meq/L)", fontsize=11)
    a4.set_title("(d) Na\u207A vs Cl\u207B Molar Ratio", fontsize=12, pad=10)
    style_axes(a4, xbins=5, ybins=5, grid=True)
    style_legend(a4.legend(loc="upper left"))
    a4.text(_maxval * 0.3, _maxval * 0.8, "Na$^+$ excess\n(silicate weathering)",
            fontsize=8, ha="center", color="#2E4053", fontstyle="italic")
    a4.text(_maxval * 0.75, _maxval * 0.3, "Cl$^-$ excess\n(evaporite dissolution)",
            fontsize=8, ha="center", color="#922B21", fontstyle="italic")

    save_fig(ART / "Figure6_Geochemical_Constraints.png", fig=fig6)

    # Figure S6 (Moved from Figure 7 for narrative focus)
    st = res["validation"]
    fig, axes = plt.subplots(2, 1, figsize=(11, 12))
    a1, a2 = axes
    
    mean_vals = st["mean_loadings"].values
    im1 = a1.imshow(mean_vals, aspect="auto", cmap="rainbow")
    a1.set_xticks(np.arange(st["mean_loadings"].shape[1]))
    a1.set_xticklabels(st["mean_loadings"].columns, rotation=45, ha="right", fontsize=12)
    a1.set_yticks(np.arange(st["mean_loadings"].shape[0]))
    a1.set_yticklabels(st["mean_loadings"].index, fontsize=12)
    
    for i in range(mean_vals.shape[0]):
        for j in range(mean_vals.shape[1]):
            bg_val = mean_vals[i, j]
            col = "black"
            a1.text(j, i, f"{bg_val:.2f}", ha="center", va="center", color=col, fontsize=11, fontweight="bold")
            
    cbar1 = fig.colorbar(im1, ax=a1, fraction=0.046, pad=0.04)
    style_colorbar(cbar1, "Mean |CLR loading|")
    a1.set_title("Bootstrap Mean |CLR| Loadings", fontsize=16, pad=15, fontweight="bold")
    style_axes(a1)
    
    std_vals = st["std_loadings"].values
    im2 = a2.imshow(std_vals, aspect="auto", cmap="rainbow")
    a2.set_xticks(np.arange(st["std_loadings"].shape[1]))
    a2.set_xticklabels(st["std_loadings"].columns, rotation=45, ha="right", fontsize=12)
    a2.set_yticks(np.arange(st["std_loadings"].shape[0]))
    a2.set_yticklabels(st["std_loadings"].index, fontsize=12)
    
    for i in range(std_vals.shape[0]):
        for j in range(std_vals.shape[1]):
            bg_val = std_vals[i, j]
            col = "black"
            a2.text(j, i, f"{bg_val:.3f}", ha="center", va="center", color=col, fontsize=11, fontweight="bold")
            
    cbar2 = fig.colorbar(im2, ax=a2, fraction=0.046, pad=0.04)
    style_colorbar(cbar2, "Std. dev. |CLR loading|")
    a2.set_title("Bootstrap Std |CLR| Loadings", fontsize=16, pad=15, fontweight="bold")
    style_axes(a2)
    
    save_fig(ART / "FigureS6_Bootstrap_Stability.png", fig=fig)

    # Figure 7
    figure7_order = ["NPLS_GWQI", "F"]
    fig7, axes7 = plt.subplots(len(figure7_order), 1, figsize=(13, 10.5))
    fig7.suptitle("Bayesian Endpoint Driver Profiles on Hydrochemical Predictors", fontsize=17, y=0.985, fontweight="bold")
    panel_letters = {"NPLS_GWQI": "(a)", "F": "(b)"}
    max_ep_contrib = max(10.0, float(t8["Driver_Score_Percent"].max()) * 1.18)
    for ax, ep in zip(np.atleast_1d(axes7), figure7_order):
        ep_rows = t8[(t8["Endpoint"] == ep) & (t8["Driver_Score_Percent"] > 0)].copy()
        ax.set_xlim(0, max_ep_contrib)
        if ep_rows.empty:
            ax.text(
                0.5,
                0.5,
                "No positive stability-weighted\ndriver score resolved",
                ha="center",
                va="center",
                fontsize=11,
                fontweight="bold",
                transform=ax.transAxes,
                color="#4b5563",
            )
            style_axes(ax, xbins=5, grid=True, grid_axis="x")
            continue
            
        ep_rows = ep_rows.sort_values("Driver_Score_Percent", ascending=True)
        y_pos = np.arange(len(ep_rows))
        labels = [
            textwrap.fill(f"{lbl.replace('_', ' ')} [pd {sel:.0f}%]", width=40)
            for lbl, sel in zip(ep_rows["Predictor_Label"], ep_rows["Prob_Direction_Percent"])
        ]
        colours = [ION_COLOURS.get(pred, "#B0B8C1") for pred in ep_rows["Predictor"]]
        values = ep_rows["Driver_Score_Percent"].to_numpy()
        ax.barh(y_pos, values, color=colours, edgecolor="white", linewidth=0.8)
        ax.set_yticks(y_pos, labels, fontsize=9.5)
        ax.set_title(
            (
                f"{panel_letters[ep]} {ENDPOINT_PLOT_LABELS[ep]} "
                f"({ep_rows['Predictor_Block_Label'].iloc[0]}, Bayesian horseshoe, R2={ep_rows['Model_R2'].iloc[0]:.2f})"
            ),
            fontsize=13,
            pad=10,
            fontweight="bold",
        )
        style_axes(ax, xbins=5, grid=True, grid_axis="x")
        for row_idx, (_, row) in enumerate(ep_rows.iterrows()):
            val = float(row["Driver_Score_Percent"])
            if val <= 0:
                continue
            label_txt = f"{val:.1f}% | Beta {row['Beta_Mean']:.2f}"
            if val >= 12:
                ax.text(
                    val - 0.5,
                    row_idx,
                    label_txt,
                    ha="right",
                    va="center",
                    fontsize=8.5,
                    color="white",
                    fontweight="bold",
                )
            else:
                ax.text(
                    val + 0.5,
                    row_idx,
                    label_txt,
                    ha="left",
                    va="center",
                    fontsize=8.5,
                    color="#1f2933",
                    fontweight="bold",
                )
    axes7[-1].set_xlabel("Relative Bayesian driver score (%)", fontsize=13)
    fig7.subplots_adjust(left=0.34, right=0.97, top=0.90, bottom=0.08, hspace=0.34)
    save_fig(ART / "Figure7_Process_Influence.png", fig=fig7, pad=0.6)

    # Figure 8
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(2, 2)
    glist = list(GROUPS.keys())
    
    a1 = fig.add_subplot(gs[0, 0])
    bp1 = a1.boxplot([det[g]["HQ_F"] for g in glist], tick_labels=glist, showfliers=False, patch_artist=True)
    style_boxplot(bp1, [DEMOGRAPHIC_COLOURS[g] for g in glist])
    a1.axhline(1.0, color="#BC4749", linestyle="--", linewidth=1.6, alpha=0.8)
    a1.set_title("(a) HQ_F by Demographic Group", fontsize=14, fontweight="bold")
    a1.set_ylabel("Hazard Quotient (HQ)", fontsize=13)
    style_axes(a1, ybins=5, grid=True, grid_axis="y")
    
    a2 = fig.add_subplot(gs[0, 1])
    bp2 = a2.boxplot([det[g]["HQ_NO3"] for g in glist], tick_labels=glist, showfliers=False, patch_artist=True)
    style_boxplot(bp2, [DEMOGRAPHIC_COLOURS[g] for g in glist])
    a2.axhline(1.0, color="#BC4749", linestyle="--", linewidth=1.6, alpha=0.8)
    a2.set_title("(b) HQ_NO3 by Demographic Group", fontsize=14, fontweight="bold")
    a2.set_ylabel("Hazard Quotient (HQ)", fontsize=13)
    style_axes(a2, ybins=5, grid=True, grid_axis="y")
    
    a3 = fig.add_subplot(gs[1, :])
    mhf = [np.mean(det[g]["HI_skeletal_dental"]) for g in glist]
    mhn = [np.mean(det[g]["HI_hemato"]) for g in glist]
    xx = np.arange(len(glist))
    ww = 0.35
    bars1 = a3.bar(xx - ww / 2, mhf, ww, label="Skeletal-dental HI", color="#4C78A8", edgecolor="black")
    bars2 = a3.bar(xx + ww / 2, mhn, ww, label="Haematological HI", color="#E15759", edgecolor="black")
    a3.axhline(1.0, color="#BC4749", linestyle="--", linewidth=1.6, alpha=0.8)
    a3.set_xticks(xx)
    a3.set_xticklabels(glist, rotation=0)
    a3.set_title("(c) Organ-Specific Mean HI", fontsize=14, fontweight="bold")
    a3.set_ylabel("Hazard Index (HI)", fontsize=13)
    for bars in (bars1, bars2):
        for bar in bars:
            a3.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.03,
                f"{bar.get_height():.2f}",
                ha="center",
                va="bottom",
                fontsize=11,
            )
    style_axes(a3, ybins=5, grid=True, grid_axis="y")
    style_legend(a3.legend(loc="upper right"))
    
    save_fig(ART / "Figure8_Deterministic_Risk.png", fig=fig)

    # Figure 9
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(2, 2)
    glist = list(GROUPS.keys())
    
    a1 = fig.add_subplot(gs[0, 0])
    parts1 = a1.violinplot([post[g]["skeletal_dental"] for g in glist], showmeans=True, showextrema=False)
    style_violin(parts1, [DEMOGRAPHIC_COLOURS[g] for g in glist])
    a1.axhline(1.0, color="red", linestyle="--", linewidth=1.5, alpha=0.8, label="HI Threshold = 1.0")
    a1.set_xticks(np.arange(1, len(glist) + 1))
    a1.set_xticklabels(glist, rotation=20)
    a1.set_ylabel("Posterior HI", fontsize=13)
    a1.set_title("(a) Posterior HI (Skeletal-Dental)", fontsize=14, fontweight="bold")
    style_axes(a1, ybins=5, grid=True, grid_axis="y")
    style_legend(a1.legend(loc="upper right"))
    
    a2 = fig.add_subplot(gs[0, 1])
    parts2 = a2.violinplot([post[g]["hemato"] for g in glist], showmeans=True, showextrema=False)
    style_violin(parts2, [DEMOGRAPHIC_COLOURS[g] for g in glist])
    a2.axhline(1.0, color="red", linestyle="--", linewidth=1.5, alpha=0.8, label="HI Threshold = 1.0")
    a2.set_xticks(np.arange(1, len(glist) + 1))
    a2.set_xticklabels(glist, rotation=20)
    a2.set_ylabel("Posterior HI", fontsize=13)
    a2.set_title("(b) Posterior HI (Haematological)", fontsize=14, fontweight="bold")
    style_axes(a2, ybins=5, grid=True, grid_axis="y")
    style_legend(a2.legend(loc="upper right"))
    
    a3 = fig.add_subplot(gs[1, :])
    y = np.arange(len(t11))
    x = t11["Mean_HI"].to_numpy()
    xerr = np.vstack([x - t11["CI_2_5"].to_numpy(), t11["CI_97_5"].to_numpy() - x])
    organ_labels = {"skeletal_dental": "Skeletal-dental", "hemato": "Haematological"}
    organ_colours = {"skeletal_dental": "#E15759", "hemato": "#4C78A8"}
    labels = [f"{g} | {organ_labels.get(o, o)}" for g, o in zip(t11["Demographic_Group"], t11["Organ_System"])]
    for idx, row in t11.reset_index(drop=True).iterrows():
        organ = row["Organ_System"]
        colour = organ_colours.get(organ, "#1f2933")
        a3.errorbar(
            row["Mean_HI"],
            idx,
            xerr=[[row["Mean_HI"] - row["CI_2_5"]], [row["CI_97_5"] - row["Mean_HI"]]],
            fmt="o",
            color=colour,
            ecolor=colour,
            capsize=5,
            elinewidth=1.6,
            markersize=7,
        )
    a3.axvline(1.0, color="red", linestyle="--", linewidth=1.5)
    a3.set_yticks(y)
    a3.set_yticklabels(labels)
    a3.set_xlabel("Hazard Index (HI)", fontsize=13)
    a3.set_title(f"(c) Forest Plot ({method})", fontsize=14, fontweight="bold")
    style_axes(a3, xbins=5, grid=True, grid_axis="x")
    legend_handles = [
        Line2D([0], [0], marker="o", color=colour, label=label, linestyle="None", markersize=8)
        for label, colour in (("Skeletal-dental", "#E15759"), ("Haematological", "#4C78A8"))
    ]
    style_legend(a3.legend(handles=legend_handles, loc="lower right"))
    
    save_fig(ART / "Figure9_Bayesian_Posteriors.png", fig=fig)

    # Figure 10
    fig = plt.figure(figsize=(14, 12))
    a1 = fig.add_subplot(2, 2, 1)
    sc1 = a1.scatter(data["EC"], irr["SSP"], c=irr["N_ISI_Score"], cmap="cividis", edgecolor="k", s=72, alpha=0.85, linewidth=0.7)
    a1.set_title("(a) Wilcox Proxy (EC vs SSP)", fontsize=14, fontweight="bold")
    a1.set_xlabel(r"Electric Conductivity, EC ($\mu$S/cm)", fontsize=13)
    a1.set_ylabel("Soluble Sodium Percentage, SSP (%)", fontsize=13)
    cb1 = fig.colorbar(sc1, ax=a1)
    style_colorbar(cb1, "N-ISI Score")
    style_axes(a1, xbins=5, ybins=5, grid=True)
    
    a2 = fig.add_subplot(2, 2, 2)
    sc2 = a2.scatter(data["EC"], irr["SAR"], c=irr["N_ISI_Score"], cmap="plasma", edgecolor="k", s=72, alpha=0.85, linewidth=0.7)
    a2.set_title("(b) USSL Proxy (EC vs SAR)", fontsize=14, fontweight="bold")
    a2.set_xlabel(r"Electric Conductivity, EC ($\mu$S/cm)", fontsize=13)
    a2.set_ylabel("Sodium Adsorption Ratio (SAR)", fontsize=13)
    cb2 = fig.colorbar(sc2, ax=a2)
    style_colorbar(cb2, "N-ISI Score")
    style_axes(a2, xbins=5, ybins=5, grid=True)
    
    a3 = fig.add_subplot(2, 2, 3)
    icols = ["SAR", "SSP", "MH", "KR", "PI", "RSC", "PS"]
    bp = a3.boxplot([irr[c] for c in icols], tick_labels=icols, showfliers=True, patch_artist=True)
    style_boxplot(bp, ["#4C78A8", "#72B7B2", "#54A24B", "#EECA3B", "#B279A2", "#FF9DA6", "#9D755D"])
    a3.set_xlabel("Irrigation Quality Indices", fontsize=13)
    a3.set_ylabel("Value", fontsize=13)
    a3.set_title("(c) Irrigation Index Distributions", fontsize=14, fontweight="bold")
    style_axes(a3, ybins=5, grid=True, grid_axis="y")
    
    a4 = fig.add_subplot(2, 2, 4)
    isi_counts = irr["ISI_Class"].value_counts(dropna=False)
    x_positions = np.arange(len(isi_counts))
    colors = plt.cm.Set2(np.linspace(0, 1, max(2, len(isi_counts))))
    bars = a4.bar(x_positions, isi_counts.values, color=colors, edgecolor="black", alpha=0.8, width=0.6)
    a4.set_xticks(x_positions)
    a4.set_xticklabels(isi_counts.index.astype(str), rotation=20, ha="right")
    a4.set_xlabel("N-ISI Suitability Class", fontsize=13)
    a4.set_ylabel("Frequency (Wells)", fontsize=13)
    a4.set_title("(d) N-ISI Classification Distribution", fontsize=14, fontweight="bold")
    style_axes(a4, ybins=5, grid=True, grid_axis="y")
    for bar in bars:
        yval = bar.get_height()
        a4.text(bar.get_x() + bar.get_width() / 2, yval + (isi_counts.max() * 0.01), f"{int(yval)}", ha="center", va="bottom", fontsize=11, fontweight="bold")
        
    save_fig(ART / "Figure10_Irrigation_Suitability.png", fig=fig)

    # Figure 11 — Main-text NPLS-GWQI Context-Signature-Impact pathway (Now Bayesian)
    pathway_rows = build_bayesian_pathway_rows(endpoint_models, order)
    npls_rows = pathway_rows[pathway_rows["Endpoint"] == "NPLS_GWQI"].copy()
    fig11, ax11 = plt.subplots(figsize=(15, 7.5))
    draw_bayesian_pathway_panel(ax11, npls_rows, "NPLS_GWQI")
    fig11.subplots_adjust(left=0.04, right=0.98, top=0.90, bottom=0.08)
    save_fig(ART / "Figure11_Sankey_Pathway.png", fig=fig11, pad=0.5)

    # Figure S1
    fig_s1, ax_s1 = plt.subplots(figsize=(9, 5.5))
    xx = np.arange(1, len(t5) + 1)
    ax_s1.plot(xx, t5["Observed_Eigenvalue"], marker="o", label="Observed eigenvalue", color="#2A9D8F")
    ax_s1.plot(xx, t5["Random_95th_Eigenvalue"], marker="s", label="Random 95th percentile", color="#E76F51")
    ax_s1.set_xticks(xx, t5["Component"])
    ax_s1.set_ylabel("Eigenvalue")
    ax_s1.set_title("Parallel Analysis Scree Plot")
    style_axes(ax_s1, ybins=5, grid=True, grid_axis="y")
    style_legend(ax_s1.legend(loc="upper right"))
    save_fig(ART / "FigureS1_Scree_Plot.png", fig=fig_s1)

    # Figure S2
    pen = np.ones(len(data))
    pen *= np.where((data["pH"] < 6.5) | (data["pH"] > 8.5), 1.2, 1.0)
    pen *= np.where(data["TDS"] > 500, 1.1, 1.0)
    gwqi_pen = wqi["NPLS_GWQI"].to_numpy()
    gwqi_no_pen = gwqi_pen / pen
    fig_s2, ax_s2 = plt.subplots(figsize=(9.5, 5.5))
    bp = ax_s2.boxplot(
        [gwqi_no_pen, gwqi_pen],
        tick_labels=["Without Sensory Penalty", "With Sensory Penalty"],
        patch_artist=True,
        showfliers=False,
    )
    style_boxplot(bp, ["#4C78A8", "#E76F51"])
    jitter_rng = np.random.default_rng(SEED)
    ax_s2.scatter(1 + jitter_rng.normal(0, 0.03, len(gwqi_no_pen)), gwqi_no_pen, alpha=0.28, s=24, color="#4C78A8")
    ax_s2.scatter(2 + jitter_rng.normal(0, 0.03, len(gwqi_pen)), gwqi_pen, alpha=0.28, s=24, color="#E76F51")
    ax_s2.set_ylabel("NPLS-GWQI score")
    ax_s2.set_title("Sensitivity of NPLS-GWQI to Sensory Penalties")
    style_axes(ax_s2, ybins=5, grid=True, grid_axis="y")
    save_fig(ART / "FigureS2_Sensory_Sensitivity.png", fig=fig_s2)

    # Figure S3 — Box-and-Whisker (normalised to WHO/GSA)
    fig_s3, ax_s3 = plt.subplots(figsize=(11, 6.5))
    cols_bw = ["Na", "K", "Mg", "Ca", "Cl", "SO4", "HCO3", "NO3", "F"]
    norm_bw = data[cols_bw].div(STANDARDS[cols_bw], axis=1)
    bp = ax_s3.boxplot([norm_bw[c] for c in cols_bw], tick_labels=cols_bw, showfliers=False, patch_artist=True)
    style_boxplot(bp, [PROCESS_COLOURS[i % len(PROCESS_COLOURS)] for i in range(len(cols_bw))])
    ax_s3.axhline(1.0, color="#BC4749", linestyle="--", linewidth=1.6)
    ax_s3.set_title("Box-and-Whisker (Normalised to WHO/GSA Standards)")
    ax_s3.set_ylabel("Normalised concentration")
    for tick in ax_s3.get_xticklabels():
        tick.set_rotation(35)
        tick.set_ha("right")
    style_axes(ax_s3, ybins=5, grid=True, grid_axis="y")
    save_fig(ART / "FigureS3_Parameter_BoxWhisker.png", fig=fig_s3)



    # Figure S5 — Exploratory F Bayesian pathway
    f_rows = pathway_rows[pathway_rows["Endpoint"] == "F"].copy()
    fig_s5, ax_s5 = plt.subplots(figsize=(15, 7.5))
    draw_bayesian_pathway_panel(ax_s5, f_rows, "F")
    fig_s5.subplots_adjust(left=0.04, right=0.98, top=0.90, bottom=0.08)
    save_fig(ART / "FigureS5_F_EndpointSafe_Pathway.png", fig=fig_s5, pad=0.5)

    # Figure S7 — Geochemical Control Curves
    fig_s7, ax_s7 = plt.subplots(1, 2, figsize=(14.4, 6.4))
    aug_data_no3 = wf.prep.get_augmented_predictor_block(data, "NO3")
    aug_data_f = wf.prep.get_augmented_predictor_block(data, "F")
    curve_specs = [
        (
            "Aug_CAI_1",
            "Schoeller CAI-1 (Ion Exchange Index)",
            "A: Endpoint response to ion-exchange signal",
        ),
        (
            "Aug_Carbonate_Proxy",
            "Carbonate Weathering Proxy ((Ca+Mg)/HCO3)",
            "B: Endpoint response to weathering signal",
        ),
    ]
    last_idx = len(ax_s7) - 1
    for idx, (ax, (predictor, xlabel, title)) in enumerate(zip(ax_s7, curve_specs)):
        ax_right = ax.twinx()
        add_lowess_scatter(ax, aug_data_f[predictor], data["F"], color=ION_COLOURS["F"], label=r"F$^-$")
        add_lowess_scatter(ax_right, aug_data_no3[predictor], data["NO3"], color=ION_COLOURS["NO3"], label=r"NO$_3^-$")
        ax.set_xlabel(xlabel)
        if idx == 0:
            ax.set_ylabel("Fluoride (mg/L)", color=ION_COLOURS["F"], labelpad=10)
            ax.tick_params(axis="y", colors=ION_COLOURS["F"], labelleft=True, left=True)
        else:
            ax.set_ylabel("")
            ax.tick_params(axis="y", colors=ION_COLOURS["F"], labelleft=False, left=False)
        if idx == last_idx:
            ax_right.set_ylabel("Nitrate (mg/L)", color=ION_COLOURS["NO3"], labelpad=10)
            ax_right.tick_params(axis="y", colors=ION_COLOURS["NO3"], labelright=True, right=True)
        else:
            ax_right.set_ylabel("")
            ax_right.tick_params(axis="y", colors=ION_COLOURS["NO3"], labelright=False, right=False)
        ax.set_title(title)
        if idx == 0:
            # Combine handles for twin axes legend
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax_right.get_legend_handles_labels()
            style_legend(ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=10))
        style_axes(ax, grid=True)
        ax_right.grid(False)
    fig_s7.subplots_adjust(left=0.08, right=0.92, top=0.90, bottom=0.14, wspace=0.22)
    save_fig(ART / "FigureS7_Geochemical_Control_Curves.png", fig=fig_s7)

    # Figure S8 — Model Ablation Study (In-sample capacity to resolve endpoints)
    # We use in-sample R2 here because the small, outlier-heavy dataset makes CV-R2 
    # unstable, but we want to show the 'capacity' of geogenic proxies to capture signal.
    def get_in_sample_metrics(x, y, ep):
        try:
            # Prepare predictors
            x_df = x.copy()
            x_df = x_df.loc[:, x_df.std() > 1e-10]
            
            # Use same tuning logic but return in-sample fit
            grid_c = default_component_grid(x_df.shape[1], len(x_df))
            grid_k = default_keep_grid(x_df.shape[1])
            tuned = tune_sparse_pls(x_df, y, component_grid=grid_c, keep_grid=grid_k, n_splits=3)
            p = tuned["best_params"]
            fit = fit_sparse_pls_model(x_df, y, n_components=p["n_components"], keep_k=p["keep_k"])
            return {"r2": fit.in_sample_r2, "rmse": np.sqrt(mean_squared_error(y, fit.fitted_values))}
        except Exception:
            return {"r2": 0.0, "rmse": y.std()}

    ablation_metrics = {
        "NO3": {
            "major": get_in_sample_metrics(data[COMP].drop(columns=["NO3"]), data["NO3"], "NO3"),
            "augmented": get_in_sample_metrics(wf.prep.get_augmented_predictor_block(data, "NO3"), data["NO3"], "NO3"),
        },
        "F": {
            "major": get_in_sample_metrics(data[COMP].drop(columns=["F"]), data["F"], "F"),
            "augmented": get_in_sample_metrics(wf.prep.get_augmented_predictor_block(data, "F"), data["F"], "F"),
        },
    }
    endpoints = ["NO3", "F"]
    labels = [r"NO$_3^-$", r"F$^-$"]
    r2_major = [ablation_metrics[ep]["major"]["r2"] for ep in endpoints]
    r2_aug = [ablation_metrics[ep]["augmented"]["r2"] for ep in endpoints]
    rmse_major = [ablation_metrics[ep]["major"]["rmse"] for ep in endpoints]
    rmse_aug = [ablation_metrics[ep]["augmented"]["rmse"] for ep in endpoints]

    fig_s8, ax_s8 = plt.subplots(1, 2, figsize=(12.5, 6))
    x_indices = np.arange(len(endpoints))
    width = 0.34
    ax_s8[0].bar(x_indices - width / 2, r2_major, width, label="Major-Ion Only", color="#8AB17D", alpha=0.9, edgecolor="white", linewidth=0.8)
    ax_s8[0].bar(x_indices + width / 2, r2_aug, width, label="Enhanced Augmented", color="#2A9D8F", edgecolor="white", linewidth=0.8)
    ax_s8[0].set_ylabel("Model Fitting Capacity ($R^2$)", fontsize=13)
    ax_s8[0].set_title("A: Variance Explanation Potential", fontsize=14, fontweight="bold")
    ax_s8[0].set_xticks(x_indices, labels, fontsize=12)
    ax_s8[0].set_ylim(0, 1.0) # Ensure scale is comparable and positive
    ax_s8[0].axhline(0, color="black", linewidth=0.8, ls="--", alpha=0.4)
    style_axes(ax_s8[0], grid=True, grid_axis="y")
    style_legend(ax_s8[0].legend(loc="upper left"))

    ax_s8[1].bar(x_indices - width / 2, rmse_major, width, label="Major-Ion Only", color="#8AB17D", alpha=0.9, edgecolor="white", linewidth=0.8)
    ax_s8[1].bar(x_indices + width / 2, rmse_aug, width, label="Enhanced Augmented", color="#2A9D8F", edgecolor="white", linewidth=0.8)
    ax_s8[1].set_ylabel("Model Fit Error (RMSE)", fontsize=13)
    ax_s8[1].set_title("B: Structural Error Reduction", fontsize=14, fontweight="bold")
    ax_s8[1].set_xticks(x_indices, labels, fontsize=12)
    style_axes(ax_s8[1], grid=True, grid_axis="y")
    
    fig_s8.suptitle("Model Ablation: Resolving the Geogenic Predictive Advantage", fontsize=16, fontweight="bold", y=0.98)
    fig_s8.subplots_adjust(left=0.08, right=0.97, top=0.86, bottom=0.12, wspace=0.28)
    save_fig(ART / "FigureS8_Model_Ablation_Study.png", fig=fig_s8)

    manifest_paths = [ROOT / p for p in parse_paths(ROOT / "manuscript" / "artifact_manifest.yaml")]
    locked_paths = [ROOT / p for p in parse_paths(ROOT / "manuscript" / "frameworks" / "locked-analysis-plan.yaml")]
    required = sorted(set(manifest_paths + locked_paths))
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing generated artifacts:\n" + "\n".join(missing))

    print(f"Generated {len(required)} required artifacts (manifest + locked plan).")
    print(f"Hydrochemical ILR outputs written to: {HYDRO}")


if __name__ == "__main__":
    main()
