import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from nplsgwqi.figure_style import apply_publication_style, save_figure, style_axes, style_legend
from nplsgwqi.wqi import calculate_npls_gwqi

ROOT = Path(__file__).resolve().parent
ART = ROOT / "manuscript" / "artifacts"

apply_publication_style()

# Standard limits from generation_script
STANDARDS = pd.Series({
    "pH": 8.5, "EC": 1500.0, "TDS": 1000.0, "Na": 200.0, "K": 12.0,
    "Mg": 50.0, "Ca": 75.0, "Cl": 250.0, "SO4": 250.0, "HCO3": 300.0,
    "NO3": 50.0, "F": 1.5,
})

def run_sobol_wqi():
    """Table S4: Sobol Sensitivity Analysis (Global Variance Decomposition)"""
    N_SAMPLES = 4096
    print(f"Starting Global Sobol Sensitivity Analysis (N={N_SAMPLES})...")
    
    # Define problem space based on standards
    names = STANDARDS.index.tolist()
    bounds = []
    for p in names:
        if p == 'pH':
            bounds.append([4.0, 10.0]) # V-shaped hazard range
        else:
            std = STANDARDS[p]
            bounds.append([std * 0.1, std * 2.5]) # Range from safe to bad
            
    problem = {
        'num_vars': len(names),
        'names': names,
        'bounds': bounds
    }

    # Generate Samples
    param_values = sobol_sample.sample(problem, N_SAMPLES)
    
    # Run the model
    df_samples = pd.DataFrame(param_values, columns=names)
    
    # Evaluate WQI
    # We use n_bootstraps=1 for speed; Sobol evaluates model structure sensitivity
    wqi_res, _ = calculate_npls_gwqi(df_samples, STANDARDS, n_bootstraps=1)
    Y = wqi_res['NPLS_GWQI'].values
    
    # Analyze
    Si = sobol.analyze(problem, Y, print_to_console=False)
    
    # Create Table S4
    t_s4 = pd.DataFrame({
        'Parameter': names,
        'S1': Si['S1'],
        'ST': Si['ST'],
        'S1_conf': Si['S1_conf'],
        'ST_conf': Si['ST_conf'],
        'Interaction_Index': Si['ST'] - Si['S1']
    }).sort_values('ST', ascending=False)
    
    t_s4.to_csv(ART / "TableS4_Sobol_Sensitivity.csv", index=False)
    print("Generated Table S4 (Sobol Indices).")
    
    # Figure S8: Sobol Sensitivity Indices Bar Chart with Error Bars
    fig, ax = plt.subplots(figsize=(13, 7))
    x = np.arange(len(t_s4))
    width = 0.35
    
    ax.bar(
        x - width / 2,
        t_s4["S1"],
        width,
        yerr=t_s4["S1_conf"],
        label="First-order (S1)",
        color="#4C78A8",
        alpha=0.9,
        capsize=4,
        edgecolor="white",
        linewidth=0.8,
    )
    ax.bar(
        x + width / 2,
        t_s4["ST"],
        width,
        yerr=t_s4["ST_conf"],
        label="Total-order (ST)",
        color="#1D3557",
        alpha=0.9,
        capsize=4,
        edgecolor="white",
        linewidth=0.8,
    )
    
    ax.set_xticks(x, t_s4["Parameter"], rotation=30, ha="right")
    ax.set_ylabel("Sensitivity index (variance fraction)")
    ax.set_title(f"Global Sensitivity Analysis (NPLS-GWQI Structure, N={N_SAMPLES})")
    style_axes(ax, ybins=5, grid=True, grid_axis="y")
    style_legend(ax.legend(loc="upper right"))
    save_figure(ART / "FigureS8_Sobol_Sensitivity.png", fig=fig, dpi=600)
    print("Generated Figure S8 (Interaction Map with Error Bars).")

def run_sobol_risk():
    """Table S5: Sobol Sensitivity Analysis for Health Risk (Children)"""
    N_SAMPLES = 4096
    print(f"Starting Global Sobol Sensitivity Analysis for Health Risk (N={N_SAMPLES})...")
    
    names = ['C_F', 'C_NO3', 'BW', 'IR']
    problem = {
        'num_vars': 4,
        'names': names,
        'bounds': [
            [0.1, 5.0],   # C_F
            [1.0, 100.0], # C_NO3
            [10.0, 20.0], # BW
            [0.5, 2.0]    # IR
        ]
    }

    param_values = sobol_sample.sample(problem, N_SAMPLES)
    
    # Evaluate HI for Children
    RfD_F = 0.06
    RfD_NO3 = 1.6
    
    C_F = param_values[:, 0]
    C_NO3 = param_values[:, 1]
    BW = param_values[:, 2]
    IR = param_values[:, 3]
    
    Y = ((C_F * IR) / (BW * RfD_F)) + ((C_NO3 * IR) / (BW * RfD_NO3))
    
    # Analyze
    Si = sobol.analyze(problem, Y, print_to_console=False)
    
    # Create Table S5
    t_s5 = pd.DataFrame({
        'Parameter': names,
        'S1': Si['S1'],
        'ST': Si['ST'],
        'S1_conf': Si['S1_conf'],
        'ST_conf': Si['ST_conf'],
        'Interaction_Index': Si['ST'] - Si['S1']
    }).sort_values('ST', ascending=False)
    
    t_s5.to_csv(ART / "TableS5_Risk_Sobol_Sensitivity.csv", index=False)
    print("Generated Table S5 (Risk Sobol Indices).")
    
    # Figure S9: Rigorous Variance Decomposition Pie Chart
    # Use S1 for direct contributions and (1 - sum(S1)) for interactions/residual
    s1_sum = t_s5['S1'].sum()
    interaction_resid = max(0, 1 - s1_sum)
    
    pie_data = t_s5['S1'].tolist() + [interaction_resid]
    pie_labels = t_s5['Parameter'].tolist() + ['Higher-Order Interactions']
    colors = ["#4C78A8", "#F28E2B", "#59A14F", "#E15759", "#B07AA1"]

    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    wedges, texts, autotexts = ax.pie(
        pie_data,
        labels=pie_labels,
        autopct="%1.1f%%",
        startangle=140,
        colors=colors[: len(pie_labels)],
        wedgeprops={"edgecolor": "white", "linewidth": 1.5},
        textprops={"fontsize": 12},
        pctdistance=0.78,
    )

    centre_circle = plt.Circle((0, 0), 0.68, fc="white")
    ax.add_artist(centre_circle)
    for autotext in autotexts:
        autotext.set_fontsize(11)
        autotext.set_fontweight("bold")
    ax.set_title("Variance Decomposition of Health Risk Drivers (Children)\n(S1 Individual Contribution vs. Interactions)")
    save_figure(ART / "FigureS9_Risk_Sobol_Pie.png", fig=fig, dpi=600)
    print("Generated Figure S9 (Refined Risk Sobol Visualization).")

if __name__ == "__main__":
    ART.mkdir(parents=True, exist_ok=True)
    run_sobol_wqi()
    run_sobol_risk()
