# NPLS-GWQI (Neutrosophic PLS Groundwater Quality Index)

Research-oriented Python code and manuscript sources for an integrated groundwater assessment workflow (quality indexing, compositional process discovery, Bayesian driver attribution, health risk, and irrigation suitability), implemented for a 34-sample Karaga District (Northern Ghana) case study.

This repository contains both:

- `src/nplsgwqi/`: the Python package code (no published PyPI package yet)
- `manuscript/`: manuscript sources, artifacts, and reviewer-response materials

## What’s inside

- **NPLS-GWQI**: `calculate_npls_gwqi()` computes an advanced Neutrosophic PLS-based Groundwater Quality Index using:
  - standards-based normalization (with special bipolar pH handling),
  - neutrosophic channels (Truth/Indeterminacy/Falsity),
  - VIP-style dynamic weights from a PLSRegression fit,
  - bootstrapping for weight stability,
  - optional sensory penalties (e.g., TDS, Fe).
- **Process attribution**: `run_process_contribution_model()` fits a constrained regression on process scores and returns:
  - relative contribution shares, and
  - grouped Shapley-style importance (process variability ranking).
- **Health risk screening**: `assess_health_risk()` provides deterministic (fast) and optional Bayesian (PyMC-based) risk estimates.
- **Irrigation indices**: `calculate_irrigation_indices()` computes common irrigation suitability indices (see module for details).
- **Integrated workflow**: `IntegratedWorkflow` (in `nplsgwqi.orchestrator`) ties together preprocessing, process discovery, regression, WQI, risk, irrigation, and validation.

## Manuscript + reviewer-response materials

Key files under `manuscript/`:

- `manuscript/sections/`: the manuscript section sources (`.md`)
- `manuscript/supplementary-methodology.md`: supplementary methodology (synced from `Supplementary Material.docx` with an addendum for reviewer-requested analyses)
- `manuscript/reviewer_comments.md`: verbatim reviewer-style comments extracted from `Comments on Karaga Manuscript.docx`
- `manuscript/response_to_comments.md`: comment-by-comment response with file+line anchors into the manuscript sources
- `manuscript/artifacts/`: tables/figures referenced by the manuscript
- `manuscript/analysis_workflow.md` and `manuscript/data_dictionary.md`: replication metadata (workflow diagram + variable dictionary)

## Requirements

Minimum: Python **3.10+**

Core dependencies (imported by the package):

- `numpy`
- `pandas`
- `scikit-learn`

Optional (only needed if you use Bayesian risk mode):

- `pymc`

For the reviewer-response artifact scripts in `manuscript/`, additional dependencies are typically required:

- `matplotlib`
- `arviz`
- `phreeqpython`

## Install / use locally

Because this repo doesn’t yet ship packaging metadata (`pyproject.toml`), the simplest way to run it is:

1. Clone the repository
2. Ensure the repository root is on your `PYTHONPATH`, or run your scripts from the repo root so `src/` is importable

Example (PowerShell):

```powershell
cd nplsgwqi
$env:PYTHONPATH = (Resolve-Path .\src).Path
python -c "from nplsgwqi import calculate_npls_gwqi; print(calculate_npls_gwqi)"
```

## Reproduce manuscript artifacts (local)

From the repo root:

- Baseline analysis artifacts: `python generation_script.py`
- Reviewer-requested additions (CBE tables, PHREEQC SI, benchmarking, irrigation disagreement, risk summaries): `python manuscript/reviewer_revision_artifacts.py`
- Regenerate the response-to-comments Markdown: `python manuscript/generate_response_to_comments.py`

Notes:

- Bayesian driver-model diagnostics can be compute-intensive depending on sampling settings (draws/tune/chains).
- Outputs are written under `manuscript/artifacts/` and should match the tables/figures referenced in `manuscript/sections/**/*.md`.

## Quickstart: compute NPLS-GWQI

`calculate_npls_gwqi(data, standards)` expects:

- `data`: `pandas.DataFrame` where columns are analytes (e.g., `pH`, `TDS`, `NO3`, `Fe`, …)
- `standards`: `pandas.Series` indexed by the same analyte names (guideline limits)

```python
import pandas as pd
from nplsgwqi import calculate_npls_gwqi

# Example dataset (replace with your measured concentrations)
data = pd.DataFrame(
    {
        "pH": [6.8, 7.3, 8.1],
        "TDS": [240, 610, 430],
        "NO3": [2.1, 15.3, 6.7],
        "Fe": [0.12, 0.45, 0.08],
    }
)

# Example standards (replace with your regulatory / guideline values)
standards = pd.Series(
    {
        "pH": 7.0,     # used via special bipolar normalization
        "TDS": 500.0,
        "NO3": 50.0,
        "Fe": 0.3,
    }
)

results, dynamic_weights = calculate_npls_gwqi(data=data, standards=standards, n_bootstraps=50)
print(results[["NPLS_GWQI", "Quality_Class"]])
print(dynamic_weights.sort_values(ascending=False))
```

Notes:

- Columns in `data` should match the index of `standards` (extra columns are ignored if not included in `standards`).
- Handle missing values before calling (e.g., impute/drop) to avoid model failures.

## Process contribution model

If you already have process/component scores (e.g., from a factor model / PCA-style workflow), you can estimate process contributions to an endpoint.

```python
import pandas as pd
from nplsgwqi import run_process_contribution_model

process_scores = pd.DataFrame(
    {
        "ProcessA": [0.2, 0.1, 0.7],
        "ProcessB": [0.5, 0.2, 0.4],
    }
)

endpoint = pd.Series([10.0, 11.5, 9.3], name="NO3")

out = run_process_contribution_model(process_scores=process_scores, target=endpoint, use_ridge=True)
print(out["relative_process_contributions"])
print(out["process_variability_importance"])
```

## Health risk screening (deterministic mode)

Deterministic mode is much faster and avoids MCMC sampling:

```python
import pandas as pd
from nplsgwqi import assess_health_risk

data = pd.DataFrame(
    {
        "F": [0.7, 1.2, 0.3],
        "NO3": [2.1, 15.3, 6.7],
        "Pb": [0.002, 0.006, 0.001],
        "As": [0.001, 0.003, 0.0005],
    }
)

risk = assess_health_risk(data, mode="deterministic")
print([k for k in risk.keys() if k.startswith("HI_sample_")][:5])
```

Bayesian mode (`mode="bayesian"` or `mode="both"`) requires `pymc` and can take time depending on `n_samples`, `tune`, and `chains`.

## Repository layout

```
manuscript/
  sections/
  artifacts/
  reviewer_comments.md
  response_to_comments.md
src/
  nplsgwqi/
    __init__.py
    wqi.py
    risk.py
    process_regression.py
    ...
tests/
data.csv
```

## Citation / attribution

If you use this code in academic work, please cite the associated manuscript/project where appropriate.

## Disclaimer

This code is provided for research and screening purposes. Always validate results against domain expertise, local standards, and quality assurance procedures before making operational or regulatory decisions.
