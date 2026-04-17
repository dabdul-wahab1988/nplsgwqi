# Data Dictionary (replication metadata)

This data dictionary describes the input variables in `data.csv` used to generate the manuscript analyses and artifacts.

## Input table: `data.csv`

All concentration variables are reported in mg/L unless stated otherwise.

| Variable | Units | Description |
|---|---:|---|
| `SampleID` | — | Borehole/sample identifier (string key). |
| `pH` | — | Field pH (unitless). |
| `EC` | µS/cm | Electrical conductivity (field). |
| `TDS` | mg/L | Total dissolved solids (field/lab as provided). |
| `Na` | mg/L | Sodium concentration. |
| `K` | mg/L | Potassium concentration. |
| `Mg` | mg/L | Magnesium concentration. |
| `Ca` | mg/L | Calcium concentration. |
| `Cl` | mg/L | Chloride concentration. |
| `SO4` | mg/L | Sulphate concentration (as SO₄²⁻). |
| `HCO3` | mg/L | Bicarbonate concentration (as HCO₃⁻ / alkalinity proxy). |
| `NO3` | mg/L | Nitrate concentration (as NO₃⁻). |
| `F` | mg/L | Fluoride concentration (as F⁻). |

## Derived outputs (computed from `data.csv`)

Derived tables and figures are written to `manuscript/artifacts/`. Key derived quantities include:

- Major-ion conversion to meq/L for charge-balance auditing and irrigation indices.
- Charge-balance error (CBE, %) per sample and summary statistics.
- ILR-PCA latent-process loadings (reported as CLR loading profiles).
- Bayesian endpoint-driver model diagnostics (R-hat/ESS, PPC, PSIS-LOO/WAIC) and prior/basis sensitivity summaries.
- PHREEQC saturation indices (SI) for calcite/dolomite/fluorite under a temperature band.
- Benchmarking outputs comparing NPLS-GWQI to conventional and fuzzy-only WQI baselines.
