# A Dual-Space Compositional Process Framework for Groundwater Quality Assessment, Relative Hydrochemical Driver Influence, and Probabilistic Health Risk in Karaga District, Ghana

**Revised outline aligned with the current validated code path and artifact manifest (`data.csv`: 34 groundwater samples, 12 measured hydrochemical variables plus `SampleID`).**

---

## Data Availability Summary

| Variable | Available? | Current role in the codebase |
|---|---|---|
| `SampleID` | Yes | Identifier only |
| `pH`, `EC`, `TDS` | Yes | Field/non-compositional block; used directly in WQI and augmented endpoint models |
| `Na`, `K`, `Mg`, `Ca`, `Cl`, `SO4`, `HCO3` | Yes | Major-ion chemistry; used in meq/L conversion, compositional discovery, geochemical constraints, WQI, irrigation, and endpoint modelling |
| `NO3`, `F` | Yes | Included in the hydrochemical composition; also used explicitly in health-risk assessment; `F` is part of the locked Bayesian endpoint-driver analysis |
| `Fe`, `Pb`, `As`, `Cd`, `Ni` | No | Not available in the current dataset; trace-metal risk and Fe sensory pathways are inactive in practice |
| `CO3`, `PO4`, `temperature` | No | Not available; any related code path is inactive or defaults to absent-data handling |

---

## 1. Introduction (~800-1,000 words)

### 1.1 Groundwater dependence and water-quality relevance
* Importance of groundwater for domestic use and smallholder livelihoods in Karaga District.
* Ongoing concern about geogenic fluoride enrichment and broader hydrochemical variability in Voltaian Basin settings.
* Public-health and resource-management relevance of distinguishing drinking-water quality, health risk, and irrigation suitability.

### 1.2 Limits of conventional assessment approaches
* Weaknesses of conventional arithmetic GWQI approaches when threshold ambiguity and variable reliability are ignored.
* Why raw-space multivariate interpretation is problematic for closed hydrochemical compositions.
* Why deterministic HQ/HI values alone can understate uncertainty in population-level exposure interpretation.

### 1.3 Code-supported methodological contribution
* A dual-space workflow is implemented:
  * Process space: ILR-based compositional discovery plus deterministic geochemical auditing.
  * Assessment space: NPLS-GWQI, Bayesian hydrochemical driver profiling, deterministic and Bayesian health-risk assessment, and entropy-weighted neutrosophic irrigation suitability.
* The current locked artifact set emphasises:
  * hydrochemical facies and geochemical signatures,
  * groundwater-quality structure,
  * Bayesian driver profiles for `NPLS-GWQI` and `F`,
  * deterministic and Bayesian risk for `F` and `NO3`,
  * irrigation suitability through `N-ISI`.

### 1.4 Objectives
1. To characterise groundwater hydrochemistry and evaluate analytical integrity using charge-balance diagnostics and ILR-based facies visualisation.
2. To assess groundwater quality using a neutrosophic PLS-based groundwater quality index with bootstrapped dynamic weighting and aesthetic sensitivity handling.
3. To identify latent hydrochemical signatures using ILR-PCA with parallel analysis and to interpret them using deterministic geochemical constraints.
4. To profile the relative hydrochemical drivers of `NPLS-GWQI` and fluoride using Bayesian horseshoe regression on an augmented predictor block, with supplementary comparator diagnostics.
5. To estimate deterministic and Bayesian non-carcinogenic health risk for `F` and `NO3` across three demographic groups: Adults, Children, and Teens.
6. To evaluate irrigation suitability using seven standard indices combined through entropy-weighted neutrosophic irrigation suitability indexing (`N-ISI`).

---

## 2. Study Area (~500-700 words)

### 2.1 Location, administrative setting, and settlement pattern
* Geographic position of Karaga District in northern Ghana.
* Rural settlement pattern and borehole dependence.
* Relevance of the sampled communities and site spread shown in Figure 1.

### 2.2 Climate, hydrology, and groundwater reliance
* Semi-arid climatic context, recharge seasonality, and likely dry-season concentration effects.
* Dependence on groundwater as a year-round water source.

### 2.3 Geology and hydrogeological context
* Karaga District lies within the Voltaian sedimentary basin setting represented in Figure 1.
* The mapped district is dominated by the `Bimbila Formation`, with an undifferentiated `Obosum Group` mainly in the south-southwest and discontinuous `Bunya Formation` bodies in an east-central to southeastern belt plus smaller northwestern outliers.
* Mapped faults and lithostratigraphic contrasts provide a geological basis for spatial hydrochemical variability.

### 2.4 Land use and likely anthropogenic pressures
* Agriculture, settlement activities, and local sanitary conditions as plausible modifiers of nitrate and general water quality.
* Framing of anthropogenic pressure as a secondary modifier superimposed on a geogenic hydrochemical template.

---

## 3. Materials and Methods (~1,500-2,000 words)

### 3.1 Dataset, sampling, and measured variables
* `data.csv` contains 34 groundwater samples and 12 measured hydrochemical variables plus `SampleID`.
* Field variables: `pH`, `EC`, `TDS`.
* Major dissolved constituents: `Na`, `K`, `Mg`, `Ca`, `Cl`, `SO4`, `HCO3`, `NO3`, `F`.
* Figure 1 is an external GIS/cartographic study-area asset; the remaining analytical figures are code-generated.

### 3.2 Preprocessing and stoichiometric conversion (`preprocessing.py`)
* Conversion of relevant ions from mg/L to meq/L using valence-based factors.
* Charge-balance error calculation from summed cations and anions.
* Zero/non-positive replacement using positive-floor or lognormal-style substitution where needed.
* Variable partitioning in the current workflow:
  * compositional hydrochemical block: `Ca`, `Mg`, `Na`, `K`, `HCO3`, `Cl`, `SO4`, `NO3`, `F`;
  * field/non-compositional block: `pH`, `EC`, `TDS`;
  * explicit risk endpoints: `NO3`, `F`.
* For endpoint-driver modelling, the code constructs an augmented predictor block that combines log-transformed ions, field variables, and derived indices/proxies such as `CAI_1`, `CAI_2`, `Na_Cl`, `Ca_HCO3`, `Mg_Ca`, `Silicate_Proxy`, and `Carbonate_Proxy`.

### 3.3 NPLS-GWQI (`wqi.py`)
* Standard-based normalisation of the measured water-quality variables.
* Bipolar hazard treatment for `pH` rather than a simple upper-threshold penalty.
* Construction of neutrosophic truth, indeterminacy, and falsity channels.
* Bootstrapped neutrosophic VIP-based dynamic weights.
* Aesthetic sensitivity handling through `TDS` and, when present, `Fe`; in the current dataset the `Fe` path is inactive.
* Final `NPLS-GWQI` classification into `Excellent`, `Good`, `Poor`, `Very Poor`, and `Unsuitable`.

### 3.4 Compositional hydrochemical discovery (`compositional.py`, `process_discovery.py`)
* ILR transformation of the nine-part hydrochemical composition.
* Parallel analysis for component retention using repeated random eigenvalue comparison.
* PCA in ILR space with robust projection logic when the sample-to-feature ratio allows.
* Back-transformation of components into CLR loadings for geochemical interpretation.
* Conservative rule-based process naming refined later by deterministic geochemical evidence.

### 3.5 Deterministic geochemical constraints (`geochemistry.py`)
* Silicate index and carbonate index.
* `Na/Cl` ratio for salinity-weathering interpretation.
* `CAI_1` and `CAI_2` as ion-exchange indicators.
* Exchange-space coordinates `Exchange_X` and `Exchange_Y`.
* Gibbs-domain assignment using reconstructed TDS and ionic ratios.
* Use of these indices to constrain or refine hydrochemical signature interpretation rather than to claim single-cause mechanisms.

### 3.6 Bayesian endpoint-driver modelling (`bayesian_regression.py`, `generation_script.py`)
* The current locked endpoint-driver analysis is Bayesian and covers `NPLS_GWQI` and `F`.
* Predictors come from the augmented hydrochemical block rather than from the older main-text split-score process contribution design.
* Predictors and targets are standardised before Bayesian horseshoe regression in PyMC.
* Main outputs are posterior mean coefficients, `95%` interval bounds, probability of direction (`pd`), feature selection by directional certainty, fitted values, and driver-score percentages.
* Supporting diagnostics include:
  * a sparse PLS comparator for `F`,
  * model diagnostics in Table 9,
  * supplementary ablation analysis comparing major-ion-only and enhanced augmented predictor sets.

### 3.7 Deterministic and Bayesian health-risk assessment (`risk.py`, `generation_script.py`)
* Deterministic ingestion-only non-carcinogenic risk is computed for `F` and `NO3`.
* The implemented demographic groups are:
  * Adults: `BW = 60 kg`, `IR = 2.5 L/day`
  * Children: `BW = 15 kg`, `IR = 1.0 L/day`
  * Teens: `BW = 45 kg`, `IR = 1.5 L/day`
* Deterministic outputs include `HQ_F`, `HQ_NO3`, `HI_skeletal_dental`, and `HI_hemato`.
* Bayesian risk uses hierarchical uncertainty in body weight and ingestion rate and returns posterior means, `95%` credible intervals, and `P(HI > 1)` for each organ-group combination.
* Arsenic and dermal pathways exist in code but are inactive for the present dataset and contaminants.

### 3.8 Irrigation suitability assessment (`irrigation.py`)
* Calculation of seven irrigation indices in meq/L space:
  * `SAR`, `SSP`, `MH`, `KR`, `PI`, `RSC`, `PS`.
* Entropy-based objective weighting from the index hazard matrix.
* Neutrosophic irrigation scoring and classification into `Unsuitable`, `Marginal`, `Good`, and `Excellent`.
* Main-text visualisation uses `EC` versus `SSP` as a Wilcox proxy and `EC` versus `SAR` as a USSL proxy.

### 3.9 Validation and supplementary diagnostics (`validation.py`, `sensitivity.py`)
* Bootstrap stability of absolute CLR loadings is implemented and reported as a supplementary figure.
* Endpoint model diagnostics include in-sample `R^2` for Bayesian models and nested CV/permutation statistics for the sparse PLS comparator.
* The code can also generate probabilistic sensitivity tables for WQI and ISI, but these outputs are outside the current locked artifact manifest and should not be cited unless promoted deliberately into the manuscript scope.

---

## 4. Results (~1,800-2,500 words)

### 4.1 Study-area geology and sampling distribution
* Geological map interpretation and sampling-site spread (Figure 1).
* Dominance of the `Bimbila Formation`, secondary `Bunya` belt, and southern `Obosum` tract.

### 4.2 Data integrity and hydrochemical characterisation
* Descriptive statistics and guideline exceedance screening (Tables 1-2).
* Charge-balance error distribution and cation-anion agreement (Figure 3).
* ILR hydrochemical facies structure and supporting WHO/GSA-normalised parameter spread (Figure 2; Figure S3).

### 4.3 Groundwater-quality structure and classification
* NPLS-GWQI architecture, neutrosophic component behaviour, and class distribution (Figure 4).
* Variable weighting and dominant water-quality contributors (Table 3).
* Class counts and classwise summary metrics (Table 4).
* Sensory-penalty sensitivity as supporting evidence (Figure S2).

### 4.4 Latent geochemical signatures and deterministic constraints
* Retained process count from parallel analysis (Table 5; Figure S1).
* CLR loading fingerprints and refined signature labels (Table 6; Figure 5).
* Geochemical constraint evidence from Gibbs, exchange, and `Na/Cl` diagnostics (Table 7; Figure 6).
* Bootstrap stability of the process architecture (Figure S6).

### 4.5 Bayesian hydrochemical driver profiles and pathway synthesis
* Main endpoint-driver results for `NPLS-GWQI` and `F` (Tables 8-9; Figure 7).
* Discussion of the driver structure in terms of mineralisation, ion exchange, and weathering proxies.
* Main-text Sankey pathway for `NPLS-GWQI` (Figure 11).
* Supplementary support from the `F` pathway, geochemical control curves, and ablation study (Figures S5, S7, S8).

### 4.6 Deterministic and Bayesian non-carcinogenic health risk
* Deterministic `HQ` and `HI` patterns across Adults, Children, and Teens (Table 10; Figure 8).
* Posterior risk summaries, uncertainty ranges, and exceedance probabilities (Table 11; Figure 9).
* Contrast between fluoride-driven skeletal-dental burden and nitrate-linked hematological burden.

### 4.7 Irrigation suitability
* Irrigation index statistics and entropy weights (Tables 12-13).
* `N-ISI` class distribution and proxy-plot behaviour (Figure 10).
* Interpretation of a predominantly marginal irrigation suitability profile.

---

## 5. Discussion (~1,200-1,800 words)

### 5.1 Hydrochemical evolution in the mapped Voltaian setting
* Relate facies structure and geochemical signatures to the mapped basin-fill lithologies and faults.
* Evaluate rock-weathering dominance, ion exchange, and selective salinity enrichment.

### 5.2 Meaning of the groundwater-quality structure
* Explain why most samples remain in excellent condition while a smaller subset is clearly degraded.
* Link WQI structure to conductivity, dissolved-ion enrichment, and geochemical context.

### 5.3 Interpreting the Bayesian driver profiles
* Discuss why the current locked driver narrative centres `NPLS-GWQI` and `F`, not `NO3`.
* Interpret the balance between promoting and suppressing signatures in the Bayesian models.

### 5.4 Public-health interpretation of fluoride and nitrate
* Explain the stronger fluoride burden and child vulnerability.
* Contrast deterministic and Bayesian risk narratives without overstating certainty.

### 5.5 Irrigation and groundwater-use implications
* Discuss what a mostly marginal `N-ISI` class structure means for irrigation planning.
* Relate index dispersion and class structure to practical water-use caution.

### 5.6 Methodological strengths, limits, and transferability
* Strengths: compositional treatment, neutrosophic weighting, Bayesian uncertainty handling.
* Limits: 34 samples, single campaign, absence of trace metals, inactive code paths for unavailable analytes.
* Transferability to similar Voltaian sedimentary groundwater settings.

---

## 6. Groundwater Management Implications (~350-500 words)

* Borehole screening and targeted intervention where fluoride-related risk is dominant.
* Monitoring and sanitary protection where nitrate or mixed quality stress is suspected.
* Risk communication prioritising children as the most vulnerable demographic in the current implementation.
* Irrigation advisory use of the `N-ISI` classes and proxy plots.

---

## 7. Conclusion (~300-400 words)

* Concise statement of the principal hydrochemical, groundwater-quality, health-risk, and irrigation findings.
* Clear summary of the code-supported methodological contribution.
* Practical limitations and immediate next steps for future work.

---

## Current Locked List of Tables

| Table No. | Title | Current manuscript role |
|---|---|---|
| Table 1 | Descriptive statistics of measured hydrochemical variables | Main text |
| Table 2 | WHO/GSA guideline comparison and exceedance summary | Main text |
| Table 3 | Bootstrapped N-VIP dynamic weights | Main text |
| Table 4 | `NPLS-GWQI` classification summary | Main text |
| Table 5 | Parallel analysis eigenvalues and retention decision | Main text |
| Table 6 | CLR loadings and refined hydrochemical signature labels | Main text |
| Table 7 | Deterministic geochemical constraints | Main text |
| Table 8 | Bayesian endpoint-driver summaries for `NPLS-GWQI` and `F` | Main text |
| Table 9 | Endpoint model diagnostics | Main text |
| Table 10 | Deterministic `HQ` and `HI` by demographic group | Main text |
| Table 11 | Bayesian posterior health-risk summary | Main text |
| Table 12 | Irrigation index statistics | Main text |
| Table 13 | Entropy weights and `N-ISI` class distribution | Main text |
| Table S1 | Full sample-level WQI output | Supplementary |
| Table S2 | Sample-level process scores and endpoint contribution matrices | Supplementary |

### Optional code-generated tables outside the current locked manifest
* `Table14_WQI_Probabilities.csv`
* `Table15_ISI_Probabilities.csv`

These sensitivity outputs are produced by code when available but are not part of the current locked Stage 7/9 artifact sequence.

---

## Current Locked List of Figures

| Figure No. | Title | Current panel logic |
|---|---|---|
| Figure 1 | Geological and location map of sampling sites in Karaga District | External GIS/cartographic map |
| Figure 2 | ILR hydrochemical characterisation | 2 × 2 ILR facies panels |
| Figure 3 | Charge-balance error and cation-anion agreement | Histogram + 1:1 scatter |
| Figure 4 | `NPLS-GWQI` structure | Components vs score, VIP radar, class boxplot, class distribution |
| Figure 5 | Latent hydrochemical signature fingerprint | CLR loading bar chart |
| Figure 6 | Deterministic geochemical constraints | Gibbs cations, Gibbs anions, exchange plot, `Na` vs `Cl` |
| Figure 7 | Bayesian endpoint-driver profiles | `NPLS-GWQI` and `F` panels |
| Figure 8 | Deterministic non-carcinogenic risk | `HQ_F`, `HQ_NO3`, organ-specific mean `HI` |
| Figure 9 | Bayesian posterior health-risk summaries | Skeletal-dental violin, hematological violin, forest summary |
| Figure 10 | Irrigation suitability | Wilcox proxy, USSL proxy, index boxplots, class distribution |
| Figure 11 | `NPLS-GWQI` context-signature-impact pathway | Main-text Sankey pathway |
| Figure S1 | Parallel-analysis scree plot | Supplementary |
| Figure S2 | Sensory-penalty sensitivity of `NPLS-GWQI` | Supplementary |
| Figure S3 | WHO/GSA-normalised parameter box-and-whisker plot | Supplementary |
| Figure S5 | Exploratory fluoride pathway | Supplementary |
| Figure S6 | Bootstrap stability of absolute CLR loadings | Supplementary |
| Figure S7 | Endpoint response to ion-exchange and weathering proxies | Supplementary |
| Figure S8 | Model ablation study | Supplementary |

---

## Notes for the `r2m` Workflow

* Stage 7 execution should follow the current `generation_script.py` and the validated `artifact_manifest.yaml`, not the older figure/table ordering.
* Figure 1 is external and should be treated as a locked map asset rather than a generated analytical figure.
* The current manuscript narrative should cite Figures 1-11 and Tables 1-13 as the main text sequence, with the supplementary set listed above.
* The package still contains legacy or optional paths, but manuscript drafting should follow the current locked artifact set unless the manifest, outline, and proposal are revised together.
