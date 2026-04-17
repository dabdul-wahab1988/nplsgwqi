# Tables

Note: `manuscript/frameworks/locked-analysis-plan.yaml` is not present in this workspace. The `Locked plan source` fields below therefore cite the matching IDs from `manuscript/artifact_manifest.yaml` provisionally. Only manifest-listed table artifacts are included here so the Stage 9 numbering remains aligned with the validated artifact set.

## Main Text Tables

### Table 1. Descriptive Statistics
**Artifact:** `manuscript/artifacts/Table1_Descriptive_Statistics.csv`
**Objective(s):** `OBJ-1`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-1`
**Caption:** Descriptive statistics for field parameters and major ions in the 34 sampled groundwater sources.
**Technical Description:** This table reports summary statistics for the measured groundwater chemistry, including central tendency and spread for pH, salinity indicators, and dissolved major ions. The values are computed directly from the analytical dataset for the full sample set.
**Scientific Description:** The groundwater is mildly alkaline on average but chemically heterogeneous, with strong upper-tail enrichment in `EC`, `TDS`, `Na`, `K`, `Cl`, and `F`. The distribution therefore reflects a mostly moderate background chemistry punctuated by a smaller number of highly mineralised waters.

### Table 2. WHO/GSA Compliance
**Artifact:** `manuscript/artifacts/Table2_WHO_GSA_Compliance.csv`
**Objective(s):** `OBJ-1`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-2`
**Caption:** Compliance of groundwater parameters with WHO and GSA guideline values, reported as exceedance counts and percentages.
**Technical Description:** The table lists each parameter together with its guideline threshold, exceedance count, exceedance percentage, and basic summary values for the sampled dataset. It is a parameter-by-parameter compliance audit against the adopted drinking-water standards.
**Scientific Description:** `HCO3` shows the most frequent exceedance, while `Na`, `K`, and `F` each exceed in `17.65%` of samples. `EC`, `TDS`, and `Cl` exceed less often, and `Mg` and `SO4` remain fully compliant in this dataset, indicating that non-compliance is concentrated in specific dissolved constituents rather than distributed evenly across all parameters.

### Table 3. VIP Weights
**Artifact:** `manuscript/artifacts/Table3_VIP_Weights.csv`
**Objective(s):** `OBJ-2`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-3`
**Caption:** Variable-importance weights for the hydrochemical inputs contributing to the `NPLS-GWQI` model.
**Technical Description:** This table ranks the hydrochemical predictors by their variable-importance contribution within the groundwater quality index workflow. The weights provide a normalised measure of relative influence across the candidate parameters.
**Scientific Description:** `EC`, `Cl`, `TDS`, `Ca`, `K`, and `Na` have the highest weights, showing that mineralisation and salinity-related variables dominate the water-quality signal encoded by the index.

### Table 4. GWQI Classification
**Artifact:** `manuscript/artifacts/Table4_GWQI_Classification.csv`
**Objective(s):** `OBJ-2`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-4`
**Caption:** Classwise summary of `NPLS-GWQI` scores and neutrosophic membership components for the sampled groundwater.
**Technical Description:** The table reports the number of samples in each groundwater-quality class together with mean score and classwise summaries of truth, indeterminacy, and falsity. The values summarise the final neutrosophic classification output for all 34 samples.
**Scientific Description:** Excellent water dominates the sample set (`25` of `34` samples), whereas only a small subset falls into poor, very poor, or unsuitable classes. Indeterminacy is comparatively low at the most extreme classes, indicating that the worst waters are separated clearly rather than weakly from the main body of acceptable samples.

### Table 5. Parallel Analysis
**Artifact:** `manuscript/artifacts/Table5_Parallel_Analysis.csv`
**Objective(s):** `OBJ-3`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-5`
**Caption:** Parallel-analysis results used to determine the number of retained compositional components.
**Technical Description:** The table compares observed component eigenvalues with the random reference distribution and the `95th` percentile random threshold. Components are retained when the observed eigenvalue exceeds the random benchmark.
**Scientific Description:** Only `PC1` to `PC5` satisfy the retention rule, supporting a five-component representation of the compositional hydrochemical structure.

### Table 6. CLR Loadings
**Artifact:** `manuscript/artifacts/Table6_CLR_Loadings.csv`
**Objective(s):** `OBJ-3`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-6`
**Caption:** Centred log-ratio loading structure and refined process labels for the five retained hydrochemical components.
**Technical Description:** This table lists the ion loadings associated with each retained component and the refined process name assigned after interpretation of the loading structure. The signs and magnitudes indicate how each ion contributes to each latent process.
**Scientific Description:** The component solution resolves separable signatures for `Na-Cl` enrichment, `NO3-Ca` association, `Mg-F` association, `K-F` association, and `Na-HCO3` exchange, indicating that several distinct hydrochemical regimes contribute to the observed groundwater chemistry.

### Table 7. Geochemical Constraints
**Artifact:** `manuscript/artifacts/Table7_Geochemical_Constraints.csv`
**Objective(s):** `OBJ-3`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-7`
**Caption:** Summary of deterministic geochemical diagnostic metrics, including Gibbs-field counts and exchange-weathering proxies.
**Technical Description:** The table reports the distribution of samples across Gibbs diagnostic fields together with summary statistics for `Na_Cl_Ratio`, `CAI_1`, and the carbonate-weathering index. These metrics are derived from the transformed hydrochemical dataset and used for process auditing.
**Scientific Description:** Rock-weathering is the dominant field (`30` samples), the mean `CAI_1` is strongly negative, and `Na/Cl` ratios are generally elevated. Together these metrics support rock-weathering control with pervasive ion exchange and sodium enrichment beyond a simple `Na = Cl` dissolution relationship.

### Table 8. Process Influence
**Artifact:** `manuscript/artifacts/Table8_Process_Influence.csv`
**Objective(s):** `OBJ-4`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-8`
**Caption:** Bayesian driver summaries for `NPLS-GWQI` and fluoride, reporting posterior direction, directional certainty, and normalised contribution share.
**Technical Description:** For each endpoint and predictor, the table reports the posterior mean coefficient (`beta_mean`), probability of direction (`pd`), driver score, and contribution percentage derived from the absolute coefficient magnitude. The predictors come from the augmented hydrochemical feature set.
**Scientific Description:** `EC` is the dominant driver of `NPLS-GWQI` (`43.18%` contribution), whereas the fluoride endpoint is split between `CAI_2` and `EC` at comparable magnitude but opposite sign. The ranked pattern indicates that conductivity, ion exchange, and proxy variables jointly structure the endpoint responses.

### Table 9. Model Diagnostics
**Artifact:** `manuscript/artifacts/Table9_Model_Diagnostics.csv`
**Objective(s):** `OBJ-4`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-9`
**Caption:** Comparative diagnostics for the endpoint models, including feature sparsity and in-sample or cross-validated fit summaries.
**Technical Description:** The table summarises endpoint-model diagnostics such as retained feature count, `R^2`, and related fit statistics across the Bayesian horseshoe and sparse PLS formulations. It is a model-performance audit rather than a hydrochemical result table.
**Scientific Description:** The `NPLS-GWQI` Bayesian model collapses onto a highly concentrated signal dominated by `EC` and achieves very high in-sample fit, whereas the fluoride models are weaker and less stable, especially under sparse PLS cross-validation. The fluoride endpoint should therefore be interpreted more cautiously than the groundwater-quality endpoint.

### Table 10. Deterministic Risk
**Artifact:** `manuscript/artifacts/Table10_Deterministic_Risk.csv`
**Objective(s):** `OBJ-5`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-10`
**Caption:** Deterministic hazard quotients and health-index summaries for adults, children, and teens exposed to groundwater fluoride and nitrate.
**Technical Description:** The table reports group-specific risk summaries for `HQ_F`, `HQ_NO3`, and aggregated endpoint burdens under the deterministic exposure assumptions used in the health-risk workflow. It includes summary measures such as means and upper-tail values for each demographic group.
**Scientific Description:** Children have the highest deterministic fluoride burden and are the only group with mean `HQ_F` above `1`, while nitrate hazard remains lower across all three groups. The main deterministic concern is therefore fluoride-related skeletal-dental burden in younger receptors.

### Table 11. Bayesian Risk
**Artifact:** `manuscript/artifacts/Table11_Bayesian_Risk.csv`
**Objective(s):** `OBJ-5`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-11`
**Caption:** Bayesian posterior summaries of non-carcinogenic health indices for adults, children, and teens.
**Technical Description:** This table reports posterior mean health indices, `95%` credible intervals, and exceedance probabilities `P(HI > 1)` for each demographic group and health-effect category. It propagates uncertainty through the population-level risk calculation instead of using a single deterministic estimate.
**Scientific Description:** Children again show the highest skeletal-dental burden, with the largest posterior mean and the greatest probability of exceeding `HI = 1`, whereas hematological risk stays comparatively low across all groups. The uncertainty-aware results therefore preserve the same vulnerability ordering as the deterministic analysis.

### Table 12. Irrigation Indices
**Artifact:** `manuscript/artifacts/Table12_Irrigation_Indices.csv`
**Objective(s):** `OBJ-6`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-12`
**Caption:** Descriptive summaries of the irrigation indices used to construct the neutrosophic irrigation suitability index (`N-ISI`).
**Technical Description:** The table reports summary statistics for the irrigation metrics `SAR`, `SSP`, `MH`, `KR`, `PI`, `RSC`, `PS`, and the derived `N_ISI_Score`. These are the direct numerical inputs to the irrigation suitability classification workflow.
**Scientific Description:** The irrigation indicators are strongly dispersed, with particularly large spread in `KR` and `PI`, showing that a subset of samples carries disproportionately poor irrigation characteristics against a broader background of moderate conditions.

### Table 13. NISI Classification
**Artifact:** `manuscript/artifacts/Table13_NISI_Classification.csv`
**Objective(s):** `OBJ-6`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-13`
**Caption:** Entropy weights of irrigation indices and final class distribution for the neutrosophic irrigation suitability index (`N-ISI`).
**Technical Description:** The table combines the entropy-derived weights for each irrigation index with the final count of samples in each `N-ISI` class. It therefore links variable weighting and class outcome within the same irrigation suitability framework.
**Scientific Description:** `PI`, `PS`, and `KR` contribute the greatest weight to the composite index, and the marginal class dominates the final classification (`23` samples). Irrigation suitability is therefore constrained primarily by a small set of high-leverage salinity-alkalinity indicators rather than by uniform contribution across all indices.

## Supplementary Tables

## Revision Tables (Reviewer-response additions)

### Table 16. Charge-balance error (per sample)
**Artifact:** `manuscript/artifacts/Table16_CBE_PerSample.csv`
**Objective(s):** `OBJ-1`
**Caption:** Per-sample charge-balance error (CBE, %) computed from meq/L major ions.

### Table 17. Charge-balance error (summary)
**Artifact:** `manuscript/artifacts/Table17_CBE_Summary.csv`
**Objective(s):** `OBJ-1`
**Caption:** Summary statistics for the CBE distribution, including pass fractions under ±5% and ±10%.

### Table 18. Saturation indices (PHREEQC)
**Artifact:** `manuscript/artifacts/Table18_Saturation_Indices.csv`
**Objective(s):** `OBJ-3`
**Caption:** PHREEQC saturation indices for calcite, dolomite, and fluorite under a 20/25/30 °C temperature band.

### Table 19. Bayesian driver-model diagnostics
**Artifact:** `manuscript/artifacts/Table19_Bayes_Diagnostics_DriverModels.csv`
**Objective(s):** `OBJ-4`
**Caption:** Sampler diagnostics (R-hat/ESS) and predictive diagnostics (PSIS-LOO/WAIC) for endpoint-driver models.

### Table 20. Horseshoe prior sensitivity
**Artifact:** `manuscript/artifacts/Table20_Prior_Sensitivity_Horseshoe.csv`
**Objective(s):** `OBJ-4`
**Caption:** Stability of ranked driver attributions under alternative global shrinkage scales (τ prior).

### Table 21. ILR basis sensitivity
**Artifact:** `manuscript/artifacts/Table21_ILR_Basis_Sensitivity.csv`
**Objective(s):** `OBJ-3`
**Caption:** Similarity of absolute CLR loading patterns under alternative orthonormal ILR bases.

### Table 22. WQI benchmarking
**Artifact:** `manuscript/artifacts/Table22_WQI_Benchmarking.csv`
**Objective(s):** `OBJ-2`
**Caption:** Sample-level comparison of NPLS-GWQI versus conventional arithmetic WQI and a fuzzy-only baseline.

### Table 23. Conventional WQI weight uncertainty
**Artifact:** `manuscript/artifacts/Table23_ConventionalWQI_WeightUncertainty.csv`
**Objective(s):** `OBJ-2`
**Caption:** Class probabilities for the conventional WQI under Monte Carlo weight perturbation.

### Table 24. Irrigation indices (per sample)
**Artifact:** `manuscript/artifacts/Table24_Irrigation_PerSample.csv`
**Objective(s):** `OBJ-6`
**Caption:** Per-sample irrigation indices with neutrosophic and classical categories for comparison.

### Table 25. Irrigation disagreement summary
**Artifact:** `manuscript/artifacts/Table25_Irrigation_Disagreement.csv`
**Objective(s):** `OBJ-6`
**Caption:** Disagreement rates between N-ISI and SAR/SSP/RSC-based classical screens.

### Table 26. Risk uncertainty decomposition
**Artifact:** `manuscript/artifacts/Table26_Risk_Uncertainty_Decomposition.csv`
**Objective(s):** `OBJ-5`
**Caption:** Variance-share summaries separating concentration-driven from exposure-rate-driven uncertainty.

### Table 27. Risk exceedance probabilities (multiple thresholds)
**Artifact:** `manuscript/artifacts/Table27_Risk_Threshold_Probabilities.csv`
**Objective(s):** `OBJ-5`
**Caption:** Exceedance probabilities for HI thresholds 0.5, 1.0, and 2.0 by demographic group.

### Table S1. Full WQI
**Artifact:** `manuscript/artifacts/TableS1_Full_WQI.csv`
**Objective(s):** `OBJ-2`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-S1`
**Caption:** Sample-level `NPLS-GWQI` scores, classes, and neutrosophic membership components for all analysed groundwater samples.
**Technical Description:** This supplementary table lists the full per-sample output of the groundwater quality workflow, including sample identifier, `NPLS-GWQI` score, final class, and the associated truth, indeterminacy, and falsity values.
**Scientific Description:** The sample-level results confirm that degraded groundwater quality is concentrated in a limited subset of sites, while most samples remain within the excellent class and preserve a low-score profile.

### Table S2. Sample Contributions
**Artifact:** `manuscript/artifacts/TableS2_Sample_Contributions.csv`
**Objective(s):** `OBJ-3`, `OBJ-4`
**Locked plan source:** Provisional replacement for missing locked plan; `manuscript/artifact_manifest.yaml` -> `TAB-S2`
**Caption:** Sample-specific latent process scores and endpoint contribution terms across the retained geochemical and Bayesian driver framework.
**Technical Description:** The table reports per-sample scores for `Process_1` to `Process_5` together with the sample-level contribution terms generated for endpoint modelling. It is the detailed trace table connecting latent process structure to individual sample influence patterns.
**Scientific Description:** The supplementary detail shows that both latent-process expression and endpoint contribution are sample dependent, allowing the most influential waters and the most active process regimes to be traced directly rather than inferred only from aggregated summaries.
