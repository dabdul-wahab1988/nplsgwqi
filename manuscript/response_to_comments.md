# Response to Reviewer Comments

This document responds comment-by-comment to `manuscript/reviewer_comments.md` (verbatim anchor).

Statuses: `Addressed`, `Partially addressed`, `Not possible with available data`.

## C01

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C01` | This manuscript tackles an important and impactful problem; groundwater quality and fluoride-related health risk in semi-arid sedimentary systems, using a thoughtfully integrated framework that addresses compositional closure, weight uncertainty in indexing, and probabilistic risk. The conceptual contribution is meaningful: the structural coupling of compositional process discovery with uncertainty-aware assessment is a step forward from stacked, uncoupled pipelines frequently seen in the literature. The case study is policy-relevant, and the probabilistic framing for vulnerable groups is valuable. | Thank you. No change required beyond incorporating the substantive revisions listed below. | N/A (commendation; no manuscript change) | Addressed |

## C02

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C02` | However, several elements require significant strengthening before this work meets the standard for publication in Scientific Reports. Chief among these are (i) limited sample size and single-season design without robust cross-validation or sensitivity analyses in the Bayesian layer; (ii) insufficient benchmarking of the proposed indices and models against established alternatives to demonstrate practical benefit; (iii) missing analytical QA/QC details crucial for data credibility; and (iv) incomplete mechanistic support for fluoride mobilization in the absence of speciation/saturation modeling. The writing is generally clear, but some essential implementation details are relegated to the Supplementary or remain unspecified in the main text. | Implemented reviewer-requested strengthening where supported: explicit CBE reporting, PHREEQC SI analysis, WQI benchmarking under weight uncertainty, expanded Bayesian diagnostics reporting, and added reproducibility metadata (workflow + data dictionary). Limitations are stated where required metadata are unavailable. | `manuscript/sections/03-methodology/subsection-02-preprocessing.md#L7` (CBE thresholding and use.); `manuscript/sections/03-methodology/subsection-05-geochemistry.md#L9` (Speciation and saturation-index support (PHREEQC).); `manuscript/sections/03-methodology/subsection-03-npls-gwqi.md#L13` (Benchmarking against conventional alternatives.); `manuscript/sections/03-methodology/subsection-06-bayesian.md#L31` (Sampler configuration and diagnostics.); `manuscript/sections/03-methodology/subsection-09-validation.md#L9` (analysis workflow diagram and a concise data dictionary); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C03

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C03` | Detailed comments are given below: | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C04

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C04` | Technical soundness evaluation | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C05

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C05` | The overall architecture is coherent and addresses real methodological shortcomings (closure, fixed weights, deterministic risk). However, the model evaluation layer is underdeveloped: lack of PSIS-LOO/WAIC, posterior predictive checks, sensitivity to priors, and external validation undermines confidence in the reported driver attributions and exceedance probabilities. | Added PSIS-LOO/WAIC reporting and PPC plots for the endpoint-driver layer and referenced them in Methods/Supplementary (Table 19; Figure S10). | `manuscript/sections/03-methodology/subsection-06-bayesian.md#L31` (Sampler configuration and diagnostics.); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C06

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C06` | The horseshoe prior is appropriate, but with p not negligible relative to n (even after exclusions), model identifiability and stability should be scrutinized via cross-validation, shrinkage diagnostics (effective number of non-zero coefficients), and sensitivity to the global scale prior. | Added PSIS-LOO/WAIC and prior-scale sensitivity (Table 19–Table 20). External validation is not possible with the available data (single n=34 snapshot, no holdout set), so results are qualified accordingly. | `manuscript/sections/03-methodology/subsection-06-bayesian.md#L31` (Sampler configuration and diagnostics.); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Partially addressed |

## C07

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C07` | ILR-PCA with parallel analysis is suitable; back-interpretation in CLR space must be executed carefully (loadings in ILR do not trivially map to raw-ion effects). Please elaborate the back-transformation and interpretation pipeline to ensure process labels are not artefacts of the chosen SBP. | Expanded the ILR→CLR back-interpretation narrative and added basis/SBP sensitivity analysis (Table 21) to demonstrate stability of dominant loading patterns. | `manuscript/sections/03-methodology/subsection-04-compositional.md#L13` (SBP/basis sensitivity.); `manuscript/sections/04-results/subsection-04-latent-processes.md#L9` (Saturation-index evaluation); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C08

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C08` | Experimental evaluation assessment | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C09

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C09` | Field procedures: Pumping for ~5 minutes may be insufficient unless tied to well volume; please justify purging criteria. Report trip blanks, field duplicates, and any inter-lab comparison if available. | Clarified the ~5-minute purge protocol and in-situ field-parameter measurement. Well volumes/stabilization logs and trip blank/duplicate metadata were not available in the provided materials and are therefore stated as limitations. | `manuscript/sections/03-methodology/subsection-01-dataset.md#L9` (Representativeness and metadata limitations.) | Partially addressed |

## C10

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C10` | Analytical QA/QC: Provide explicit CBE acceptance threshold (±5% or ±10%), percentage of samples within range, and distribution of CBE. Report LOD/LOQ for all analytes and handling of any censored values. | Specified CBE thresholds (±10% primary; ±5% reported) and reported pass fractions and distribution (Table 16–Table 17; Figure 3). LOD/LOQ metadata and censoring flags were not available and are stated as a limitation. | `manuscript/sections/03-methodology/subsection-02-preprocessing.md#L7` (CBE thresholding and use.); `manuscript/sections/04-results/subsection-02-data-integrity.md#L3` (Data integrity was screened); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Partially addressed |

## C11

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C11` | Seasonal coverage: With only dry-season data, emphasize that exceedance probabilities are conditional on this season. Consider scenario analysis incorporating plausible wet-season dilution or draw on nearby longitudinal datasets to inform prior distributions. | Reinforced that exceedance probabilities are conditional on dry-season sampling and added an explicit seasonality limitation statement. Wet-season scenario analysis is not defensible without additional seasonal data or external longitudinal datasets. | `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Partially addressed |

## C12

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C12` | Comparison with related work | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C13

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C13` | The paper’s central novelty is integration rather than invention of new algorithms. Prior studies have applied compositional methods (ILR/CLR), fuzzy/neutrosophic indices, and Bayesian risk, but typically in isolation. This manuscript advances the field by structurally coupling these components; however, it would benefit from head-to-head comparisons with (i) deterministic WQIs under Monte Carlo weight perturbation, (ii) fuzzy vs neutrosophic indices, and (iii) standard regression vs horseshoe for driver profiling to quantify practical gains. | Added head-to-head benchmarking vs conventional arithmetic WQI and a fuzzy-only baseline, plus Monte Carlo weight-uncertainty analysis (Table 22–Table 23; Figure S11). | `manuscript/sections/03-methodology/subsection-03-npls-gwqi.md#L13` (Benchmarking against conventional alternatives.); `manuscript/sections/04-results/subsection-03-gwqi.md#L7` (Benchmarking against conventional indexing alternatives); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C14

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C14` | Expand on how this framework materially changes decisions compared to a conventional pipeline in the same study area (e.g., re-ranking of sites at the potability/irrigation interface, different remediation priorities). | Expanded decision relevance by reporting where class assignments diverge under conventional vs neutrosophic indexing and under weight uncertainty, and by stating operational implications for boundary/degraded sites. | `manuscript/sections/04-results/subsection-03-gwqi.md#L7` (Benchmarking against conventional indexing alternatives); `manuscript/sections/06-management-implications/section.md#L7` (Benchmarking against conventional WQI baselines) | Addressed |

## C15

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C15` | Discussion of broader impact and significance | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C16

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C16` | The demographic disaggregation in probabilistic risk assessment is valuable for policy targeting. The approach is potentially transferable across Voltaian and other semi-arid sedimentary contexts where fluoride risk is endemic. | Retained demographic disaggregation and strengthened multi-threshold exceedance reporting (Table 27) to support policy targeting. | `manuscript/sections/04-results/subsection-06-health-risk.md#L7` (To support multi-threshold decision use) | Addressed |

## C17

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C17` | The framework could serve as a template for regulators if accompanied by transparent, open-source code and minimal data requirements; otherwise, method complexity may hinder uptake. Clarify data and code availability, including the neutrosophic_pls package and Bayesian modeling scripts. | Added explicit code availability language and provided workflow + data dictionary materials to support regulator-facing reproducibility. | `manuscript/sections/07-conclusion/section.md#L37` (**Code Availability:**); `manuscript/analysis_workflow.md#L1` (# Analysis Workflow); `manuscript/data_dictionary.md#L1` (# Data Dictionary) | Addressed |

## C18

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C18` | Reproducibility and transparency | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C19

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C19` | Many critical steps are in the Supplementary Material or external repositories; for a methods-centric paper, more implementation specifics should be in the main text (e.g., neutrosophic membership functions, aesthetic penalty definition, bootstrapping scheme, PCA robustness settings, Bayesian sampler settings, convergence diagnostics: R-hat, ESS). | Moved/added key implementation specifics into the Methods and Supplementary addendum (diagnostics, benchmarking, basis sensitivity, PHREEQC SI) rather than leaving them implicit. | `manuscript/sections/03-methodology/subsection-06-bayesian.md#L31` (Sampler configuration and diagnostics.); `manuscript/sections/03-methodology/subsection-09-validation.md#L9` (analysis workflow diagram and a concise data dictionary); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C20

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C20` | Provide an analysis workflow diagram and a data dictionary. Share input data (anonymized) with metadata to enable replication. | Added an analysis workflow diagram and a concise data dictionary to support replication from the included `data.csv`. | `manuscript/analysis_workflow.md#L1` (# Analysis Workflow); `manuscript/data_dictionary.md#L1` (# Data Dictionary); `manuscript/sections/03-methodology/subsection-09-validation.md#L9` (analysis workflow diagram and a concise data dictionary) | Addressed |

## C21

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C21` | Data quality and sampling strategy | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C22

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C22` | Clarify the representativeness of WVI-maintained public boreholes versus private wells/shallow hand-dug wells. Describe spatial coverage (per formation), well depths, screen intervals, and proximity to potential nitrate sources. | Clarified the sampling frame as WVI-maintained public boreholes and stated that well depth/screen interval and proximity-to-source metadata are unavailable in the provided materials, limiting representativeness inference. | `manuscript/sections/03-methodology/subsection-01-dataset.md#L9` (Representativeness and metadata limitations.) | Partially addressed |

## C23

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C23` | Report temperature and redox proxies (e.g., DO, Eh) if available; these are relevant to nitrogen species and carbonate equilibria. | DO/Eh and temperature fields are not present in `data.csv`. For saturation indices, we report a temperature sensitivity band (20/25/30 °C) and disclose this assumption; redox-proxy reporting is not possible without additional measurements. | `manuscript/sections/03-methodology/subsection-05-geochemistry.md#L9` (Speciation and saturation-index support (PHREEQC).) | Not possible with available data |

## C24

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C24` | Statistical methodology | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C25

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C25` | For Bayesian models, include: number of chains, warmup and sampling iterations, adaptation parameters, R-hat and ESS, posterior predictive checks (PPC plots), LOO model diagnostics, and sensitivity to prior scales (e.g., alternative τ priors, regularized horseshoe). | Reported Bayesian sampler diagnostics, PPC, PSIS-LOO/WAIC, and prior sensitivity summaries (Table 19–Table 20; Figure S10; Supplementary S5.2–S5.3). | `manuscript/sections/03-methodology/subsection-06-bayesian.md#L31` (Sampler configuration and diagnostics.); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C26

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C26` | For NPLS-GWQI, compare dynamic weights to expert-assigned weights and quantify variability across bootstraps (e.g., confidence intervals for N-VIP). Report how indeterminacy I and falsity F affect classification stability. | Dynamic weights and classification stability are already reported; we further added benchmarking vs a conventional equal-weight WQI and a fuzzy-only baseline under weight uncertainty (Table 22–Table 23). A direct comparison to expert-assigned weights is not possible because a defensible expert-weight vector was not provided. | `manuscript/sections/04-results/subsection-03-gwqi.md#L7` (Benchmarking against conventional indexing alternatives); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Partially addressed |

## C27

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C27` | Hydrogeochemical interpretation and geochemical modeling | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C28

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C28` | Interpretations of fluoride mobilization would be strengthened by speciation and saturation indices (calcite, dolomite, fluorite, fluorapatite), and ion activity corrections. If fluorite is undersaturated and calcite is near saturation, ion exchange and carbonate equilibria can be better disentangled. | Added PHREEQC saturation indices for calcite/dolomite/fluorite with temperature sensitivity (Table 18; Figure S9). Fluorapatite SI cannot be computed without phosphate data and is disclosed as a limitation. | `manuscript/sections/03-methodology/subsection-05-geochemistry.md#L9` (Speciation and saturation-index support (PHREEQC).); `manuscript/sections/04-results/subsection-04-latent-processes.md#L9` (Saturation-index evaluation); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Partially addressed |

## C29

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C29` | Consider plotting Na-normalized or charge-balanced diagrams, and address whether halite dissolution is negligible given Na/Cl trends. | Halite-dissolution plausibility is addressed via Na/Cl diagnostics; Na-normalized diagram additions were not implemented in this revision and are noted as a potential future enhancement. | `manuscript/sections/04-results/subsection-04-latent-processes.md#L9` (Saturation-index evaluation) | Partially addressed |

## C30

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C30` | If “conductivity enrichment” is shorthand for evaporation/concentration, align terminology with Gibbs fields and mass-balance reasoning. | Aligned terminology by interpreting conductivity enrichment explicitly as evapoconcentration/mineralisation where relevant. | `manuscript/sections/04-results/subsection-04-latent-processes.md#L9` (Saturation-index evaluation) | Addressed |

## C31

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C31` | Risk assessment modeling | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C32

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C32` | Specify exposure pathways modeled, parameter distributions used (ingestion rates, body weights, exposure frequency/duration), and their sources (local vs literature). Provide tornado plots and contribution-to-variance to highlight dominant uncertainty drivers. | Clarified ingestion-only pathway choice and pointed to Supplementary distributions; added an uncertainty-driver decomposition (Table 26) as a compact analogue to tornado-style uncertainty ranking under the current model structure. | `manuscript/sections/03-methodology/subsection-07-health-risk.md#L9` (Pathways, thresholds, and uncertainty drivers.); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Partially addressed |

## C33

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C33` | Clarify the “skeletal-dental hazard threshold” definition (HQ>1 vs specific benchmarks), and report probabilities for multiple cutoffs. If nitrate dermal exposure is negligible, justify excluding it or show it is immaterial. | Clarified that HI>1 is the primary hazard threshold and added multi-threshold exceedance probabilities (Table 27). Justified ingestion-only modelling and explicitly reported which pathways are excluded. | `manuscript/sections/03-methodology/subsection-07-health-risk.md#L9` (Pathways, thresholds, and uncertainty drivers.); `manuscript/sections/04-results/subsection-06-health-risk.md#L7` (To support multi-threshold decision use) | Addressed |

## C34

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C34` | Irrigation suitability | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C35

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C35` | List the seven parameters used, their standards or reference ranges, and the entropy weighting mechanics. Compare neutrosophic ISI classifications to classical SAR/Na%/RSC-based categorizations, including disagreement analysis. | Ensured the seven irrigation indices are explicitly listed; added per-sample reporting and disagreement summary vs classical SAR/SSP/RSC criteria (Table 24–Table 25). | `manuscript/sections/03-methodology/subsection-08-irrigation.md#L15` (compared against classical); `manuscript/sections/04-results/subsection-07-irrigation.md#L5` (Per-sample irrigation indices); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C36

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C36` | Questions for Clarification | Noted (comment section header). | N/A (section header in comments) | Addressed |

## C37

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C37` | What explicit charge-balance error (CBE) threshold did you use, and what fraction of samples passed it? Please provide the CBE distribution. | Provided explicit CBE thresholds, pass fractions, and distribution summary (Table 16–Table 17; Figure 3). | `manuscript/sections/03-methodology/subsection-02-preprocessing.md#L7` (CBE thresholding and use.); `manuscript/sections/04-results/subsection-02-data-integrity.md#L3` (Data integrity was screened) | Addressed |

## C38

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C38` | How were detection limits (LOD/LOQ) for each analyte handled? Were any measurements censored, and if so, how were they treated in ILR and Bayesian analyses? | LOD/LOQ and censoring metadata are not present in the supplied dataset materials, so censored-value handling cannot be specified; this is reported explicitly as a limitation. | `manuscript/sections/03-methodology/subsection-02-preprocessing.md#L7` (CBE thresholding and use.) | Not possible with available data |

## C39

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C39` | How did you determine that a ~5-minute purge was sufficient for each borehole? Can you report well volumes or stabilization criteria (pH/EC/temperature) used to confirm representativeness? | Clarified the stated purge practice (~5 minutes) and noted that per-well volumes and stabilization logs are not available, limiting further justification. | `manuscript/sections/03-methodology/subsection-01-dataset.md#L9` (Representativeness and metadata limitations.) | Partially addressed |

## C40

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C40` | Please detail the ILR to CLR back-interpretation of loadings and how SBP choices affect process labeling. Did you test alternative SBPs to confirm stability of conclusions? | Added an explicit ILR→CLR back-interpretation explanation and basis sensitivity summary (Table 21). | `manuscript/sections/03-methodology/subsection-04-compositional.md#L13` (SBP/basis sensitivity.); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C41

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C41` | For the horseshoe models, please provide sampler diagnostics (R-hat, ESS), posterior predictive checks, PSIS-LOO or WAIC, and prior sensitivity analyses. How many predictors retained non-negligible mass, and were results robust to alternative shrinkage priors? | Added sampler diagnostics, PPC, PSIS-LOO/WAIC, and prior sensitivity outputs (Table 19–Table 20; Figure S10). Diagnostics are reported transparently, and weaker stability in the fluoride model is explicitly qualified in interpretation. | `manuscript/sections/03-methodology/subsection-06-bayesian.md#L31` (Sampler configuration and diagnostics.); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Partially addressed |

## C42

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C42` | Can you benchmark NPLS-GWQI against at least one conventional WQI and a fuzzy index under Monte Carlo weight uncertainty to quantify classification stability and practical decision differences? | Implemented benchmarking against conventional and fuzzy-only WQI baselines and quantified conventional-WQI class stability under Monte Carlo weight uncertainty (Table 22–Table 23; Figure S11). | `manuscript/sections/03-methodology/subsection-03-npls-gwqi.md#L13` (Benchmarking against conventional alternatives.); `manuscript/sections/04-results/subsection-03-gwqi.md#L7` (Benchmarking against conventional indexing alternatives); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C43

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C43` | Which exposure pathways and parameter distributions were used for the hierarchical Bayesian risk model, and are they locally derived? Can you provide exceedance probabilities for multiple HQ thresholds and a contribution-to-variance breakdown? | Clarified exposure pathways (ingestion-only), provided multi-threshold exceedance probabilities (Table 27), and added an uncertainty-driver decomposition (Table 26). Parameter distributions are documented in the Supplementary methodology. | `manuscript/sections/03-methodology/subsection-07-health-risk.md#L9` (Pathways, thresholds, and uncertainty drivers.); `manuscript/sections/04-results/subsection-06-health-risk.md#L7` (To support multi-threshold decision use); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Addressed |

## C44

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C44` | Did you perform speciation/saturation modeling (e.g., PHREEQC) to test fluoride mobilization hypotheses (fluorite saturation, calcite/dolomite buffering)? If not, can you add this analysis or justify its omission? | Added PHREEQC saturation-index analysis for key carbonate/fluoride phases (Table 18; Figure S9) and disclosed limitations (temperature not in `data.csv`; phosphate not measured for fluorapatite SI). | `manuscript/sections/03-methodology/subsection-05-geochemistry.md#L9` (Speciation and saturation-index support (PHREEQC).); `manuscript/supplementary-methodology.md#L936` (# S5. Reviewer-requested additions) | Partially addressed |

## C45

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C45` | Please clarify the regional affiliation of Karaga District (Northern Region vs North East Region) and correct any inconsistencies. | Standardised the manuscript to state that Karaga District is in the Northern Region of Ghana and removed the inconsistent North East Region phrasing. | `manuscript/sections/02-study-area/subsection-01-location.md#L3` (Northern Region) | Addressed |

## C46

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C46` | What seven parameters were included in N-ISI, and how do neutrosophic classifications compare with classical SAR/Na%/RSC decision criteria? | Explicitly listed the seven N-ISI indices and added a classical-comparison disagreement summary and per-sample reporting (Table 24–Table 25). | `manuscript/sections/03-methodology/subsection-08-irrigation.md#L15` (compared against classical); `manuscript/sections/04-results/subsection-07-irrigation.md#L5` (Per-sample irrigation indices) | Addressed |

## C47

| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |
|---|---|---|---|---|
| `C47` | The CRMs listed (BCR 398/399) are typically for matrices other than groundwater major ions. Can you clarify the specific CRMs and check standards used for ion chromatography and cation measurements? | Removed specific CRM catalogue codes that cannot be verified from the provided materials and restated QA/QC more conservatively as use of appropriate standards and certified reference materials; detailed CRM identifiers remain unavailable. | `manuscript/sections/03-methodology/subsection-01-dataset.md#L9` (Representativeness and metadata limitations.) | Partially addressed |
