# Reviewer Comments (verbatim)

Source: Comments on Karaga Manuscript.docx

Do not edit wording here; this file is the accuracy anchor for responses.

## C01

This manuscript tackles an important and impactful problem; groundwater quality and fluoride-related health risk in semi-arid sedimentary systems, using a thoughtfully integrated framework that addresses compositional closure, weight uncertainty in indexing, and probabilistic risk. The conceptual contribution is meaningful: the structural coupling of compositional process discovery with uncertainty-aware assessment is a step forward from stacked, uncoupled pipelines frequently seen in the literature. The case study is policy-relevant, and the probabilistic framing for vulnerable groups is valuable.

## C02

However, several elements require significant strengthening before this work meets the standard for publication in Scientific Reports. Chief among these are (i) limited sample size and single-season design without robust cross-validation or sensitivity analyses in the Bayesian layer; (ii) insufficient benchmarking of the proposed indices and models against established alternatives to demonstrate practical benefit; (iii) missing analytical QA/QC details crucial for data credibility; and (iv) incomplete mechanistic support for fluoride mobilization in the absence of speciation/saturation modeling. The writing is generally clear, but some essential implementation details are relegated to the Supplementary or remain unspecified in the main text.

## C03

Detailed comments are given below:

## C04

Technical soundness evaluation

## C05

The overall architecture is coherent and addresses real methodological shortcomings (closure, fixed weights, deterministic risk). However, the model evaluation layer is underdeveloped: lack of PSIS-LOO/WAIC, posterior predictive checks, sensitivity to priors, and external validation undermines confidence in the reported driver attributions and exceedance probabilities.

## C06

The horseshoe prior is appropriate, but with p not negligible relative to n (even after exclusions), model identifiability and stability should be scrutinized via cross-validation, shrinkage diagnostics (effective number of non-zero coefficients), and sensitivity to the global scale prior.

## C07

ILR-PCA with parallel analysis is suitable; back-interpretation in CLR space must be executed carefully (loadings in ILR do not trivially map to raw-ion effects). Please elaborate the back-transformation and interpretation pipeline to ensure process labels are not artefacts of the chosen SBP.

## C08

Experimental evaluation assessment

## C09

Field procedures: Pumping for ~5 minutes may be insufficient unless tied to well volume; please justify purging criteria. Report trip blanks, field duplicates, and any inter-lab comparison if available.

## C10

Analytical QA/QC: Provide explicit CBE acceptance threshold (±5% or ±10%), percentage of samples within range, and distribution of CBE. Report LOD/LOQ for all analytes and handling of any censored values.

## C11

Seasonal coverage: With only dry-season data, emphasize that exceedance probabilities are conditional on this season. Consider scenario analysis incorporating plausible wet-season dilution or draw on nearby longitudinal datasets to inform prior distributions.

## C12

Comparison with related work

## C13

The paper’s central novelty is integration rather than invention of new algorithms. Prior studies have applied compositional methods (ILR/CLR), fuzzy/neutrosophic indices, and Bayesian risk, but typically in isolation. This manuscript advances the field by structurally coupling these components; however, it would benefit from head-to-head comparisons with (i) deterministic WQIs under Monte Carlo weight perturbation, (ii) fuzzy vs neutrosophic indices, and (iii) standard regression vs horseshoe for driver profiling to quantify practical gains.

## C14

Expand on how this framework materially changes decisions compared to a conventional pipeline in the same study area (e.g., re-ranking of sites at the potability/irrigation interface, different remediation priorities).

## C15

Discussion of broader impact and significance

## C16

The demographic disaggregation in probabilistic risk assessment is valuable for policy targeting. The approach is potentially transferable across Voltaian and other semi-arid sedimentary contexts where fluoride risk is endemic.

## C17

The framework could serve as a template for regulators if accompanied by transparent, open-source code and minimal data requirements; otherwise, method complexity may hinder uptake. Clarify data and code availability, including the neutrosophic_pls package and Bayesian modeling scripts.

## C18

Reproducibility and transparency

## C19

Many critical steps are in the Supplementary Material or external repositories; for a methods-centric paper, more implementation specifics should be in the main text (e.g., neutrosophic membership functions, aesthetic penalty definition, bootstrapping scheme, PCA robustness settings, Bayesian sampler settings, convergence diagnostics: R-hat, ESS).

## C20

Provide an analysis workflow diagram and a data dictionary. Share input data (anonymized) with metadata to enable replication.

## C21

Data quality and sampling strategy

## C22

Clarify the representativeness of WVI-maintained public boreholes versus private wells/shallow hand-dug wells. Describe spatial coverage (per formation), well depths, screen intervals, and proximity to potential nitrate sources.

## C23

Report temperature and redox proxies (e.g., DO, Eh) if available; these are relevant to nitrogen species and carbonate equilibria.

## C24

Statistical methodology

## C25

For Bayesian models, include: number of chains, warmup and sampling iterations, adaptation parameters, R-hat and ESS, posterior predictive checks (PPC plots), LOO model diagnostics, and sensitivity to prior scales (e.g., alternative τ priors, regularized horseshoe).

## C26

For NPLS-GWQI, compare dynamic weights to expert-assigned weights and quantify variability across bootstraps (e.g., confidence intervals for N-VIP). Report how indeterminacy I and falsity F affect classification stability.

## C27

Hydrogeochemical interpretation and geochemical modeling

## C28

Interpretations of fluoride mobilization would be strengthened by speciation and saturation indices (calcite, dolomite, fluorite, fluorapatite), and ion activity corrections. If fluorite is undersaturated and calcite is near saturation, ion exchange and carbonate equilibria can be better disentangled.

## C29

Consider plotting Na-normalized or charge-balanced diagrams, and address whether halite dissolution is negligible given Na/Cl trends.

## C30

If “conductivity enrichment” is shorthand for evaporation/concentration, align terminology with Gibbs fields and mass-balance reasoning.

## C31

Risk assessment modeling

## C32

Specify exposure pathways modeled, parameter distributions used (ingestion rates, body weights, exposure frequency/duration), and their sources (local vs literature). Provide tornado plots and contribution-to-variance to highlight dominant uncertainty drivers.

## C33

Clarify the “skeletal-dental hazard threshold” definition (HQ>1 vs specific benchmarks), and report probabilities for multiple cutoffs. If nitrate dermal exposure is negligible, justify excluding it or show it is immaterial.

## C34

Irrigation suitability

## C35

List the seven parameters used, their standards or reference ranges, and the entropy weighting mechanics. Compare neutrosophic ISI classifications to classical SAR/Na%/RSC-based categorizations, including disagreement analysis.

## C36

Questions for Clarification

## C37

What explicit charge-balance error (CBE) threshold did you use, and what fraction of samples passed it? Please provide the CBE distribution.

## C38

How were detection limits (LOD/LOQ) for each analyte handled? Were any measurements censored, and if so, how were they treated in ILR and Bayesian analyses?

## C39

How did you determine that a ~5-minute purge was sufficient for each borehole? Can you report well volumes or stabilization criteria (pH/EC/temperature) used to confirm representativeness?

## C40

Please detail the ILR to CLR back-interpretation of loadings and how SBP choices affect process labeling. Did you test alternative SBPs to confirm stability of conclusions?

## C41

For the horseshoe models, please provide sampler diagnostics (R-hat, ESS), posterior predictive checks, PSIS-LOO or WAIC, and prior sensitivity analyses. How many predictors retained non-negligible mass, and were results robust to alternative shrinkage priors?

## C42

Can you benchmark NPLS-GWQI against at least one conventional WQI and a fuzzy index under Monte Carlo weight uncertainty to quantify classification stability and practical decision differences?

## C43

Which exposure pathways and parameter distributions were used for the hierarchical Bayesian risk model, and are they locally derived? Can you provide exceedance probabilities for multiple HQ thresholds and a contribution-to-variance breakdown?

## C44

Did you perform speciation/saturation modeling (e.g., PHREEQC) to test fluoride mobilization hypotheses (fluorite saturation, calcite/dolomite buffering)? If not, can you add this analysis or justify its omission?

## C45

Please clarify the regional affiliation of Karaga District (Northern Region vs North East Region) and correct any inconsistencies.

## C46

What seven parameters were included in N-ISI, and how do neutrosophic classifications compare with classical SAR/Na%/RSC decision criteria?

## C47

The CRMs listed (BCR 398/399) are typically for matrices other than groundwater major ions. Can you clarify the specific CRMs and check standards used for ion chromatography and cation measurements?
