## 3.8 Irrigation Suitability Assessment

Irrigation suitability was evaluated using seven standard hydrochemical indices computed in meq/L space: the sodium adsorption ratio (SAR), the soluble sodium percentage (SSP), the magnesium hazard (MH), the Kelly’s ratio (KR), the permeability index (PI), the residual sodium carbonate (RSC), and the potential salinity (PS). These indices collectively characterise the sodicity hazard, permeability risk, residual alkalinity, and generalised salinity impact of groundwater on soils and crops. Details for the formulation of the Irrigation Suitability Index are reported in the Supplementary methodology S3.

Objective weights for each irrigation index were derived using the entropy-weighting method, which assigns greater weight to indices showing larger relative dispersion in the dataset, thereby avoiding the subjectivity associated with expert-assigned fixed weights <sup>38,39</sup>.

Neutrosophic irrigation suitability indexing (N-ISI) was then applied to aggregate the seven entropy-weighted indices into a composite suitability score using the same neutrosophic truth-indeterminacy-falsity architecture as outline in the Supplementary methodology S3:

$$
\text{N-ISI}_{s} = 100 \cdot \sum_{j = 1}^{m}{w_{j}\eta_{j,s}}
$$

where $w_{j}$ is the normalised entropy weight, $\eta_{j,s}$ is the per-index suitability score derived from the truth and falsity components. N-ISI classification assigned each sample to one of four irrigation suitability classes: Excellent, Good, Marginal, or Unsuitable (see Supplementary methodology S3.6).

To improve transparency and enable comparison with established decision frameworks, N-ISI classifications were compared against classical SAR, SSP, and RSC-based categorisations, and an overall disagreement summary was computed (Table 25). The per-sample irrigation indices underpinning both the neutrosophic and classical classifications are reported in Table 24.
