## 3.6 Bayesian Endpoint-Driver Modelling

The relative influence of hydrochemical predictors on NPLS-GWQI and on fluoride concentration was quantified using Bayesian horseshoe regression implemented in PyMC. Let $y_{s}^{*}$ denote the standardized endpoint value for sample $s$ and $x_{s}^{*}$ the corresponding standardized predictor vector from the augmented hydrochemical block. The model was specified as

$$
y_{s}^{*} \sim \mathcal{N}\left( \mu_{s},\sigma \right),\quad\quad\mu_{s} = \alpha + {\mathbf{x}_{\mathbf{s}}^{*}}^{\top}\mathbf{\beta}
$$

with priors

$$
\alpha \sim \mathcal{N}(0,1),\quad\quad\sigma \sim \text{HalfNormal}(1),
$$

$$
\beta_{j} \mid \lambda_{j},\tau \sim \mathcal{N}\left( 0,\tau^{2}\lambda_{j}^{2} \right),\quad\quad\lambda_{j} \sim \text{Half-Cauchy}(0,1),\quad\quad\tau \sim \text{Half-Cauchy}(0,1).
$$

The horseshoe prior was selected because its global-local shrinkage structure pulls weak, noise-like coefficients towards zero while allowing a small number of informative predictors to retain substantial posterior mass, which is appropriate given the size of the augmented predictor block relative to the 34-sample dataset <sup>37</sup>.

For both endpoint models, predictors were drawn from an augmented hydrochemical block comprising log-transformed major-ion concentrations, field variables (pH, EC, and TDS), linear ion-exchange and weathering proxies ($CAI_{1}$, $CAI_{2}$, Silicate Proxy, and Carbonate Proxy), and log-transformed geochemical ratios (Na/Cl, Ca/HCO<sub>3</sub>, and Mg/Ca). For the fluoride endpoint model, fluoride itself was excluded from the major-ion predictor block.

Posterior sampling was performed with the NUTS algorithm in PyMC. Model outputs included posterior mean coefficients, 95% equal-tail credible intervals, and the probability of direction for each predictor, defined as

$$
\text{pd}\left( \beta_{j} \right) = \max\left\{ \Pr\left( \beta_{j} > 0\mid y \right),\ \ \Pr\left( \beta_{j} < 0\mid y \right) \right\} \times 100.
$$

This quantity represents the posterior probability that the coefficient has the reported sign. Normalised driver contribution percentages were computed from the absolute posterior mean coefficients, and predictors with pd \> 95% were treated as directionally supported. Model fit was summarised using an in-sample pseudo-$R^{2}$ computed from posterior-mean fitted values transformed back to the original response scale. A supplementary sparse PLS model was fitted for fluoride as an external comparator, and its cross-validation and permutation statistics are reported alongside the Bayesian results.

**Sampler configuration and diagnostics.** To address identifiability concerns when $p$ is not negligible relative to $n$, NUTS sampling was performed using multiple chains with extended warmup and a conservative target acceptance rate. Convergence and effective-sample-size diagnostics (R-hat and ESS) were computed for all regression parameters, and posterior predictive checks (PPC) were used to assess whether replicated data drawn from the posterior reproduce the observed endpoint distribution. Out-of-sample predictive adequacy was approximated using PSIS-LOO and WAIC computed from the model log-likelihood. A prior-sensitivity analysis was conducted by repeating the fits under alternative global shrinkage scales (τ prior scale) and quantifying the stability of ranked driver attributions across priors (Table 20). Full diagnostics, including PSIS-LOO Pareto-$k$ summaries and PPC plots, are reported in Table 19 and Figure S10, and are used to qualify driver interpretation in Section 4.5.
