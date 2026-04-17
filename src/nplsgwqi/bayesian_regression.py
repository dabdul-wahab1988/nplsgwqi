import numpy as np
import pandas as pd
import pymc as pm
from sklearn.preprocessing import StandardScaler

from .compositional import clr_transform

def _replace_nonpositive(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        pos = out.loc[out[col] > 0, col]
        floor = pos.min() / np.sqrt(2) if not pos.empty else 1e-6
        out.loc[out[col] <= 0, col] = floor
    return out

def prepare_clr_predictors(comp_data: pd.DataFrame, endpoint_name: str | None = None, include_scale: bool = True) -> pd.DataFrame:
    predictors = comp_data.copy()
    if endpoint_name in predictors.columns:
        predictors = predictors.drop(columns=[endpoint_name])
    predictors = _replace_nonpositive(predictors)
    
    clr_df = clr_transform(predictors)
    if include_scale:
        scale_val = np.mean(np.log(predictors), axis=1)
        clr_df["Total_Mineralization_Scale"] = scale_val
        
    return clr_df

def run_bayesian_endpoint_model(
    comp_data: pd.DataFrame,
    target: pd.Series,
    endpoint_name: str,
    *,
    include_scale: bool = True,
    use_augmented: bool = False,
    draws: int = 2000,
    tune: int = 1000,
    chains: int = 2,
    cores: int = 1,
    random_state: int = 42,
    log_likelihood: bool = True,
    posterior_predictive: bool = False,
    tau_beta: float = 1.0,
    lam_beta: float = 1.0,
    target_accept: float = 0.95,
    max_treedepth: int | None = None,
    likelihood: str = "normal",
    studentt_nu: float = 4.0,
    regularized_horseshoe: bool = False,
    slab_scale: float = 2.0,
    shrinkage_family: str = "cauchy",
) -> dict:
    if use_augmented:
        # We assume comp_data is ALREADY the augmented block from get_augmented_predictor_block
        x_df = comp_data.copy()
        block_label = "Augmented Field+Compositional"
    else:
        x_df = prepare_clr_predictors(
            comp_data, 
            endpoint_name=endpoint_name if endpoint_name in comp_data.columns else None,
            include_scale=include_scale
        )
        block_label = "Scale-Aware CLR chemistry"
        if endpoint_name in comp_data.columns:
            block_label = f"{endpoint_name}-safe {block_label}"
    
    # Scale predictors and target for robust sampling
    x_scaler = StandardScaler()
    X_scaled = x_scaler.fit_transform(x_df)
    
    y_scaler = StandardScaler()
    y_scaled = y_scaler.fit_transform(target.to_numpy().reshape(-1, 1)).ravel()

    features = list(x_df.columns)
    n_features = len(features)
    
    with pm.Model() as bayes_model:
        # Horseshoe / regularized-horseshoe prior for sparsity (non-centered parameterization).
        if shrinkage_family == "normal":
            tau = pm.HalfNormal("tau", sigma=tau_beta)
            lam = pm.HalfNormal("lam", sigma=lam_beta, shape=n_features)
        elif shrinkage_family == "cauchy":
            tau = pm.HalfCauchy("tau", beta=tau_beta)
            lam = pm.HalfCauchy("lam", beta=lam_beta, shape=n_features)
        else:
            raise ValueError("shrinkage_family must be 'cauchy' or 'normal'")

        beta_raw = pm.Normal("beta_raw", mu=0, sigma=1, shape=n_features)
        if regularized_horseshoe:
            # Regularized horseshoe slab (stabilizes heavy tails and reduces divergences under collinearity/outliers).
            c = pm.HalfNormal("slab_scale", sigma=slab_scale)
            lam_tilde = pm.math.sqrt((c**2 * lam**2) / (c**2 + (tau**2) * (lam**2)))
            beta = pm.Deterministic("beta", beta_raw * tau * lam_tilde)
        else:
            beta = pm.Deterministic("beta", beta_raw * tau * lam)
        alpha = pm.Normal('alpha', mu=0, sigma=1)
        
        # Expected value
        mu = alpha + pm.math.dot(X_scaled, beta)
        
        # Error term
        sigma = pm.HalfNormal('sigma', sigma=1)
        
        # Likelihood (robust option helps with influential points / high Pareto-k).
        if likelihood == "normal":
            Y_obs = pm.Normal('Y_obs', mu=mu, sigma=sigma, observed=y_scaled)
        elif likelihood == "studentt":
            Y_obs = pm.StudentT("Y_obs", nu=float(studentt_nu), mu=mu, sigma=sigma, observed=y_scaled)
        else:
            raise ValueError("likelihood must be 'normal' or 'studentt'")
        
        # Sample
        sample_kwargs = dict(
            draws=draws,
            tune=tune,
            chains=chains,
            target_accept=target_accept,
            random_seed=random_state,
            progressbar=False,
            cores=int(cores),
            return_inferencedata=True,
            idata_kwargs={"log_likelihood": bool(log_likelihood)},
        )
        if max_treedepth is not None:
            sample_kwargs["max_treedepth"] = int(max_treedepth)

        trace = pm.sample(**sample_kwargs)

        if posterior_predictive:
            trace = pm.sample_posterior_predictive(
                trace,
                var_names=["Y_obs"],
                random_seed=random_state,
                progressbar=False,
                extend_inferencedata=True,
            )

    fit_config = {
        "draws": int(draws),
        "tune": int(tune),
        "chains": int(chains),
        "cores": int(cores),
        "random_state": int(random_state),
        "target_accept": float(target_accept),
        "max_treedepth": int(max_treedepth) if max_treedepth is not None else None,
        "likelihood": str(likelihood),
        "studentt_nu": float(studentt_nu),
        "regularized_horseshoe": bool(regularized_horseshoe),
        "slab_scale": float(slab_scale),
        "shrinkage_family": str(shrinkage_family),
        "tau_beta": float(tau_beta),
        "lam_beta": float(lam_beta),
        "log_likelihood": bool(log_likelihood),
        "posterior_predictive": bool(posterior_predictive),
        "use_augmented": bool(use_augmented),
        "include_scale": bool(include_scale),
    }
    
    # Extract posteriors for coefficients
    beta_post = trace.posterior['beta'].values.reshape(-1, n_features)
    beta_mean = np.mean(beta_post, axis=0)
    beta_hdi_lower = np.percentile(beta_post, 2.5, axis=0)
    beta_hdi_upper = np.percentile(beta_post, 97.5, axis=0)
    
    # Probability of direction (pd) - how strictly positive/negative is the effect?
    pd_scores = np.maximum(np.mean(beta_post > 0, axis=0), np.mean(beta_post < 0, axis=0)) * 100.0
    
    # Identify "selected" features (e.g., pd > 95%)
    selected_mask = pd_scores > 95.0
    selected_features = [f for i, f in enumerate(features) if selected_mask[i]]
    
    # Map back to original scale (approximate for display/residuals)
    alpha_post = trace.posterior['alpha'].values.flatten()
    # X_scaled is (n_samples, n_features), beta_post.T is (n_features, n_draws)
    # The dot product is (n_samples, n_draws). alpha_post is (n_draws,)
    mu_post = alpha_post[None, :] + (X_scaled @ beta_post.T)
    y_pred_scaled = np.mean(mu_post, axis=1)
    y_pred = y_scaler.inverse_transform(y_pred_scaled.reshape(-1, 1)).ravel()
    
    # Calculate a pseudo-R2 (Bayesian R2 is complex, using standard R2 on posterior mean)
    from sklearn.metrics import r2_score
    r2_mean = r2_score(target, y_pred)
    
    # Calculate driver score (relative absolute mean effect)
    abs_mean = np.abs(beta_mean)
    driver_score = (abs_mean / np.sum(abs_mean)) * 100.0
    
    contrib = pd.DataFrame(X_scaled * beta_mean, columns=features, index=target.index)
    contrib["Fitted"] = y_pred
    contrib["Observed"] = target.values
    contrib["Residual"] = target.values - y_pred
    
    return {
        "endpoint": endpoint_name,
        "transform": "Bayesian_CLR+Scale",
        "predictor_columns": features,
        "predictor_block_label": block_label,
        "fit_config": fit_config,
        "selected_features": selected_features,
        "in_sample_r2": float(r2_mean),
        "beta_mean": pd.Series(beta_mean, index=features),
        "beta_hdi_lower": pd.Series(beta_hdi_lower, index=features),
        "beta_hdi_upper": pd.Series(beta_hdi_upper, index=features),
        "prob_direction": pd.Series(pd_scores, index=features),
        "driver_score": pd.Series(driver_score, index=features).sort_values(ascending=False),
        "sample_level_contributions": contrib,
        "fitted_values": pd.Series(y_pred, index=target.index),
        "trace": trace
    }
