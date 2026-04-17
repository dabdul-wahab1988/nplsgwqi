import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
from typing import Dict, Any, List
from .wqi import calculate_npls_gwqi
from .irrigation import calculate_irrigation_indices, assess_irrigation_suitability

def run_probabilistic_sensitivity(
    data: pd.DataFrame, 
    standards: pd.Series, 
    n_samples: int = 500, 
    error_pct: float = 0.05,
    random_seed: int = 42
) -> Dict[str, pd.DataFrame]:
    """
    Performs Bayesian Measurement Error Simulation for N-WQI and N-ISI.
    
    This simulates 5% analytical error in laboratory measurements and calculates
    the probability of each sample belonging to a specific quality class.
    """
    # 1. Select parameters involved in WQI and Irrigation
    # Note: We focus on the numeric parameters used in calculations
    calc_vars = list(set(standards.index.tolist() + ['Ca', 'Mg', 'Na', 'K', 'Cl', 'SO4', 'HCO3', 'pH']))
    # Filter for variables actually present in data
    calc_vars = [v for v in calc_vars if v in data.columns]
    
    obs_data = data[calc_vars].values
    n_obs, n_feats = obs_data.shape
    
    print(f"Starting Bayesian Sensitivity Analysis (Analytical Error = {error_pct*100}%)...")
    
    # 2. PyMC Measurement Error Model
    # We model the 'True' concentration as a latent variable
    with pm.Model() as model:
        # Prior for True Concentration (Truncated Normal to keep positive)
        # We use a relatively informed prior centered at the observation
        true_conc = pm.TruncatedNormal(
            "true_conc", 
            mu=obs_data, 
            sigma=obs_data * error_pct, 
            lower=0, 
            shape=(n_obs, n_feats)
        )
        
        # Likelihood (The observation is a noisy version of the true value)
        pm.Normal("obs", mu=true_conc, sigma=obs_data * error_pct, observed=obs_data)
        
        # Sampling
        trace = pm.sample(n_samples, tune=500, chains=2, random_seed=random_seed, progressbar=False)
        
    # Extract posterior samples of the 'True' concentrations
    # Shape: (draws, samples, features)
    posterior_draws = trace.posterior["true_conc"].values
    posterior_draws = posterior_draws.reshape(-1, n_obs, n_feats)
    
    total_draws = posterior_draws.shape[0]
    
    # 3. Propagate Uncertainty
    wqi_classes_matrix = []
    isi_classes_matrix = []
    
    # Stoichiometric conversion factors (mg/L -> meq/L)
    STOI_FACTORS = {
        'Ca': 2 / 40.078, 'Mg': 2 / 24.305, 'Na': 1 / 22.990, 'K': 1 / 39.098,
        'Cl': 1 / 35.45, 'HCO3': 1 / 61.016, 'SO4': 2 / 96.06, 'NO3': 1 / 62.004,
        'CO3': 2 / 60.01, 'F': 1 / 18.998
    }
    
    # To save time, we calculate indices for each draw
    for i in range(total_draws):
        # original mg/L data for WQI
        draw_df = pd.DataFrame(posterior_draws[i], columns=calc_vars, index=data.index)
        
        # Calculate N-WQI (requires mg/L)
        wqi_res, _ = calculate_npls_gwqi(draw_df[standards.index], standards, n_bootstraps=1)
        wqi_classes_matrix.append(wqi_res['Quality_Class'].values)
        
        # Calculate N-ISI (requires meq/L)
        try:
            # Convert draw_df to meq for irrigation
            draw_meq = draw_df.copy()
            for ion, factor in STOI_FACTORS.items():
                if ion in draw_meq.columns:
                    draw_meq[ion] *= factor
            
            irr_indices = calculate_irrigation_indices(draw_meq)
            isi_res = assess_irrigation_suitability(irr_indices)
            isi_classes_matrix.append(isi_res['ISI_Class'].values)
        except Exception:
            continue

    # 4. Aggregate Probabilities
    def calculate_probs(class_matrix, classes_order):
        df_classes = pd.DataFrame(class_matrix)
        probs_list = []
        for col in df_classes.columns:
            counts = df_classes[col].value_counts(normalize=True)
            probs = {c: counts.get(c, 0.0) * 100 for c in classes_order}
            # Add most probable class
            probs['Most_Probable_Class'] = counts.idxmax()
            probs['Confidence_Pct'] = counts.max() * 100
            probs_list.append(probs)
        return pd.DataFrame(probs_list, index=data.index)

    wqi_order = ['Excellent', 'Good', 'Poor', 'Very Poor', 'Unsuitable']
    isi_order = ['Excellent', 'Good', 'Marginal', 'Unsuitable']
    
    wqi_probs = calculate_probs(wqi_classes_matrix, wqi_order)
    isi_probs = calculate_probs(isi_classes_matrix, isi_order)
    
    return {
        'WQI_Probabilities': wqi_probs,
        'ISI_Probabilities': isi_probs
    }
