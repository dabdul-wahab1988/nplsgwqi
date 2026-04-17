import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression

def _calculate_nvip(pls_model, x_tif: np.ndarray, x_combined: pd.DataFrame) -> tuple:
    """
    Computes exact channel-decomposed Neutrosophic VIP.
    Because sklearn's PLSRegression auto-scales the input data by default,
    we must securely scale our x_tif projection matrices using the 
    exact mean and scale learned by the model from x_combined.
    """
    # Model parameters
    w = pls_model.x_weights_  # (n_features, n_components)
    q = pls_model.y_loadings_ # (n_components, n_targets)
    
    # Sklearn scaling factors from the fit
    x_mean = pls_model._x_mean
    x_std = pls_model._x_std
    
    # Extract channels
    T = x_tif[:, :, 0]
    I = x_tif[:, :, 1]
    F = x_tif[:, :, 2]
    
    n_features = x_tif.shape[1]
    
    # We must center/scale the individual channels similarly to how x_combined was scaled.
    # Note: X_combined = T + I + F. If we scale T by x_std, then (T+I+F)/x_std = X_scaled
    # Wait, T_c + I_c + F_c must equal the scaled X_combined. 
    # Let's allocate the mean: T_c = (T - T.mean) / x_std
    T_c = (T - T.mean(axis=0)) / x_std
    I_c = (I - I.mean(axis=0)) / x_std
    F_c = (F - F.mean(axis=0)) / x_std
    
    # Projections onto the learned PLS weights
    t_T = T_c @ w  # (n_samples, n_components)
    t_I = I_c @ w
    t_F = F_c @ w
    
    q_sq = np.sum(q**2, axis=0)  # (n_components,)
    
    # Channel Sums of Squares
    ss_T = np.sum(t_T**2, axis=0) * q_sq
    ss_I = np.sum(t_I**2, axis=0) * q_sq
    ss_F = np.sum(t_F**2, axis=0) * q_sq
    
    total_ss = np.sum(ss_T + ss_I + ss_F)
    if total_ss < 1e-12:
        total_ss = 1e-12
        
    w_sq = w ** 2
    
    vip_t_sq = n_features * (w_sq @ ss_T) / total_ss
    vip_i_sq = n_features * (w_sq @ ss_I) / total_ss
    vip_f_sq = n_features * (w_sq @ ss_F) / total_ss
    
    vip_t = np.sqrt(np.maximum(vip_t_sq, 0))
    vip_i = np.sqrt(np.maximum(vip_i_sq, 0))
    vip_f = np.sqrt(np.maximum(vip_f_sq, 0))
    
    # Aggregate VIP using the L2 norm
    aggregate_vip = np.sqrt(vip_t**2 + vip_i**2 + vip_f**2)
    return aggregate_vip, vip_t, vip_i, vip_f

def calculate_npls_gwqi(data: pd.DataFrame, standards: pd.Series, n_bootstraps: int = 50) -> tuple:
    """
    Advanced Neutrosophic PLS-based Groundwater Quality Index.
    
    Improves upon standard PLS-WQI by:
    1. Using Bipolar V-Shaped Normalization for pH (acidity and alkalinity).
    2. Using Neutrosophic Indeterminacy (I) for reliability-based weighting.
    3. Using Truth (T) and Falsity (F) to refine the sub-index scores.
    4. Applying Bootstrapping to ensure VIP weight stability.
    """
    # 1. Normalize by standards
    # Default is simple division (lower is better)
    Y = data.div(standards, axis=1)
    
    # UPGRADE: Bipolar pH Normalization
    # pH hazard is zero at 7.0 and 1.0 at thresholds (6.5 or 8.5)
    if 'pH' in data.columns:
        ph_vals = data['pH']
        # Acidic hazard: (7.0 - pH) / (7.0 - 6.5)
        acid_haz = (7.0 - ph_vals) / 0.5
        # Alkaline hazard: (pH - 7.0) / (8.5 - 7.0)
        alk_haz = (ph_vals - 7.0) / 1.5
        # Combine: take the max hazard (bipolar)
        Y['pH'] = np.maximum(acid_haz, alk_haz).clip(lower=0.0)

    # 2. Neutrosophic logic components
    # Truth (T): Degree of compliance (High when Y << 1)
    T = np.exp(-1.5 * Y)
    
    # Falsity (F): Degree of violation (High when Y > 1)
    F = 1 / (1 + np.exp(-10 * (Y - 1)))
    
    # NEW: Aesthetic/Sensory Penalty (TDS, Fe, etc.)
    # Note: pH penalty is now handled mathematically in Y['pH']
    sensory_penalty = pd.Series(1.0, index=data.index)
    if 'TDS' in data.columns:
        # Penalize if TDS > 500 (Aesthetic threshold)
        tds_penalty = np.where(data['TDS'] > 500, 1.1, 1.0)
        sensory_penalty *= tds_penalty
    if 'Fe' in data.columns:
        # Fe > 0.3 causes staining and metallic taste
        fe_penalty = np.where(data['Fe'] > 0.3, 1.25, 1.0)
        sensory_penalty *= fe_penalty

    # Indeterminacy (I): Ambiguity near threshold (Gaussian at Y=1)
    I = np.exp(-0.5 * ((Y - 1) / 0.2) ** 2)
    
    # 3. Create True Neutrosophic Channels
    # We modulate the raw physical concentration (data) by its neutrosophic states.
    # This ensures the PLS variance represents the physical quantity's structural impact.
    X_T = data * T
    X_I = data * I
    X_F = data * F
    
    # Create the 3D x_tif multi-channel matrix and the unified predictor
    x_tif = np.stack([X_T.values, X_I.values, X_F.values], axis=-1)
    X_combined = X_T + X_I + X_F
    
    # Bootstrapped N-VIP scores for stability
    vip_accumulator = []
    n_samples = data.shape[0]
    
    for _ in range(n_bootstraps):
        boot_idx = np.random.choice(np.arange(n_samples), size=n_samples, replace=True)
        X_b = X_combined.iloc[boot_idx]
        Y_b = Y.iloc[boot_idx]
        x_tif_b = x_tif[boot_idx]
        
        try:
            pls = PLSRegression(n_components=min(3, data.shape[1]), scale=True)
            pls.fit(X_b, Y_b)
            n_vip_agg, _, _, _ = _calculate_nvip(pls, x_tif_b, X_b)
            vip_accumulator.append(n_vip_agg)
        except:
            continue
            
    if not vip_accumulator:
        # Fallback to single fit if bootstrap fails
        pls = PLSRegression(n_components=min(3, data.shape[1]), scale=True)
        pls.fit(X_combined, Y)
        n_vip_agg, _, _, _ = _calculate_nvip(pls, x_tif, X_combined)
        mean_vips = n_vip_agg
    else:
        mean_vips = np.mean(vip_accumulator, axis=0)
        
    dynamic_weights = mean_vips / np.sum(mean_vips)
    
    # 4. Neutrosophic Aggregation
    # Instead of raw Y * 100, we use (1 + F - T) as a penalty-aware scalar
    # This ensures that high violations (F -> 1, T -> 0) weigh more than simple ratios
    neutro_scalar = (1 + F - T)
    sub_indices = Y * 100 * neutro_scalar
    
    gwqi = sub_indices.dot(dynamic_weights) * sensory_penalty
    
    results = pd.DataFrame({
        'NPLS_GWQI': gwqi,
        'Quality_Class': pd.cut(gwqi, bins=[0, 50, 100, 200, 300, np.inf], 
                                labels=['Excellent', 'Good', 'Poor', 'Very Poor', 'Unsuitable']),
        'Truth_Mean': T.mean(axis=1),
        'Indeterminacy_Mean': I.mean(axis=1),
        'Falsity_Mean': F.mean(axis=1)
    })
    
    return results, pd.Series(dynamic_weights, index=data.columns, name='Dynamic_Weights')
