import numpy as np
import pandas as pd

def _calculate_entropy_weights(df: pd.DataFrame) -> pd.Series:
    """
    Calculates objective weights using the Information Entropy method.
    Higher variability in a parameter leads to a higher weight.
    """
    # 1. Min-max normalization
    df_norm = (df - df.min()) / (df.max() - df.min() + 1e-10)
    
    # 2. Calculate probability matrix
    p = df_norm / (df_norm.sum() + 1e-10)
    
    # 3. Calculate Entropy (e)
    n = len(df)
    k = 1.0 / np.log(n) if n > 1 else 1.0
    # Add small constant to log to avoid -inf
    e = -k * (p * np.log(p + 1e-10)).sum(axis=0)
    
    # 4. Degree of Diversification (d)
    d = 1.0 - e
    
    # 5. Weights (w)
    w = d / (d.sum() + 1e-10)
    return pd.Series(w, index=df.columns)

def calculate_irrigation_indices(data_meq: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates standard and advanced irrigation suitability indices using ion concentrations in meq/L.
    """
    indices = pd.DataFrame(index=data_meq.index)
    
    # 1. Sodium Adsorption Ratio (SAR)
    indices['SAR'] = data_meq['Na'] / np.sqrt((data_meq['Ca'] + data_meq['Mg']) / 2)
    
    # 2. Soluble Sodium Percentage (SSP)
    total_cations = data_meq['Ca'] + data_meq['Mg'] + data_meq['Na'] + data_meq['K']
    indices['SSP'] = (data_meq['Na'] + data_meq['K']) * 100 / total_cations
    
    # 3. Magnesium Hazard (MH)
    indices['MH'] = (data_meq['Mg'] * 100) / (data_meq['Ca'] + data_meq['Mg'])
    
    # 4. Kelly's Ratio (KR)
    indices['KR'] = data_meq['Na'] / (data_meq['Ca'] + data_meq['Mg'])
    
    # 5. Permeability Index (PI)
    indices['PI'] = (data_meq['Na'] + np.sqrt(data_meq['HCO3'])) * 100 / (data_meq['Ca'] + data_meq['Mg'] + data_meq['Na'])

    # 6. NEW: Residual Sodium Carbonate (RSC)
    # RSC = (HCO3 + CO3) - (Ca + Mg)
    co3 = data_meq['CO3'] if 'CO3' in data_meq.columns else 0
    indices['RSC'] = (data_meq['HCO3'] + co3) - (data_meq['Ca'] + data_meq['Mg'])
    
    # 7. NEW: Potential Salinity (PS)
    # PS = Cl + 0.5 * SO4
    indices['PS'] = data_meq['Cl'] + 0.5 * data_meq['SO4']
    
    # Classifications
    indices['SAR_Class'] = pd.cut(indices['SAR'], bins=[0, 10, 18, 26, np.inf], 
                                  labels=['Excellent', 'Good', 'Doubtful', 'Unsuitable'])
    indices['SSP_Class'] = pd.cut(indices['SSP'], bins=[0, 20, 40, 60, 80, 100], 
                                  labels=['Excellent', 'Good', 'Permissible', 'Doubtful', 'Unsuitable'])
    indices['MH_Class'] = np.where(indices['MH'] > 50, 'Unsuitable', 'Suitable')
    indices['KR_Class'] = np.where(indices['KR'] > 1, 'Unsuitable', 'Suitable')
    indices['PI_Class'] = pd.cut(indices['PI'], bins=[0, 25, 75, np.inf], 
                                 labels=['Unsuitable', 'Good', 'Excellent'])
    indices['RSC_Class'] = pd.cut(indices['RSC'], bins=[-np.inf, 1.25, 2.5, np.inf], 
                                  labels=['Safe', 'Marginal', 'Unsuitable'])
    
    return indices

def assess_irrigation_suitability(indices: pd.DataFrame) -> pd.DataFrame:
    """
    Advanced Neutrosophic Irrigation Suitability Index (N-ISI).
    Uses Entropy weights and Truth/Indeterminacy/Falsity logic.
    """
    # 1. Define thresholds for each index (S_j) based on FAO/USSL guidelines
    thresholds = {
        'SAR': 10.0, 'SSP': 60.0, 'MH': 50.0, 'KR': 1.0, 
        'PI': 75.0, 'RSC': 1.25, 'PS': 5.0
    }
    
    # 2. Calculate Normalized Hazard (Y)
    # Y = 1 means exactly at the threshold. Y > 1 is violation.
    Y = pd.DataFrame(index=indices.index)
    for col, S_j in thresholds.items():
        if col == 'PI': 
            # Permeability Index: Higher is Better. 
            # Standard Doneen classes: >75 (Excellent), 25-75 (Good), <25 (Unsuitable)
            # We map 75 -> Y=1, 100 -> Y=0.
            Y[col] = (100.0 - indices[col]) / (100.0 - S_j)
        else:
            # Other indices: Lower is Better.
            Y[col] = indices[col] / S_j
            
    # Clamp hazards at zero (negative hazard = perfect compliance)
    Y = Y.clip(lower=0.0)
            
    # 3. Neutrosophic Components
    # Truth (T): Compliance. Calibrated so T=0.5 at Y=1.
    T = np.exp(-0.693 * Y)
    
    # Falsity (F): Violation. Sigmoid centered at Y=1.
    F = 1 / (1 + np.exp(-10 * (Y - 1)))
    
    # Indeterminacy (I): Ambiguity near threshold.
    I = np.exp(-0.5 * ((Y - 1) / 0.2) ** 2)
    
    # 4. Objective Weights via Entropy (on Hazard Matrix Y)
    # We weight parameters based on their discriminatory hazard variance
    weights = _calculate_entropy_weights(Y)
    
    # 5. Aggregation
    # Suitability Score = Weighted sum of (1 + T - F) / 2
    # Result is 1.0 if perfectly compliant (T=1, F=0), 0.5 at threshold (T=0.5, F=0.5), 0.0 if failing.
    n_score = (1 + T - F) / 2
    isi_score = (n_score * weights).sum(axis=1) * 100
    
    # Bounding the score strictly between [0, 100]
    isi_score = isi_score.clip(lower=0.0, upper=100.0)
    
    results = pd.DataFrame({
        'N_ISI_Score': isi_score,
        'ISI_Class': pd.cut(isi_score, bins=[-np.inf, 50, 70, 85, np.inf], 
                            labels=['Unsuitable', 'Marginal', 'Good', 'Excellent']),
        'Indeterminacy_Mean': I.mean(axis=1),
        'Entropy_Weights': [weights.to_dict()] * len(indices)
    })
    
    return results
