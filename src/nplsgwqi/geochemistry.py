import numpy as np
import pandas as pd

def calculate_geochemical_constraints(data_meq: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates deterministic hydrogeochemical indices to constrain process interpretation.
    """
    constraints = pd.DataFrame(index=data_meq.index)
    
    # 1. Total Dissolved Cations (TZ+) and Anions (TZ-)
    cat_cols = [c for c in ['Ca', 'Mg', 'Na', 'K'] if c in data_meq.columns]
    ani_cols = [c for c in ['Cl', 'SO4', 'HCO3', 'NO3', 'CO3', 'F'] if c in data_meq.columns]
    tz_plus = data_meq[cat_cols].sum(axis=1).replace(0, 1e-6)
    tz_minus = data_meq[ani_cols].sum(axis=1).replace(0, 1e-6)
    
    # 2. Silicate vs Carbonate Weathering Proxies
    # (Na + K) / TZ+ -> If near 1, Silicate weathering is dominant
    na_k_sum = data_meq[[c for c in ['Na', 'K'] if c in data_meq.columns]].sum(axis=1)
    constraints['Silicate_Index'] = na_k_sum / tz_plus
    
    # (Ca + Mg) / HCO3 -> If near 0.5, Carbonate weathering is dominant
    ca_mg_sum = data_meq[[c for c in ['Ca', 'Mg'] if c in data_meq.columns]].sum(axis=1)
    constraints['Carbonate_Index'] = ca_mg_sum / data_meq['HCO3']
    
    # 3. Na/Cl Ratio -> If > 1, Albite weathering is likely; if approx 1, Halite dissolution
    constraints['Na_Cl_Ratio'] = data_meq['Na'] / data_meq['Cl']
    
    # 4. Ion Exchange Indices (Schoeller CAI)
    # CAI 1 = [Cl - (Na + K)] / Cl
    cl_minus_nak = (data_meq['Cl'] - na_k_sum)
    constraints['CAI_1'] = cl_minus_nak / data_meq['Cl']
    # CAI 2 = [Cl - (Na + K)] / (SO4 + HCO3 + CO3 + NO3)
    cai2_denom_cols = [c for c in ['SO4', 'HCO3', 'CO3', 'NO3'] if c in data_meq.columns]
    denom = data_meq[cai2_denom_cols].sum(axis=1).replace(0, 1e-6)
    constraints['CAI_2'] = cl_minus_nak / denom
    
    # 5. Cation Exchange Value (CEV)
    # (Ca + Mg) - (HCO3 + SO4) vs (Na - Cl)
    constraints['Exchange_X'] = (data_meq['Na'] - data_meq['Cl'])
    constraints['Exchange_Y'] = (data_meq['Ca'] + data_meq['Mg']) - (data_meq['HCO3'] + data_meq['SO4'])
    
    # 6. Gibbs Classification
    # We calculate TDS_calc by converting meq/L back to mg/L using exact equivalent weights
    # Equivalent Weight = Atomic Weight / Valence
    from .preprocessing import Preprocessor
    weights = Preprocessor.STOI_VALENCE_WEIGHTS
    
    tds_components = pd.DataFrame(index=data_meq.index)
    for ion, props in weights.items():
        if ion in data_meq.columns:
            tds_components[ion] = data_meq[ion] * props['eq_weight']
            
    tds_calc = tds_components.sum(axis=1)
    na_ratio = data_meq['Na'] / (data_meq['Na'] + data_meq['Ca'])
    
    conditions = [
        (tds_calc < 50),
        (tds_calc >= 50) & (tds_calc <= 1000),
        (tds_calc > 1000)
    ]
    choices = ['Precipitation', 'Rock-Weathering', 'Evaporation']
    constraints['Gibbs_Dominance'] = np.select(conditions, choices, default='Mixed')
    constraints['TDS_Calculated'] = tds_calc
    
    return constraints

def audit_process_names(loadings_clr: pd.DataFrame, constraints: pd.DataFrame) -> dict:
    """
    Uses geochemical constraints to assign conservative hydrochemical signature labels.
    """
    refined_names = {}
    
    # Use robust central tendency because ratio-style indices can be highly skewed.
    med_carbonate = constraints["Carbonate_Index"].replace([np.inf, -np.inf], np.nan).median()
    med_na_cl = constraints["Na_Cl_Ratio"].replace([np.inf, -np.inf], np.nan).median()
    med_cai = constraints["CAI_1"].replace([np.inf, -np.inf], np.nan).median()
    
    for proc in loadings_clr.index:
        loadings = loadings_clr.loc[proc]
        top_pos = loadings.nlargest(2).index.tolist()
        
        # Conservative interpretation logic: assign signatures/associations rather
        # than definitive mineral-reaction labels unless the evidence is unusually strong.
        if "Na" in top_pos and "HCO3" in top_pos:
            if med_cai < 0:
                refined_names[proc] = "Na-HCO3 Exchange Signature"
            else:
                refined_names[proc] = "Na-HCO3 Alkaline Signature"

        elif "Na" in top_pos and "Cl" in top_pos:
            if med_na_cl > 1.2:
                refined_names[proc] = "Na-Cl Enrichment Signature"
            elif 0.8 <= med_na_cl <= 1.2:
                refined_names[proc] = "Na-Cl Salinity Signature"
            else:
                refined_names[proc] = "Cl-Dominant Salinity Signature"

        elif "Ca" in top_pos and "HCO3" in top_pos:
            if 0.4 <= med_carbonate <= 0.6:
                refined_names[proc] = "Ca-HCO3 Carbonate Signature"
            else:
                refined_names[proc] = "Ca-HCO3 Weathering Signature"

        else:
            refined_names[proc] = f"{'-'.join(top_pos)} Association"
            
    return refined_names
