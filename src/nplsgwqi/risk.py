import numpy as np
import pandas as pd
import pymc as pm

def calculate_sample_risk_deterministic(data: pd.DataFrame, toxref: dict, bw_mean: np.ndarray, ir_mean: np.ndarray, ed_mean: np.ndarray, sa_mean: np.ndarray, et_mean: np.ndarray) -> dict:
    """Calculates deterministic sample-level risk estimates."""
    contaminants = [c for c in data.columns if c in toxref and toxref[c]['RfD_oral'] != np.inf]
    sample_risk = {}
    
    # We will use Adults (index 0) for sample-level base reporting by default
    bw = bw_mean[0]
    ir = ir_mean[0]
    ed = ed_mean[0]
    sa = sa_mean[0]
    et = et_mean[0]
    
    ef = 365.0
    at = ed * 365.0
    
    f_ing = (ir * ef * ed) / (at * bw)
    s_derm = (et * sa * 1.0 * ef * ed) / (at * bw)
    
    for c in contaminants:
        tref = toxref[c]
        c_val = data[c]
        
        hq_ing = (c_val * f_ing) / tref['RfD_oral']
        if tref['RfD_derm'] != np.inf and tref['Kp'] > 0:
            hq_der = (c_val * tref['Kp'] * 1e-3 / tref['RfD_derm']) * s_derm
            hq_total = hq_ing + hq_der
        else:
            hq_total = hq_ing
            
        sample_risk[f'HQ_sample_{c}'] = hq_total
            
    # NEW: Aggregating into Organ-Specific HI for deterministic samples
    organ_sets = set(org for c in contaminants for org in toxref[c]['organs'])
    for org in organ_sets:
        hi_org = pd.Series(0.0, index=data.index)
        for c in contaminants:
            if org in toxref[c].get('organs', []):
                hi_org += sample_risk[f'HQ_sample_{c}']
        sample_risk[f'HI_sample_{org}'] = hi_org
        
    return sample_risk
def assess_health_risk(
    data: pd.DataFrame,
    toxref: dict = None,
    n_samples: int = 2000,
    mode: str = 'both',
    random_seed: int | None = None,
    chains: int = 2,
    tune: int = 500,
    target_accept: float = 0.9,
) -> dict:
    """
    Hierarchical Bayesian Multivariate Probabilistic Risk Assessment (HBMPRA)
    and Deterministic Sample-Level Risk.
    """
    if toxref is None:
        toxref = {
            'F': {'RfD_oral': 0.06, 'RfD_derm': np.inf, 'Kp': 0.0, 'organs': ['skeletal_dental']},
            'NO3': {'RfD_oral': 1.6, 'RfD_derm': np.inf, 'Kp': 0.0, 'organs': ['hemato']},
            'Pb': {'RfD_oral': 3.5e-3, 'RfD_derm': 5.25e-4, 'Kp': 1e-4, 'organs': ['neuro', 'nephro', 'hemato']},
            'As': {'RfD_oral': 3.0e-4, 'RfD_derm': 1.23e-4, 'Kp': 1e-3, 'organs': ['hepato', 'cardiovascular', 'derm']},
            'Ca': {'RfD_oral': np.inf, 'RfD_derm': np.inf, 'Kp': 0.0, 'organs': []},
            'pH': {'RfD_oral': np.inf, 'RfD_derm': np.inf, 'Kp': 0.0, 'organs': []},
        }

    groups = ["Adults", "Children", "Teens"]
    n_groups = len(groups)
    
    bw_mean = np.array([60.0, 15.0, 45.0])
    ir_mean = np.array([2.5, 1.0, 1.5])
    ed_mean = np.array([30.0, 6.0, 16.0])
    sa_mean = np.array([18000.0, 6600.0, 12000.0])
    et_mean = np.array([0.58, 1.0, 0.7])
    
    contaminants = [c for c in data.columns if c in toxref and toxref[c]['RfD_oral'] != np.inf]
    if not contaminants:
        return {}
        
    # NEW: Speciation Logic (Probabilistic)
    # If As is present, we model the As(III) fraction
    # As(III) is ~60x more toxic. Default RfD in toxref is often for As(total) or As(V).
    # We will split As into As3 and As5 for the risk calculation.
    if 'As' in contaminants:
        as_idx = contaminants.index('As')
        contaminants.pop(as_idx)
        contaminants.extend(['As_III', 'As_V'])
        
        # Define specific RfDs for species if not already in toxref
        toxref['As_III'] = {'RfD_oral': 3.0e-4, 'RfD_derm': 1.23e-4, 'Kp': 1e-3, 'organs': ['hepato', 'cardiovascular', 'derm']}
        toxref['As_V'] = {'RfD_oral': 1.0e-3, 'RfD_derm': 4.1e-4, 'Kp': 1e-3, 'organs': ['hepato', 'cardiovascular', 'derm']}

    result = {}
    
    if mode in ['deterministic', 'both']:
        # For deterministic, assume 20% As(III) as a conservative estimate if no Eh
        as_3_ratio = 0.2
        data_det = data.copy()
        if 'As' in data.columns:
            data_det['As_III'] = data['As'] * as_3_ratio
            data_det['As_V'] = data['As'] * (1 - as_3_ratio)
            
        sample_risk = calculate_sample_risk_deterministic(data_det, toxref, bw_mean, ir_mean, ed_mean, sa_mean, et_mean)
        result.update(sample_risk)
        
    if mode in ['bayesian', 'both']:
        with pm.Model() as risk_model:
            CV_BW = 0.21
            sigma_log_bw = np.sqrt(np.log(1.0 + CV_BW**2))
            mu_log_bw = np.log(bw_mean) - 0.5 * np.log(1.0 + CV_BW**2)
            z_log_bw = pm.Normal("z_log_bw", mu=0.0, sigma=1.0, shape=n_groups)
            BW_g = pm.Deterministic("BW_g", pm.math.exp(mu_log_bw + z_log_bw * sigma_log_bw))
            
            IR_perkg_med = ir_mean / bw_mean 
            sigma_log_ir = 0.6
            mu_log_ir = np.log(IR_perkg_med)
            z_log_ir = pm.Normal("z_log_ir", mu=0.0, sigma=1.0, shape=n_groups)
            IR_perkg_g = pm.Deterministic("IR_perkg_g", pm.math.exp(mu_log_ir + z_log_ir * sigma_log_ir))
            
            EF_days = pm.Data("EF_days", np.array([365.0, 365.0, 365.0]))
            EF_v = EF_days / 365.0
            
            ED_v = pm.Data("ED_v", ed_mean)
            AT_nc = ED_v * 365.0
            
            SA_v = pm.Data("SA_v", sa_mean)
            
            hq_dict = {}
            for c in contaminants:
                tref = toxref[c]
                
                # Special handling for split species in Bayesian model
                if c in ['As_III', 'As_V']:
                    raw_c = 'As'
                    as_total = pm.TruncatedNormal(f"C_{raw_c}", mu=np.mean(data[raw_c]), sigma=np.std(data[raw_c]) + 1e-6, lower=0.0)
                    
                    # Probabilistic Speciation: Use a Beta prior for the fraction
                    # Alpha=2, Beta=8 gives a mean of 0.2 (20% As III) with uncertainty
                    as3_frac = pm.Beta("as3_frac", alpha=2, beta=8)
                    if c == 'As_III':
                        C = as_total * as3_frac
                    else:
                        C = as_total * (1 - as3_frac)
                else:
                    mu_c = np.mean(data[c])
                    sigma_c = np.std(data[c]) + 1e-6
                    C = pm.TruncatedNormal(f"C_{c}", mu=mu_c, sigma=sigma_c, lower=0.0)
                
                F_ing_nc = (IR_perkg_g * EF_v * ED_v) / (AT_nc / 365.0)
                CDI_ing = pm.Deterministic(f"CDI_ing_{c}", C * F_ing_nc)
                HQ_ing = CDI_ing / tref['RfD_oral']
                
                if tref['RfD_derm'] != np.inf and tref['Kp'] > 0:
                    S_group = (et_mean * SA_v * 1.0 * EF_days * ED_v) / (AT_nc * BW_g)
                    HQ_der = (C * tref['Kp'] * 1e-3 / tref['RfD_derm']) * S_group
                    HQ_total = pm.Deterministic(f"HQ_{c}", HQ_ing + HQ_der)
                else:
                    HQ_total = pm.Deterministic(f"HQ_{c}", HQ_ing)
                    
                hq_dict[c] = HQ_total
                
            organ_sets = set(org for c in contaminants for org in toxref[c]['organs'])
            for org in organ_sets:
                org_hqs = [hq_dict[c] for c in contaminants if org in toxref[c]['organs']]
                pm.Deterministic(f"HI_{org}", pm.math.sum(org_hqs, axis=0))
                
            trace = pm.sample(
                draws=n_samples,
                tune=tune,
                chains=chains,
                cores=1,
                target_accept=target_accept,
                random_seed=random_seed,
                progressbar=False,
                return_inferencedata=True,
                idata_kwargs={"log_likelihood": False},
            )
        
        for i, g in enumerate(groups):
            for c in contaminants:
                hq_post = trace.posterior[f"HQ_{c}"].values.reshape(-1, n_groups)[:, i]
                result[f'HQ_{g}_{c}_mean'] = float(np.mean(hq_post))
                result[f'HQ_{g}_{c}_posterior'] = hq_post
                
            for org in organ_sets:
                hi_post = trace.posterior[f"HI_{org}"].values.reshape(-1, n_groups)[:, i]
                result[f'HI_{g}_{org}_mean'] = float(np.mean(hi_post))
                result[f'HI_{g}_{org}_95_CI'] = [float(np.percentile(hi_post, 2.5)), float(np.percentile(hi_post, 97.5))]
                result[f'HI_{g}_{org}_prob_exceed_1'] = float(np.mean(hi_post > 1.0))
                result[f'HI_{g}_{org}_posterior'] = hi_post
                
        result['pymc_trace'] = trace
                
    return result
