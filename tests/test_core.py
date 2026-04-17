import pandas as pd
import numpy as np
from nplsgwqi.core import calculate_npls_gwqi, run_process_contribution_model, assess_health_risk, generate_artifacts, calculate_irrigation_indices
from nplsgwqi.process_discovery import run_process_discovery
from nplsgwqi.preprocessing import Preprocessor
from nplsgwqi.irrigation import assess_irrigation_suitability

def get_dummy_data():
    np.random.seed(42)
    # 50 samples, 4 parameters
    data = pd.DataFrame(np.random.lognormal(mean=1.5, sigma=0.5, size=(50, 4)), 
                        columns=['pH', 'Ca', 'Mg', 'Na'])
    standards = pd.Series({'pH': 8.5, 'Ca': 75.0, 'Mg': 30.0, 'Na': 200.0})
    return data, standards

def test_irrigation():
    data, _ = get_dummy_data()
    # Add dummy HCO3, Cl, SO4 for irrigation indices
    data['HCO3'] = np.random.lognormal(2, 0.5, 50)
    data['Cl'] = np.random.lognormal(2, 0.5, 50)
    data['SO4'] = np.random.lognormal(2, 0.5, 50)
    data['K'] = np.random.lognormal(0.5, 0.2, 50)
    
    # Calculate meq/L for irrigation
    prep = Preprocessor(main_compositional_vars=['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4'])
    data_meq = prep.to_meq(data)
    
    indices = calculate_irrigation_indices(data_meq)
    assert 'SAR' in indices.columns
    assert 'RSC' in indices.columns
    assert 'PS' in indices.columns
    
    suitability = assess_irrigation_suitability(indices)
    assert 'N_ISI_Score' in suitability.columns
    assert len(suitability) == 50

def test_normalization():
    data, standards = get_dummy_data()
    assert len(data) == 50
    assert len(standards) == 4

def test_pls():
    data, standards = get_dummy_data()
    gwqi_results, weights = calculate_npls_gwqi(data, standards)
    assert 'NPLS_GWQI' in gwqi_results.columns
    assert len(weights) == 4
    # Weights should sum to 1
    assert np.isclose(weights.sum(), 1.0)

def test_process_contribution_model():
    data, _ = get_dummy_data()
    # Create positive data for CoDA
    comp_data = data[['Ca', 'Mg', 'Na']].copy()
    
    discovery = run_process_discovery(comp_data, n_components=2)
    scores = discovery['scores']
    
    target = data['pH']
    res = run_process_contribution_model(scores, target)
    
    assert 'coefficients' in res
    assert 'relative_process_contributions' in res
    # 2 processes * 2 regimes = 4 poles
    assert len(res['coefficients']) == 4
    assert len(res['relative_process_contributions']) == 5 # 4 poles + 1 unexplained

def test_bayesian_risk():
    np.random.seed(42)
    df = pd.DataFrame({
        'F': np.random.normal(1.5, 0.2, 50),
        'Pb': np.random.normal(0.01, 0.005, 50),
        'As': np.random.normal(0.005, 0.001, 50)
    })
    # Fast test with deterministic mode only
    res = assess_health_risk(df, n_samples=100, mode='deterministic')
    assert 'HQ_sample_F' in res
    assert len(res['HQ_sample_F']) == 50

def test_pipeline():
    from src.nplsgwqi.orchestrator import generate_artifacts
    generate_artifacts('manuscript/artifacts')
    # Artifacts generated successfully


