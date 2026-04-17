import pandas as pd
import numpy as np
from nplsgwqi.orchestrator import IntegratedWorkflow

# Mock data (increased samples for stability)
data = pd.DataFrame({
    'pH': [7.2, 6.4, 8.6, 7.5, 7.8, 6.9, 7.1, 8.2],
    'TDS': [400, 600, 1200, 450, 520, 480, 510, 800],
    'Ca': [40, 20, 80, 35, 42, 38, 41, 60],
    'Mg': [12, 6, 24, 10, 14, 11, 13, 18],
    'Na': [30, 100, 200, 25, 35, 28, 31, 150],
    'K': [2, 5, 10, 2, 3, 2, 2.5, 7],
    'HCO3': [200, 150, 400, 180, 220, 190, 210, 300],
    'Cl': [20, 150, 300, 25, 40, 30, 35, 200],
    'SO4': [10, 50, 100, 15, 20, 12, 18, 60],
    'NO3': [5, 45, 10, 2, 8, 4, 6, 20],
    'Fe': [0.1, 0.4, 0.05, 0.05, 0.1, 0.05, 0.15, 0.25],
    'As': [0.005, 0.02, 0.001, 0.002, 0.004, 0.003, 0.002, 0.015],
    'Pb': [0.001, 0.005, 0.01, 0.002, 0.001, 0.003, 0.002, 0.008]
})

standards = pd.Series({
    'pH': 8.5, 'TDS': 500, 'Ca': 75, 'Mg': 30, 'Na': 200, 'K': 12,
    'HCO3': 300, 'Cl': 250, 'SO4': 250, 'NO3': 45, 'Fe': 0.3,
    'As': 0.01, 'Pb': 0.01
})

# Setup workflow
main_comp = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3']
endpoints = ['TDS', 'pH', 'As', 'Pb']

workflow = IntegratedWorkflow(main_comp_vars=main_comp, endpoint_vars=endpoints, standards=standards)

# Run
results = workflow.run(data, n_components=2)

print("Workflow successfully executed.")
print("\n--- Irrigation Results (Advanced N-ISI) ---")
print(results['irrigation'][['SAR', 'SSP', 'MH', 'RSC', 'PS', 'N_ISI_Score', 'ISI_Class']])

print("\n--- WQI Results (with Sensory Penalty) ---")
print(results['wqi'][['NPLS_GWQI', 'Quality_Class']])

print("\n--- Health Risk (Deterministic with As Speciation) ---")
print(f"Hazard Index (Neurotoxicity) for Sample 1: {results['risk'].get('HI_sample_neuro', 'N/A')}")
