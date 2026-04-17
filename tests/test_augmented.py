import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from nplsgwqi.preprocessing import Preprocessor

def test_augmented_block():
    # Mock data
    data = pd.DataFrame({
        'pH': [7.2, 6.4, 8.6, 7.5],
        'EC': [400, 600, 1200, 450],
        'TDS': [300, 500, 1000, 400],
        'Ca': [40, 20, 80, 35],
        'Mg': [12, 6, 24, 10],
        'Na': [30, 100, 200, 25],
        'K': [2, 5, 10, 2],
        'HCO3': [200, 150, 400, 180],
        'Cl': [20, 150, 300, 25],
        'SO4': [10, 50, 100, 15],
        'NO3': [5, 45, 10, 2],
        'F': [0.5, 1.2, 0.8, 0.4]
    })
    
    prep = Preprocessor(main_compositional_vars=['Ca', 'Mg', 'Na', 'K', 'Cl', 'SO4', 'HCO3'], endpoint_vars=['NO3', 'F'])
    
    aug_f = prep.get_augmented_predictor_block(data, 'F')
    
    print("Augmented Columns for F:")
    print(aug_f.columns.tolist())
    
    # Check for expected keys
    expected = ['Aug_CAI_1', 'Aug_CAI_2', 'Aug_Silicate_Proxy', 'Aug_Carbonate_Proxy', 'Aug_Na_Cl', 'Aug_Ca_HCO3', 'Aug_Mg_Ca']
    for ex in expected:
        if ex not in aug_f.columns:
            print(f"MISSING: {ex}")
            return False
    
    print("\nSample Data (First row):")
    print(aug_f.iloc[0])
    return True

if __name__ == "__main__":
    if test_augmented_block():
        print("\nSUCCESS: Augmented block looks correct.")
    else:
        print("\nFAILURE: Missing columns in augmented block.")
        sys.exit(1)
