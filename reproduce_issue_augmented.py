import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from nplsgwqi.bayesian_regression import run_bayesian_endpoint_model

def test_stability_augmented():
    data = pd.read_csv('data.csv')
    comp = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4"] # Ions
    # We will add pH and EC
    
    for ep in ["NO3", "F"]:
        print(f"\nEvaluating Stability for {ep} (Augmented)...")
        y = data[ep]
        
        # Build augmented block
        X = data[comp + ["pH", "EC"]]
        
        # 5-fold CV
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        for train_idx, test_idx in kf.split(data):
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
            
            try:
                # We need to handle pH/EC correctly. 
                # For now, let's just pass them to run_bayesian_endpoint_model
                # It will treat them as compositional, which is a bit hacky but let's see if signal appears.
                res = run_bayesian_endpoint_model(
                    X_train, y_train, ep,
                    draws=500, tune=200, chains=1
                )
                print(f"  Fold In-Sample R2: {res['in_sample_r2']:.4f}")
                if res['selected_features']:
                    print(f"  Selected: {res['selected_features']}")
                
            except Exception as e:
                print(f"  Error in fold: {e}")

if __name__ == "__main__":
    test_stability_augmented()
