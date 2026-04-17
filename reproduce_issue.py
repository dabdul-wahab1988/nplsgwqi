import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from nplsgwqi.bayesian_regression import run_bayesian_endpoint_model

def test_stability():
    data = pd.read_csv('data.csv')
    comp = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F"]
    
    for ep in ["NO3", "F"]:
        print(f"\nEvaluating Stability for {ep}...")
        y = data[ep]
        X = data[comp]
        
        # 5-fold CV
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        cv_scores = []
        for train_idx, test_idx in kf.split(data):
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
            
            # Use a simpler model for fast CV check (Ridge)
            # Or run the actual Bayesian model if we have time
            # Let's run the actual one to see why it fails
            try:
                res = run_bayesian_endpoint_model(
                    X_train, y_train, ep,
                    draws=500, tune=200, chains=1 # Fast for testing
                )
                # Prediction on test set
                # We need a predict function. bayesian_regression.py doesn't have one that works with new data easily.
                # Let's use the mean posterior coefficients
                from sklearn.preprocessing import StandardScaler
                from nplsgwqi.bayesian_regression import prepare_clr_predictors
                
                x_train_clr = prepare_clr_predictors(X_train, ep)
                x_test_clr = prepare_clr_predictors(X_test, ep)
                
                scaler = StandardScaler().fit(x_train_clr)
                X_test_scaled = scaler.transform(x_test_clr)
                
                beta = res['beta_mean']
                # alpha is not returned in the dict, but we can infer it or modify the code
                # For now, let's just see in-sample R2 first
                print(f"  Fold In-Sample R2: {res['in_sample_r2']:.4f}")
                
            except Exception as e:
                print(f"  Error in fold: {e}")

if __name__ == "__main__":
    test_stability()
