import numpy as np
import pandas as pd

def closure(data: pd.DataFrame, kappa: float = 1.0) -> pd.DataFrame:
    """Applies closure operation (rescaling to sum to kappa)."""
    return data.div(data.sum(axis=1), axis=0) * kappa

def clr_transform(data: pd.DataFrame) -> pd.DataFrame:
    """Centered Log-Ratio transform."""
    gm = np.exp(np.mean(np.log(data), axis=1))
    return np.log(data.div(gm, axis=0))

def ilr_transform(data: pd.DataFrame, V: np.ndarray = None) -> pd.DataFrame:
    """Isometric Log-Ratio transform using a given basis matrix V."""
    D = data.shape[1]
    if V is None:
        # Default Helmert contrast matrix if V is not provided
        V = np.zeros((D, D - 1))
        for i in range(D - 1):
            v = np.zeros(D)
            v[:i+1] = 1.0 / np.sqrt((i+1)*(i+2))
            v[i+1] = -np.sqrt((i+1) / (i+2))
            V[:, i] = v
    
    clr_data = clr_transform(data).values
    ilr_data = clr_data @ V
    return pd.DataFrame(ilr_data, index=data.index, columns=[f'ILR_{i+1}' for i in range(D-1)])

def split_score_regimes(scores: pd.DataFrame) -> pd.DataFrame:
    """
    Splits scores into positive and negative regimes.
    """
    pos_scores = np.maximum(scores, 0)
    neg_scores = np.maximum(-scores, 0)
    
    pos_scores.columns = [f"{col}(+)" for col in scores.columns]
    neg_scores.columns = [f"{col}(-)" for col in scores.columns]
    
    return pd.concat([pos_scores, neg_scores], axis=1)
