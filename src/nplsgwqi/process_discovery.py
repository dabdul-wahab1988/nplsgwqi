import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.covariance import MinCovDet
from .compositional import ilr_transform


def helmert_basis(n_parts: int) -> np.ndarray:
    V = np.zeros((n_parts, n_parts - 1))
    for i in range(n_parts - 1):
        v = np.zeros(n_parts)
        v[: i + 1] = 1.0 / np.sqrt((i + 1) * (i + 2))
        v[i + 1] = -np.sqrt((i + 1) / (i + 2))
        V[:, i] = v
    return V

def parallel_analysis(data: pd.DataFrame, n_iterations: int = 100) -> int:
    """Selects the number of components using Parallel Analysis."""
    n_samples, n_features = data.shape
    obs_pca = PCA().fit(data)
    obs_eig = obs_pca.explained_variance_
    
    rnd_eig = []
    for _ in range(n_iterations):
        rnd_data = np.random.normal(size=(n_samples, n_features))
        rnd_pca = PCA().fit(rnd_data)
        rnd_eig.append(rnd_pca.explained_variance_)
    
    # Compare observed eigenvalues to 95th percentile of random ones
    crit_eig = np.percentile(np.array(rnd_eig), 95, axis=0)
    # Ensure at least 1 component
    count = int(np.sum(obs_eig > crit_eig))
    return max(1, count)

class ProcessDiscoveryEngine:
    def __init__(self, n_components: int = None, method: str = 'unsupervised', robust: bool = True):
        self.n_components = n_components
        self.method = method
        self.robust = robust
        self.pca_model = None
        self.loadings_clr = None
        self.ilr_basis_ = None
        self.used_robust_projection_ = False
        self.columns_ = None

    def fit(self, comp_data: pd.DataFrame, endpoints: pd.DataFrame = None):
        if self.method == 'unsupervised':
            self.columns_ = list(comp_data.columns)
            self.ilr_basis_ = helmert_basis(comp_data.shape[1])
            ilr_data = ilr_transform(comp_data, V=self.ilr_basis_)

            if self.n_components is None:
                self.n_components = parallel_analysis(ilr_data)

            if self.robust:
                # Robust PCA using Minimum Covariance Determinant (MCD)
                # Note: MinCovDet requires n_samples > n_features
                if ilr_data.shape[0] > ilr_data.shape[1] + 1:
                    mcd = MinCovDet().fit(ilr_data)
                    self.pca_model = PCA(n_components=self.n_components)
                    self.pca_model.fit(mcd.covariance_)
                    self.used_robust_projection_ = True
                else:
                    # Fallback to standard PCA if samples are too few
                    self.pca_model = PCA(n_components=self.n_components)
                    self.pca_model.fit(ilr_data)
                    self.used_robust_projection_ = False
            else:
                self.pca_model = PCA(n_components=self.n_components)
                self.pca_model.fit(ilr_data)
                self.used_robust_projection_ = False

            # Convert components back to clr space for interpretability
            # The PCA components are in ILR space. To get them in CLR space:
            # Loadings_CLR = PCA_Components @ V.T
            loadings_clr_val = self.pca_model.components_ @ self.ilr_basis_.T

            self.loadings_clr = pd.DataFrame(
                loadings_clr_val,
                columns=comp_data.columns,
                index=[f'Process_{i+1}' for i in range(self.n_components)]
            )
            return self
        else:
            # Endpoint-aware or other methods
            return self.fit(comp_data)

    def transform(self, comp_data: pd.DataFrame) -> pd.DataFrame:
        if self.pca_model is None or self.ilr_basis_ is None or self.columns_ is None:
            raise ValueError("ProcessDiscoveryEngine must be fit before calling transform().")
        comp_aligned = comp_data.loc[:, self.columns_]
        ilr_data = ilr_transform(comp_aligned, V=self.ilr_basis_)
        if self.used_robust_projection_:
            scores_ilr = ilr_data.values @ self.pca_model.components_.T
        else:
            scores_ilr = self.pca_model.transform(ilr_data)
        return pd.DataFrame(
            scores_ilr,
            index=comp_data.index,
            columns=[f'Process_{i+1}' for i in range(self.n_components)]
        )

    def fit_transform(self, comp_data: pd.DataFrame, endpoints: pd.DataFrame = None) -> pd.DataFrame:
        return self.fit(comp_data, endpoints=endpoints).transform(comp_data)
            
    def suggest_process_names(self):
        """
        Rule-based process naming based on dominant CLR loadings.
        """
        names = {}
        if self.loadings_clr is None:
            return names
            
        for proc in self.loadings_clr.index:
            loadings = self.loadings_clr.loc[proc]
            top_positive = loadings[loadings > 0].nlargest(2).index.tolist()
            top_negative = loadings[loadings < 0].nsmallest(2).index.tolist()
            
            if 'Na' in top_positive and 'Cl' in top_positive:
                names[proc] = 'Salinization / Evaporation'
            elif 'Ca' in top_positive and 'HCO3' in top_positive:
                names[proc] = 'Carbonate Weathering'
            elif 'NO3' in top_positive or 'Cl' in top_positive:
                names[proc] = 'Anthropogenic Loading'
            elif 'F' in top_positive:
                names[proc] = 'Fluoride Mobilization'
            else:
                names[proc] = f'Mixed Process ({", ".join(top_positive)} vs {", ".join(top_negative)})'
                
        return names

def run_process_discovery(comp_data: pd.DataFrame, n_components: int = None, robust: bool = True) -> dict:
    engine = ProcessDiscoveryEngine(n_components=n_components, robust=robust)
    scores = engine.fit_transform(comp_data)
    names = engine.suggest_process_names()
    
    return {
        'scores': scores,
        'loadings_clr': engine.loadings_clr,
        'suggested_names': names,
        'explained_variance': engine.pca_model.explained_variance_ratio_
    }
