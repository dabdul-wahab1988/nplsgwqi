import numpy as np
import pandas as pd

class Preprocessor:
    # Stoichiometric conversion factors (Valence / Atomic Weight)
    # and Equivalent Weights (Atomic Weight / Valence)
    # Factor = 1 / Eq_Weight
    STOI_VALENCE_WEIGHTS = {
        'Ca': {'factor': 2 / 40.078, 'eq_weight': 20.039},
        'Mg': {'factor': 2 / 24.305, 'eq_weight': 12.1525},
        'Na': {'factor': 1 / 22.990, 'eq_weight': 22.990},
        'K': {'factor': 1 / 39.098, 'eq_weight': 39.098},
        'Cl': {'factor': 1 / 35.45, 'eq_weight': 35.45},
        'HCO3': {'factor': 1 / 61.016, 'eq_weight': 61.016},
        'SO4': {'factor': 2 / 96.06, 'eq_weight': 48.03},
        'NO3': {'factor': 1 / 62.004, 'eq_weight': 62.004},
        'CO3': {'factor': 2 / 60.01, 'eq_weight': 30.005},
        'F': {'factor': 1 / 18.998, 'eq_weight': 18.998}
    }

    def __init__(self, main_compositional_vars: list, trace_vars: list = None, non_compositional_vars: list = None, endpoint_vars: list = None):
        self.main_compositional_vars = main_compositional_vars
        self.trace_vars = trace_vars or []
        self.non_compositional_vars = non_compositional_vars or []
        self.endpoint_vars = endpoint_vars or []
        
    def to_meq(self, data: pd.DataFrame) -> pd.DataFrame:
        """Converts mg/L to meq/L for relevant ions."""
        df_meq = data.copy()
        for ion, props in self.STOI_VALENCE_WEIGHTS.items():
            if ion in df_meq.columns:
                df_meq[ion] *= props['factor']
        return df_meq

    def calculate_charge_balance(self, data_meq: pd.DataFrame) -> pd.Series:
        """Calculates Charge Balance Error (CBE) in %."""
        cations = ['Ca', 'Mg', 'Na', 'K']
        anions = ['Cl', 'HCO3', 'SO4', 'NO3', 'CO3', 'F']
        
        sum_cat = data_meq[[c for c in cations if c in data_meq.columns]].sum(axis=1)
        sum_ani = data_meq[[a for a in anions if a in data_meq.columns]].sum(axis=1)
        
        cbe = ((sum_cat - sum_ani) / (sum_cat + sum_ani)) * 100
        return cbe

    def partition_variables(self, data: pd.DataFrame) -> dict:
        """
        Classifies variables into blocks.
        """
        partitions = {
            'main_comp': data[[c for c in self.main_compositional_vars if c in data.columns]],
            'trace_comp': data[[c for c in self.trace_vars if c in data.columns]],
            'non_comp': data[[c for c in self.non_compositional_vars if c in data.columns]],
            'endpoints': data[[c for c in self.endpoint_vars if c in data.columns]]
        }
        return partitions

    def replace_zeros(self, data: pd.DataFrame, method='lognormal', lod_dict: dict = None) -> pd.DataFrame:
        """
        Improved zero/nondetect handling.
        """
        df = data.copy()
        for col in df.columns:
            if not pd.api.types.is_numeric_dtype(df[col]):
                continue
            mask = df[col] <= 0
            if mask.any():
                if method == 'lognormal':
                    lod = lod_dict.get(col, df[col][df[col] > 0].min()) if lod_dict else df[col][df[col] > 0].min()
                    df.loc[mask, col] = lod / np.sqrt(2)
                else:
                    delta = 0.01
                    df.loc[mask, col] = delta
        return df

    def handle_censor_nondetect(self, data: pd.DataFrame, lod_dict: dict) -> pd.DataFrame:
        df = data.copy()
        for col, lod in lod_dict.items():
            if col in df.columns:
                df.loc[df[col] < lod, col] = lod / np.sqrt(2)
        return df

    def get_augmented_predictor_block(self, data: pd.DataFrame, endpoint_name: str) -> pd.DataFrame:
        """
        Generates a high-signal predictor block by combining:
        1. Major ions (Log-transformed)
        2. Field parameters (pH, EC)
        3. Geochemical ratios and indices (Na/Cl, CAI, Weathering proxies)
        """
        # Convert to meq for index calculations
        data_meq = self.to_meq(data)
        
        # Ensure no zeros in raw and meq data for ratios
        df_raw = self.replace_zeros(data.select_dtypes(include=[np.number]))
        df_meq = self.replace_zeros(data_meq.select_dtypes(include=[np.number]))
        
        # 1. Compositional Block (Log-transformed raw concentrations)
        ions = ['Ca', 'Mg', 'Na', 'K', 'Cl', 'SO4', 'HCO3', 'NO3', 'F']
        if endpoint_name in ions:
            ions.remove(endpoint_name)
        X_ions = np.log(df_raw[ions])
            
        # 2. Field Parameters
        field = [c for c in ['pH', 'EC', 'TDS'] if c in df_raw.columns]
        X_field = df_raw[field]
        
        # 3. Geochemical Ratios and Indices (Calculated in meq space)
        # Ion Exchange Proxies
        df_meq['CAI_1'] = (df_meq['Cl'] - (df_meq['Na'] + df_meq['K'])) / df_meq['Cl']
        
        cai2_denom_cols = [c for c in ['SO4', 'HCO3', 'CO3', 'NO3'] if c in df_meq.columns]
        denom_cai2 = df_meq[cai2_denom_cols].sum(axis=1).replace(0, 1e-6)
        df_meq['CAI_2'] = (df_meq['Cl'] - (df_meq['Na'] + df_meq['K'])) / denom_cai2
        
        # Weathering and Salinity Proxies
        df_meq['Na_Cl'] = df_meq['Na'] / df_meq['Cl']
        df_meq['Ca_HCO3'] = df_meq['Ca'] / df_meq['HCO3']
        df_meq['Mg_Ca'] = df_meq['Mg'] / df_meq['Ca']
        
        # Silicate vs Carbonate Weathering
        tz_plus_cols = [c for c in ['Ca', 'Mg', 'Na', 'K'] if c in df_meq.columns]
        tz_plus = df_meq[tz_plus_cols].sum(axis=1).replace(0, 1e-6)
        df_meq['Silicate_Proxy'] = (df_meq['Na'] + df_meq['K']) / tz_plus
        df_meq['Carbonate_Proxy'] = (df_meq['Ca'] + df_meq['Mg']) / df_meq['HCO3']

        indices = ['CAI_1', 'CAI_2', 'Na_Cl', 'Ca_HCO3', 'Mg_Ca', 'Silicate_Proxy', 'Carbonate_Proxy']
        
        # Handle potential negatives/zeros in indices before logging where appropriate
        # For CAI, we keep them linear as they are often negative (ion exchange direction)
        # For ratios, we use log to linearize
        X_indices_linear = df_meq[['CAI_1', 'CAI_2', 'Silicate_Proxy', 'Carbonate_Proxy']]
        X_ratios_log = np.log(df_meq[['Na_Cl', 'Ca_HCO3', 'Mg_Ca']])
        
        # Combine
        X = pd.concat([X_ions, X_field, X_indices_linear, X_ratios_log], axis=1)
        
        # Clean up column names for the model
        X.columns = [f"Aug_{c}" for c in X.columns]
        return X

