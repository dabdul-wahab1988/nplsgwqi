import os
import pandas as pd
from .preprocessing import Preprocessor
from .process_discovery import run_process_discovery
from .process_regression import run_process_contribution_model
from .wqi import calculate_npls_gwqi
from .risk import assess_health_risk
from .irrigation import calculate_irrigation_indices, assess_irrigation_suitability
from .validation import bootstrap_process_stability, cross_validate_regression, check_endpoint_leakage, report_model_diagnostics
from .geochemistry import calculate_geochemical_constraints, audit_process_names

class IntegratedWorkflow:
    def __init__(self, main_comp_vars, endpoint_vars, standards, toxref=None, use_ridge: bool = True):
        self.prep = Preprocessor(main_compositional_vars=main_comp_vars, endpoint_vars=endpoint_vars)
        self.standards = standards
        self.toxref = toxref
        self.use_ridge = use_ridge

    def run(self, data: pd.DataFrame, n_components: int = None, robust: bool = True):
        # 1. Preprocessing (Stoichiometry and Charge Balance)
        data_meq = self.prep.to_meq(data)
        cbe = self.prep.calculate_charge_balance(data_meq)
        
        # 2. Geochemical Constraints (NEW: Deterministic Audit)
        constraints = calculate_geochemical_constraints(data_meq)

        partitions = self.prep.partition_variables(data)
        comp_data = self.prep.replace_zeros(partitions['main_comp'])
        endpoints = partitions['endpoints']

        # 3. Process Discovery
        discovery_res = run_process_discovery(comp_data, n_components=n_components, robust=robust)
        
        # 4. Refine Process Names using Geochemical Audit
        discovery_res['suggested_names'] = audit_process_names(discovery_res['loadings_clr'], constraints)
        n_comp_actual = discovery_res['scores'].shape[1]

        # 5. Process Regression with Leakage Detection
        regression_results = {}
        for ep in endpoints.columns:
            if check_endpoint_leakage(self.prep.main_compositional_vars, ep):
                sub_comp = comp_data.drop(columns=[ep])
                sub_discovery = run_process_discovery(sub_comp, n_components=n_comp_actual, robust=robust)
                sub_discovery['suggested_names'] = audit_process_names(sub_discovery['loadings_clr'], constraints)
                regression_results[ep] = run_process_contribution_model(sub_discovery['scores'], endpoints[ep], use_ridge=self.use_ridge)
                regression_results[ep]['basis_name'] = f"{ep}_safe"
                regression_results[ep]['basis_label'] = f"{ep}-safe hydrochemical basis"
                regression_results[ep]['component_name_map'] = sub_discovery['suggested_names']
            else:
                regression_results[ep] = run_process_contribution_model(discovery_res['scores'], endpoints[ep], use_ridge=self.use_ridge)
                regression_results[ep]['basis_name'] = 'full'
                regression_results[ep]['basis_label'] = 'Full hydrochemical basis'
                regression_results[ep]['component_name_map'] = discovery_res['suggested_names']

        # 6. WQI (Raw Space)
        wqi_res, weights = calculate_npls_gwqi(data[self.standards.index].copy(), self.standards)
        regression_results['NPLS_GWQI'] = run_process_contribution_model(discovery_res['scores'], wqi_res['NPLS_GWQI'], use_ridge=self.use_ridge)
        regression_results['NPLS_GWQI']['basis_name'] = 'full'
        regression_results['NPLS_GWQI']['basis_label'] = 'Full hydrochemical basis'
        regression_results['NPLS_GWQI']['component_name_map'] = discovery_res['suggested_names']

        # 7. Risk (Raw Space)
        risk_res = assess_health_risk(data, toxref=self.toxref, mode='deterministic')

        # 8. Irrigation Suitability (NEW)
        irr_indices = calculate_irrigation_indices(data_meq)
        irr_suitability_res = assess_irrigation_suitability(irr_indices)
        irr_results = pd.concat([irr_indices, irr_suitability_res], axis=1)

        # 9. Validation and Diagnostics
        stability = bootstrap_process_stability(comp_data, n_components=n_comp_actual, n_bootstraps=20)
        diagnostics = report_model_diagnostics(regression_results)

        return {
            'charge_balance_error': cbe,
            'geochemical_constraints': constraints,
            'discovery': discovery_res,
            'regression': regression_results,
            'wqi': wqi_res,
            'wqi_weights': weights,
            'risk': risk_res,
            'irrigation': irr_results,
            'validation': stability,
            'diagnostics': diagnostics
        }

def generate_artifacts(output_dir: str = 'manuscript/artifacts') -> None:
    os.makedirs(output_dir, exist_ok=True)
    # Figure 1 is provided externally (QGIS)
