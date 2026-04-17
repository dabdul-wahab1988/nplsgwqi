from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = ROOT / "manuscript"


def parse_reviewer_comments(md_path: Path) -> dict[str, str]:
    text = md_path.read_text(encoding="utf-8")
    parts = re.split(r"^##\s+(C\d+)\s*$", text, flags=re.M)
    out: dict[str, str] = {}
    for i in range(1, len(parts), 2):
        cid = parts[i].strip()
        body = parts[i + 1].strip("\n").strip()
        out[cid] = body
    if not out:
        raise ValueError(f"No comments parsed from {md_path}")
    return out


def find_line(path: Path, needle: str) -> int:
    lines = path.read_text(encoding="utf-8").splitlines()
    for idx, line in enumerate(lines, 1):
        if needle in line:
            return idx
    raise ValueError(f"Needle not found in {path}: {needle!r}")


def esc_table_cell(s: str) -> str:
    s = s.replace("|", "\\|")
    s = s.replace("\r\n", "\n").replace("\r", "\n")
    s = s.replace("\n\n", "<br><br>").replace("\n", "<br>")
    return s


def loc(path: str, line: int, note: str) -> str:
    return f"`{path}#L{line}` ({note})"


def main() -> None:
    comments = parse_reviewer_comments(MANUSCRIPT_DIR / "reviewer_comments.md")

    anchors = {
        "cbe_methods": ("manuscript/sections/03-methodology/subsection-02-preprocessing.md", "CBE thresholding and use."),
        "cbe_results": ("manuscript/sections/04-results/subsection-02-data-integrity.md", "Data integrity was screened"),
        "sampling_qaqc": ("manuscript/sections/03-methodology/subsection-01-dataset.md", "Representativeness and metadata limitations."),
        "ilr_backinterp": ("manuscript/sections/03-methodology/subsection-04-compositional.md", "SBP/basis sensitivity."),
        "phreeqc_methods": ("manuscript/sections/03-methodology/subsection-05-geochemistry.md", "Speciation and saturation-index support (PHREEQC)."),
        "bayes_methods": ("manuscript/sections/03-methodology/subsection-06-bayesian.md", "Sampler configuration and diagnostics."),
        "bayes_reloo": ("manuscript/sections/03-methodology/subsection-06-bayesian.md", "`loo_elpd_reloo`"),
        "benchmark_methods": ("manuscript/sections/03-methodology/subsection-03-npls-gwqi.md", "Benchmarking against conventional alternatives."),
        "benchmark_results": ("manuscript/sections/04-results/subsection-03-gwqi.md", "Benchmarking against conventional indexing alternatives"),
        "latent_results": ("manuscript/sections/04-results/subsection-04-latent-processes.md", "Saturation-index evaluation"),
        "bayes_results": ("manuscript/sections/04-results/subsection-05-bayesian-drivers.md", "Bayesian diagnostic and validation outputs"),
        "risk_results": ("manuscript/sections/04-results/subsection-06-health-risk.md", "To support multi-threshold decision use"),
        "risk_methods": ("manuscript/sections/03-methodology/subsection-07-health-risk.md", "Pathways, thresholds, and uncertainty drivers."),
        "irr_results": ("manuscript/sections/04-results/subsection-07-irrigation.md", "Per-sample irrigation indices"),
        "irr_methods": ("manuscript/sections/03-methodology/subsection-08-irrigation.md", "compared against classical"),
        "region_fix": ("manuscript/sections/02-study-area/subsection-01-location.md", "Northern Region"),
        "repro_methods": ("manuscript/sections/03-methodology/subsection-09-validation.md", "analysis workflow diagram and a concise data dictionary"),
        "workflow": ("manuscript/analysis_workflow.md", "# Analysis Workflow"),
        "dictionary": ("manuscript/data_dictionary.md", "# Data Dictionary"),
        "supp_addendum": ("manuscript/supplementary-methodology.md", "# S5. Reviewer-requested additions"),
        "supp_reloo": ("manuscript/supplementary-methodology.md", "`loo_elpd_reloo`"),
        "cbe_sensitivity": ("manuscript/supplementary-methodology.md", "recomputed the Bayesian risk summary under a screening rule"),
        "code_avail": ("manuscript/sections/07-conclusion/section.md", "**Code Availability:**"),
        "mgmt_decisions": ("manuscript/sections/06-management-implications/section.md", "Benchmarking against conventional WQI baselines"),
    }

    # Resolve anchors to line numbers
    anchor_lines: dict[str, str] = {}
    for key, (rel_path, needle) in anchors.items():
        line = find_line(ROOT / rel_path, needle)
        anchor_lines[key] = loc(rel_path, line, needle)

    def L(*keys: str) -> str:
        return "; ".join(anchor_lines[k] for k in keys)

    responses: dict[str, dict[str, str]] = {
        "C01": {"status": "Addressed", "loc": "N/A (commendation; no manuscript change)", "resp": "Thank you. No change required beyond incorporating the substantive revisions listed below."},
        "C02": {"status": "Addressed", "loc": L("cbe_methods", "phreeqc_methods", "benchmark_methods", "bayes_methods", "repro_methods", "supp_addendum"), "resp": "Implemented reviewer-requested strengthening where supported: explicit CBE reporting, PHREEQC SI analysis, WQI benchmarking under weight uncertainty, expanded Bayesian diagnostics reporting, and added reproducibility metadata (workflow + data dictionary). Limitations are stated where required metadata are unavailable."},
        "C03": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C04": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C05": {"status": "Addressed", "loc": L("bayes_methods", "bayes_reloo", "supp_reloo"), "resp": "Added PSIS-LOO/WAIC reporting and PPC plots for the endpoint-driver layer and referenced them in Methods/Supplementary (Table 19; Figure S10). Because PSIS-LOO flagged an influential observation (Pareto-k > 1), we report an exact-LOO correction for the flagged point (reloo-style) as `loo_elpd_reloo` in Table 19."},
        "C06": {"status": "Partially addressed", "loc": L("bayes_methods", "bayes_reloo", "supp_reloo"), "resp": "Added PSIS-LOO/WAIC, an exact-LOO correction for PSIS failures (reloo-style; Table 19 `loo_elpd_reloo`), and prior-scale sensitivity (Table 20). External validation is not possible with the available data (single n=34 snapshot, no holdout set), so results are qualified accordingly."},
        "C07": {"status": "Addressed", "loc": L("ilr_backinterp", "latent_results", "supp_addendum"), "resp": "Expanded the ILR→CLR back-interpretation narrative and added basis/SBP sensitivity analysis (Table 21) to demonstrate stability of dominant loading patterns."},
        "C08": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C09": {"status": "Partially addressed", "loc": L("sampling_qaqc"), "resp": "Clarified the ~5-minute purge protocol and in-situ field-parameter measurement. Well volumes/stabilization logs and trip blank/duplicate metadata were not available in the provided materials and are therefore stated as limitations."},
        "C10": {"status": "Partially addressed", "loc": L("cbe_methods", "cbe_results", "cbe_sensitivity"), "resp": "Specified CBE thresholds (±10% primary; ±5% reported) and reported pass fractions and distribution (Table 16–Table 17; Figure 3). LOD/LOQ metadata and censoring flags were not available and are stated as a limitation. To demonstrate that keeping all n=34 samples does not materially drive risk conclusions, we report a CBE-screen sensitivity of the Bayesian risk summary (Table 28)."},
        "C11": {"status": "Partially addressed", "loc": L("supp_addendum"), "resp": "Reinforced that exceedance probabilities are conditional on dry-season sampling and added an explicit seasonality limitation statement. Wet-season scenario analysis is not defensible without additional seasonal data or external longitudinal datasets."},
        "C12": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C13": {"status": "Addressed", "loc": L("benchmark_methods", "benchmark_results", "supp_addendum"), "resp": "Added head-to-head benchmarking vs conventional arithmetic WQI and a fuzzy-only baseline, plus Monte Carlo weight-uncertainty analysis (Table 22–Table 23; Figure S11)."},
        "C14": {"status": "Addressed", "loc": L("benchmark_results", "mgmt_decisions"), "resp": "Expanded decision relevance by reporting where class assignments diverge under conventional vs neutrosophic indexing and under weight uncertainty, and by stating operational implications for boundary/degraded sites."},
        "C15": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C16": {"status": "Addressed", "loc": L("risk_results"), "resp": "Retained demographic disaggregation and strengthened multi-threshold exceedance reporting (Table 27) to support policy targeting."},
        "C17": {"status": "Addressed", "loc": L("code_avail", "workflow", "dictionary"), "resp": "Added explicit code availability language and provided workflow + data dictionary materials to support regulator-facing reproducibility."},
        "C18": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C19": {"status": "Addressed", "loc": L("bayes_methods", "repro_methods", "supp_addendum"), "resp": "Moved/added key implementation specifics into the Methods and Supplementary addendum (diagnostics, benchmarking, basis sensitivity, PHREEQC SI) rather than leaving them implicit."},
        "C20": {"status": "Addressed", "loc": L("workflow", "dictionary", "repro_methods"), "resp": "Added an analysis workflow diagram and a concise data dictionary to support replication from the included `data.csv`."},
        "C21": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C22": {"status": "Partially addressed", "loc": L("sampling_qaqc"), "resp": "Clarified the sampling frame as WVI-maintained public boreholes and stated that well depth/screen interval and proximity-to-source metadata are unavailable in the provided materials, limiting representativeness inference."},
        "C23": {"status": "Not possible with available data", "loc": L("phreeqc_methods"), "resp": "DO/Eh and temperature fields are not present in `data.csv`. For saturation indices, we report a temperature sensitivity band (20/25/30 °C) and disclose this assumption; redox-proxy reporting is not possible without additional measurements."},
        "C24": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C25": {"status": "Addressed", "loc": L("bayes_methods", "bayes_reloo", "supp_reloo"), "resp": "Reported Bayesian sampler configuration/diagnostics (R-hat, ESS), PPC, and PSIS-LOO/WAIC summaries (Table 19; Figure S10; Supplementary S5.2). Because PSIS-LOO flagged a Pareto-k > 1 point, we also report an exact-LOO correction (reloo-style) as `loo_elpd_reloo` in Table 19. Prior-scale sensitivity is reported in Table 20 (Supplementary S5.3)."},
        "C26": {"status": "Partially addressed", "loc": L("benchmark_results", "supp_addendum"), "resp": "Dynamic weights and classification stability are already reported; we further added benchmarking vs a conventional equal-weight WQI and a fuzzy-only baseline under weight uncertainty (Table 22–Table 23). A direct comparison to expert-assigned weights is not possible because a defensible expert-weight vector was not provided."},
        "C27": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C28": {"status": "Partially addressed", "loc": L("phreeqc_methods", "latent_results", "supp_addendum"), "resp": "Added PHREEQC saturation indices for calcite/dolomite/fluorite with temperature sensitivity (Table 18; Figure S9). Fluorapatite SI cannot be computed without phosphate data and is disclosed as a limitation."},
        "C29": {"status": "Partially addressed", "loc": L("latent_results"), "resp": "Halite-dissolution plausibility is addressed via Na/Cl diagnostics; Na-normalized diagram additions were not implemented in this revision and are noted as a potential future enhancement."},
        "C30": {"status": "Addressed", "loc": L("latent_results"), "resp": "Aligned terminology by interpreting conductivity enrichment explicitly as evapoconcentration/mineralisation where relevant."},
        "C31": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C32": {"status": "Partially addressed", "loc": L("risk_methods", "supp_addendum"), "resp": "Clarified ingestion-only pathway choice and pointed to Supplementary distributions; added an uncertainty-driver decomposition (Table 26) as a compact analogue to tornado-style uncertainty ranking under the current model structure."},
        "C33": {"status": "Addressed", "loc": L("risk_methods", "risk_results"), "resp": "Clarified that HI>1 is the primary hazard threshold and added multi-threshold exceedance probabilities (Table 27). Justified ingestion-only modelling and explicitly reported which pathways are excluded."},
        "C34": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C35": {"status": "Addressed", "loc": L("irr_methods", "irr_results", "supp_addendum"), "resp": "Ensured the seven irrigation indices are explicitly listed; added per-sample reporting and disagreement summary vs classical SAR/SSP/RSC criteria (Table 24–Table 25)."},
        "C36": {"status": "Addressed", "loc": "N/A (section header in comments)", "resp": "Noted (comment section header)."},
        "C37": {"status": "Addressed", "loc": L("cbe_methods", "cbe_results", "cbe_sensitivity"), "resp": "Provided explicit CBE thresholds, pass fractions, and distribution summary (Table 16–Table 17; Figure 3). Added a sensitivity check showing the Bayesian risk summary changes only modestly under a |CBE| ≤ 10% screening rule (Table 28), supporting the choice to treat CBE as a transparent integrity diagnostic rather than a hard exclusion criterion."},
        "C38": {"status": "Not possible with available data", "loc": L("cbe_methods"), "resp": "LOD/LOQ and censoring metadata are not present in the supplied dataset materials, so censored-value handling cannot be specified; this is reported explicitly as a limitation."},
        "C39": {"status": "Partially addressed", "loc": L("sampling_qaqc"), "resp": "Clarified the stated purge practice (~5 minutes) and noted that per-well volumes and stabilization logs are not available, limiting further justification."},
        "C40": {"status": "Addressed", "loc": L("ilr_backinterp", "supp_addendum"), "resp": "Added an explicit ILR→CLR back-interpretation explanation and basis sensitivity summary (Table 21)."},
        "C41": {"status": "Partially addressed", "loc": L("bayes_methods", "bayes_reloo", "supp_reloo"), "resp": "Added sampler diagnostics, PPC, PSIS-LOO/WAIC, and prior sensitivity outputs (Table 19–Table 20; Figure S10). Where PSIS-LOO failed (Pareto-k > 1), we computed an exact-LOO correction (reloo-style) and report it as `loo_elpd_reloo` in Table 19. Diagnostics are reported transparently, and weaker stability in the fluoride driver model is explicitly qualified in interpretation."},
        "C42": {"status": "Addressed", "loc": L("benchmark_methods", "benchmark_results", "supp_addendum"), "resp": "Implemented benchmarking against conventional and fuzzy-only WQI baselines and quantified conventional-WQI class stability under Monte Carlo weight uncertainty (Table 22–Table 23; Figure S11)."},
        "C43": {"status": "Addressed", "loc": L("risk_methods", "risk_results", "supp_addendum"), "resp": "Clarified exposure pathways (ingestion-only), provided multi-threshold exceedance probabilities (Table 27), and added an uncertainty-driver decomposition (Table 26). Parameter distributions are documented in the Supplementary methodology."},
        "C44": {"status": "Partially addressed", "loc": L("phreeqc_methods", "supp_addendum"), "resp": "Added PHREEQC saturation-index analysis for key carbonate/fluoride phases (Table 18; Figure S9) and disclosed limitations (temperature not in `data.csv`; phosphate not measured for fluorapatite SI)."},
        "C45": {"status": "Addressed", "loc": L("region_fix"), "resp": "Standardised the manuscript to state that Karaga District is in the Northern Region of Ghana and removed the inconsistent North East Region phrasing."},
        "C46": {"status": "Addressed", "loc": L("irr_methods", "irr_results"), "resp": "Explicitly listed the seven N-ISI indices and added a classical-comparison disagreement summary and per-sample reporting (Table 24–Table 25)."},
        "C47": {"status": "Partially addressed", "loc": L("sampling_qaqc"), "resp": "Removed specific CRM catalogue codes that cannot be verified from the provided materials and restated QA/QC more conservatively as use of appropriate standards and certified reference materials; detailed CRM identifiers remain unavailable."},
    }

    missing = sorted(set(comments) - set(responses))
    extra = sorted(set(responses) - set(comments))
    if missing or extra:
        raise SystemExit(f"Response mapping mismatch. Missing={missing} Extra={extra}")

    out_lines = [
        "# Response to Reviewer Comments",
        "",
        "This document responds comment-by-comment to `manuscript/reviewer_comments.md` (verbatim anchor).",
        "",
        "Statuses: `Addressed`, `Partially addressed`, `Not possible with available data`.",
        "",
    ]

    for cid in sorted(comments.keys(), key=lambda x: int(x[1:])):
        ctext = esc_table_cell(comments[cid])
        resp = esc_table_cell(responses[cid]["resp"])
        status = responses[cid]["status"]
        location = responses[cid]["loc"]

        out_lines += [
            f"## {cid}",
            "",
            "| Comment ID | Reviewer Comment (verbatim) | Response | Manuscript Change Location | Status |",
            "|---|---|---|---|---|",
            f"| `{cid}` | {ctext} | {resp} | {location} | {status} |",
            "",
        ]

    out_path = MANUSCRIPT_DIR / "response_to_comments.md"
    out_path.write_text("\n".join(out_lines).rstrip() + "\n", encoding="utf-8")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
