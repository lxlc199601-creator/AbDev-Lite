# AbDev-Lite v1.0

AbDev-Lite v1.0 is a local Streamlit platform for integrated antibody variable-region developability screening. It supports internal demos, batch candidate triage, result summary, leadership reporting, and future extension while keeping optional external evidence fully user controlled.

AbDev-Lite v1.0 provides rule-based and imported computational screening outputs for antibody variable-region developability assessment. Results are intended for candidate triage, reporting, and planning support only. They do not replace experimental binding, expression, stability, immunogenicity, structural, formulation, or CMC studies.

## v1.0 Feature Overview

- Sequence input cleaning and QC
- Physicochemical property calculation
- Rule-based liability scanning
- Optional IMGT numbering and CDR/FR mapping
- Liability-to-region mapping
- CDR-adjusted risk scoring
- Optional humanness/germline result import
- Candidate prioritization with A/B/C/D classes
- Formulation feature extraction and Expreso-lite adapter with rule-based fallback
- Optional structure result import and structural risk summary
- External tool adapter framework with local input package generation and manual result import
- Final integrated assessment with go/no-go suggestion and evidence completeness level
- Executive decision summary for Streamlit and HTML reporting
- Excel and HTML report export
- Streamlit v1.0 dashboard

## System Flow

Input antibody variable-region sequences
-> Sequence QC
-> Physicochemical properties
-> Liability scanning
-> IMGT CDR/FR mapping
-> Humanness import
-> Candidate ranking
-> Formulation recommendation
-> Structure result import
-> External tool result import
-> Final integrated assessment
-> Excel / HTML report

## Installation

```bash
pip install -r requirements.txt
```

AbNumber/ANARCI supports IMGT numbering when available. If it is not installed or cannot number a sequence, the rest of the pipeline still runs and records skipped or warning status rows.

## Run

```bash
streamlit run app.py
```

Upload `data/example_input.xlsx` or a compatible CSV/XLSX input, optionally upload supported result files, then click **Run Analysis**.

## Required Input Format

| Column | Meaning |
| --- | --- |
| antibody_id | Antibody identifier, such as Ab001 |
| molecule_format | mAb, BsAb, VHH, scFv, Fab, or a user label |
| chain_id | Chain identifier, such as VH, VL, VH1, VL1 |
| chain_type | heavy, light, vhh, scfv, or unknown |
| region_type | variable_region or full_length |
| sequence_scope | VH, VL, VHH, scFv, full_heavy, full_light, or related scope |
| sequence | Amino acid sequence |

## Optional Input Files

- `data/example_humanness_results.xlsx`: imported humanness/germline evidence
- `data/example_structure_results.xlsx`: imported structure summary evidence
- `data/example_external_tool_results.xlsx`: imported external tool result evidence

If an optional file is not provided, the app displays `Not provided / optional module skipped` and continues.

## Output Files

Reports are generated under `outputs/`:

- `outputs/abdev_lite_results.xlsx`
- `outputs/abdev_lite_report.html`

The Excel workbook includes `Final_Assessment` plus all upstream result sheets: cleaned input, QC, properties, numbering, region summary, liability tables, chain scores, humanness, antibody summary, candidate ranking, formulation, structure, tool registry, external run plan, external results, and external summary.

The HTML report includes Executive Summary, Final Integrated Assessment, Candidate Prioritization, Key Risk Overview, module sections, Method Summary, and Limitations and Disclaimer.

## Module Notes

- Final assessment uses transparent rules, not machine learning.
- `go_no_go_suggestion` maps A/B/C/D classes to Advance, Advance with review, Engineering review, or Deprioritize / redesign, with escalation for high structure, formulation, external, or combined humanness/CDR risk.
- `confidence_level` means evidence completeness level only. It is not experimental confidence or probability of success.
- Browser automation remains disabled.
- No sequence is uploaded automatically to third-party tools.

## Current Limitations

- No humanization mutation design.
- No real formulation conclusion.
- No experimental conclusion or CMC decision.
- No forced AlphaFold, IgFold, BioPhi, IgBLAST, or browser automation dependency.
- Imported external results require user review.
- BsAb analysis remains chain-level and does not evaluate full molecular architecture.

## Roadmap

- Broader validation datasets and regression tests
- Additional optional adapters for permitted local/API workflows
- Improved visualization for batch comparisons
- More configurable scoring thresholds
- Expanded documentation for internal deployment patterns

## GitHub Usage

Recommended release workflow:

1. Keep runtime outputs out of Git.
2. Commit source, templates, docs, example files, and `.gitkeep` placeholders.
3. Tag the stable release as `v1.0.0`.
4. Attach release notes summarizing Final Integrated Assessment, Executive Decision Summary, report polish, and documentation stabilization.

## Data Safety

Do not commit confidential antibody sequences, external tool outputs, credentials, cookies, API keys, or environment files. Runtime folders such as `outputs/`, `external_inputs/`, `external_results/`, and `external_runs/` are ignored except for `.gitkeep` placeholders.
