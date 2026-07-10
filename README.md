# AbDev-Lite

AbDev-Lite is a local lightweight Streamlit MVP for antibody variable-region developability screening. Version `MVP v0.8` analyzes sequence-level signals for VH, VL, VHH, and scFv-derived variable domains, supports optional IMGT numbering and CDR/FR liability mapping, imports optional external humanness/germline and structure prediction summary results, ranks candidates with a rule-based decision matrix, adds formulation recommendations, then exports Excel and offline HTML reports.

The project is intentionally local and lightweight. It does not call external webpages, upload sequences to third-party platforms, or run BioPhi, IgBLAST, Sapiens, OASis, TAP, IgFold, ImmuneBuilder, AlphaFold, ColabFold, docking, or molecular dynamics.

## Current Features

- Variable-region CSV/XLSX input
- Sequence cleaning and basic QC
- Physicochemical properties
- Rule-based liability scanning
- Optional IMGT numbering through AbNumber/ANARCI when available
- FR/CDR region extraction and liability-to-region mapping
- CDR-adjusted risk scoring
- Optional external humanness/germline result import
- Candidate prioritization with final priority score and A/B/C/D class
- Formulation feature extraction and recommendation
- Optional external structure prediction result import
- Structure confidence metrics
- CDR loop confidence summary
- Surface/aggregation/charge patch risk summary
- Structural risk class
- Structural risk integration into candidate ranking
- Excel report, HTML report, and Streamlit interface

## Input Format

Upload a CSV or XLSX file with these required columns:

| Column | Meaning |
| --- | --- |
| antibody_id | Antibody identifier, such as Ab001 or BsAb001 |
| molecule_format | mAb, BsAb, VHH, scFv, Fab, or another user label |
| chain_id | Chain identifier, such as VH, VL, VH1, VL1 |
| chain_type | heavy, light, vhh, scfv, or unknown |
| region_type | variable_region or full_length |
| sequence_scope | VH, VL, VHH, scFv, VH_arm1, VL_arm1, VH_arm2, VL_arm2, full_heavy, full_light |
| sequence | Amino acid sequence |

## Optional Humanness/Germline Import

AbDev-Lite can import an optional external humanness/germline assessment file in CSV or XLSX format. If it is not provided, sequence-level analysis, prioritization, formulation, and reporting still run.

Recommended columns include `antibody_id`, `chain_id`, `sequence_scope`, `humanness_tool`, `humanness_score`, `humanness_percentile`, `closest_human_germline`, `closest_species`, `v_gene`, `j_gene`, `identity_to_human_germline`, `framework_identity`, `cdr_identity`, and `humanness_notes`.

Human-likeness assessment is not equivalent to clinical immunogenicity prediction. AbDev-Lite does not run BioPhi, IgBLAST, Sapiens, OASis, or automatic humanization design.

## Optional Structure Result Import

AbDev-Lite v0.8 can import an optional external structure prediction summary file in CSV or XLSX format. If it is not provided, AbDev-Lite reports structure status as `Not Available` and does not penalize candidate ranking.

Recommended columns include `antibody_id`, `chain_id`, `sequence_scope`, `structure_tool`, `structure_model_file`, `structure_status`, `mean_plddt`, `cdr1_plddt`, `cdr2_plddt`, `cdr3_plddt`, `vh_vl_orientation_confidence`, `predicted_surface_hydrophobic_patch_score`, `predicted_aggregation_patch_score`, `predicted_charge_patch_score`, and `structural_notes`.

Structural risk interpretation is based only on imported computational metrics or user-provided annotations. AbDev-Lite v0.8 does not run IgFold, ImmuneBuilder, AlphaFold, ColabFold, docking, molecular dynamics, antigen-binding prediction, paratope prediction, or experimental structural validation.

## Output Files

After clicking **Run Analysis**, reports are written to `outputs/`:

- `outputs/abdev_lite_results.xlsx`
- `outputs/abdev_lite_report.html`

The Excel workbook contains:

- `Input_Cleaned`
- `Sequence_QC`
- `Chain_Properties`
- `Numbering_Residues`
- `Region_Summary`
- `Liability_Sites`
- `Liability_Region_Map`
- `Chain_Risk_Scores`
- `Humanness_Results`
- `Antibody_Summary`
- `Candidate_Ranking`
- `Formulation_Features`
- `Expreso_Predictions`
- `Formulation_Recommendations`
- `Structure_Results`
- `Structural_Risk_Summary`

The HTML report can be opened offline in a browser and includes executive metrics, numbering and region summaries, humanness summary, Structural Risk Summary, candidate prioritization, formulation recommendation, antibody summary cards, detailed result tables, and disclaimers.

## Installation

```bash
pip install -r requirements.txt
```

If AbNumber/ANARCI installation fails, install it separately with conda or mamba. Numbering is optional; if unavailable, the rest of the analysis still completes and numbering warnings are written to the report tables.

## Run

```bash
streamlit run app.py
```

Upload the main input file, optionally upload humanness/germline and structure result files, then click **Run Analysis**. The page provides **Download Excel Results** and **Download HTML Report** buttons.

## Example Data

```text
data/example_input.xlsx
data/example_humanness_results.xlsx
data/example_structure_results.xlsx
```

Use `example_input.xlsx` alone for a compatibility smoke test. Use `example_humanness_results.xlsx` and `example_structure_results.xlsx` to test optional v0.8 imports, merge behavior, structural risk summary, and candidate ranking integration.

## Methodology

AbDev-Lite v0.8 performs local rule-based screening and aggregation:

- Input sequence cleaning and QC
- Biopython-based physicochemical property calculation
- Rule-based liability motif scanning
- Optional AbNumber/ANARCI IMGT numbering
- IMGT-based FR/CDR extraction and liability mapping
- Rule-based sequence-level and CDR-adjusted risk scoring
- Optional external humanness/germline result import
- Optional external structure summary import
- Structural confidence, CDR loop confidence, and patch-risk aggregation
- Rule-based candidate prioritization from developability, CDR/FR mapping, optional humanness metrics, and optional structural risk metrics
- Formulation-related feature extraction and recommendation
- Antibody-level aggregation

## Current Limitations

- v0.8 does not run IgFold, AlphaFold, ColabFold, or ImmuneBuilder automatically.
- v0.8 only imports external structure summary results.
- Structural interpretation is computational and should be reviewed.
- No antigen-binding, paratope prediction, docking, or molecular dynamics is performed.
- No experimental structural validation is performed.
- Hydrophobic patch detection in the sequence-level liability scanner remains a sequence-level proxy and does not represent true structural surface patch analysis.
- BsAb analysis is chain-level only and does not evaluate chain pairing, heterodimerization, Fc engineering, linker geometry, or full molecule architecture.
- Full-length sequences are only analyzed for basic sequence-level properties in this MVP.
- Human-likeness assessment is not equivalent to clinical immunogenicity prediction.
- Candidate ranking is rule-based computational triage and does not predict experimental success.
