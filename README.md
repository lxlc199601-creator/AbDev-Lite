# AbDev-Lite

AbDev-Lite is a local lightweight Streamlit MVP for antibody variable-region developability screening. Version `MVP v0.9` keeps the v0.8 sequence QC, liability scanning, IMGT numbering, CDR/FR mapping, humanness import, candidate prioritization, formulation recommendation, optional structure import, Excel report, HTML report, and Streamlit interface, then adds an External Tool Adapter Framework.

The v0.9 external tool workflow is manual-assisted and auditable. It prepares local input packages, tracks planned external runs, imports user-provided result files, summarizes imported evidence, and integrates external high/medium risk flags into `Antibody_Summary` and `Candidate_Ranking`.

The project is intentionally local and lightweight. It does not automatically submit sequences to third-party websites, does not bypass CAPTCHA or login permissions, does not store credentials or API keys, and does not enable browser automation by default.

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
- Excel report, HTML report, and Streamlit interface
- External tool registry
- External input package generation
- Manual-assisted external tool workflow
- External result import
- External tool summary
- Integration of external tool results into `Antibody_Summary` and `Candidate_Ranking`
- Browser automation interface reserved but disabled by default

## External Tool Adapter Framework

The registry is defined in `tool_specs/external_tools.yaml` and currently includes:

- BioPhi / OASis
- IgBLAST
- TAP
- IgFold
- ImmuneBuilder
- Custom Manual Tool

Generated local external input packages are written under `external_inputs/`. User-provided external result files can be imported from CSV or XLSX. Runtime folders `external_inputs/`, `external_results/`, and `external_runs/` are ignored by Git except for `.gitkeep` files, so real project sequences and external run results should not be committed.

Browser automation is represented only by a reserved adapter interface in `src/browser_adapter.py`. It is disabled by default and does not install Playwright or Selenium, open webpages, submit sequences, bypass CAPTCHA, bypass login requirements, save credentials, or store API keys.

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

## Optional Imports

AbDev-Lite can import optional humanness/germline results, optional structure prediction summary results, and v0.9 external tool results. If any optional file is not provided, the core sequence-level analysis still runs and the corresponding output sheets keep headers and empty/default status rows.

External tool result imports use the standard `External_Tool_Results` fields, including `tool_id`, `antibody_id`, `chain_id`, `sequence_scope`, `result_metric_name`, `result_metric_value`, `result_risk_class`, and `result_interpretation`. Missing fields are allowed but are marked with warnings.

## Output Files

After clicking **Run Analysis**, reports are written to `outputs/`:

- `outputs/abdev_lite_results.xlsx`
- `outputs/abdev_lite_report.html`

The Excel workbook contains the v0.8 sheets plus:

- `Tool_Registry`
- `External_Tool_Run_Plan`
- `External_Tool_Results`
- `External_Tool_Summary`

The HTML report includes an **External Tool Integration** section with registry, run plan, imported result, antibody-level summary, and disclaimer tables.

## Installation

```bash
pip install -r requirements.txt
```

If AbNumber/ANARCI installation fails, install it separately with conda or mamba. Numbering is optional; if unavailable, the rest of the analysis still completes and numbering warnings are written to the report tables.

## Run

```bash
streamlit run app.py
```

Upload the main input file, optionally select external tools to generate input packages, optionally upload humanness/germline, structure, or external tool result files, then click **Run Analysis**. The page provides **Download Excel Results** and **Download HTML Report** buttons.

## Example Data

```text
data/example_input.xlsx
data/example_humanness_results.xlsx
data/example_structure_results.xlsx
data/example_external_tool_results.xlsx
```

Use `example_input.xlsx` alone for a compatibility smoke test. Use `example_external_tool_results.xlsx` to test v0.9 external result import and candidate ranking integration.

## Limitations

- v0.9 does not automatically submit sequences to websites.
- Browser automation is disabled by default.
- Users must comply with third-party tool terms of use.
- Users should not upload confidential antibody sequences to public tools without authorization.
- External results are imported as user-provided computational evidence and should be reviewed.
- Web automation may be considered in a later version only for tools that explicitly permit it.
- No antigen-binding prediction, paratope prediction, docking, molecular dynamics, or experimental validation is performed.
- BsAb analysis remains chain-level and does not evaluate chain pairing, heterodimerization, Fc engineering, linker geometry, or full molecule architecture.
- Human-likeness assessment is not equivalent to clinical immunogenicity prediction.
