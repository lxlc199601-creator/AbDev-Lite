# AbDev-Lite

AbDev-Lite is a local lightweight Streamlit MVP for antibody variable-region developability screening. Version `MVP v0.6` analyzes sequence-level signals for VH, VL, VHH, and scFv-derived variable domains, adds optional IMGT numbering through AbNumber/ANARCI, maps liability motifs to FR/CDR regions, imports optional external humanness/germline results, ranks candidates with a rule-based decision matrix, then exports both Excel results and an offline HTML report.

The project is intentionally local and lightweight. It does not call external webpages, BioPhi, IgBLAST, Sapiens, OASis, TAP, IgFold, AlphaFold, ColabFold, or deep learning structure models.

## Current Features

- Variable-region sequence input
- CSV/XLSX upload
- Sequence cleaning
- Basic QC
- Physicochemical properties
- Rule-based liability scanning
- IMGT numbering through AbNumber/ANARCI when available
- FR/CDR region extraction
- Liability-to-region mapping
- CDR-adjusted risk scoring
- Optional external humanness/germline result import
- Humanness risk classification
- Human germline summary
- Combined developability + humanness flag
- Candidate prioritization
- Rule-based final priority score
- A/B/C/D priority class
- Decision label
- Review reason
- Next-step recommendation
- Region-level report
- Chain-level risk scoring
- Antibody-level summary
- Excel export
- HTML report export
- Streamlit interface

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

AbDev-Lite v0.6 can import an optional external humanness/germline assessment file in CSV or XLSX format. This file is not required. If it is not provided, AbDev-Lite runs the original sequence-level developability workflow and still generates the `Humanness_Results` sheet with empty headers.

Recommended humanness/germline columns:

| Column | Meaning |
| --- | --- |
| antibody_id | Antibody identifier matching the main input |
| chain_id | Chain identifier matching the main input when available |
| sequence_scope | Sequence scope matching the main input when chain_id is not sufficient |
| humanness_tool | Source label such as BioPhi_OASis, Sapiens, IgBLAST, Manual, or Other |
| humanness_score | Optional numeric score from the external method |
| humanness_percentile | Optional 0-100 percentile |
| closest_human_germline | Closest human germline label such as IGHV3-23 or IGKV1-39 |
| closest_species | human, mouse, rabbit, camelid, unknown, or another imported label |
| v_gene | Imported V gene call |
| j_gene | Imported J gene call |
| identity_to_human_germline | Optional 0-100 identity value |
| framework_identity | Optional 0-100 framework identity value |
| cdr_identity | Optional 0-100 CDR identity value |
| humanness_notes | User or tool notes |

BioPhi/OASis/Sapiens and IgBLAST can be used outside AbDev-Lite to generate results for import, and may be connected through optional adapters in a later version. They are not mandatory dependencies in v0.6. AbDev-Lite does not automatically open websites, upload sequences to third-party platforms, or call these tools.

Humanness risk classes are conservative screening labels derived from imported percentile, framework identity, human germline identity, and closest-species fields. Human-likeness assessment is not equivalent to clinical immunogenicity prediction. v0.6 does not perform automatic humanization mutation design and does not modify user input sequences.

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

The HTML report can be opened offline in a browser and includes:

- Executive Summary
- Numbering Summary
- Region-Level Liability Summary
- Risk Distribution
- Antibody Summary Cards
- Liability Sites Table
- Liability Region Map Table
- Humanness Summary
- Humanness Results Table
- Candidate Prioritization
- Chain-Level Risk Table
- Method Summary
- Limitations and Disclaimer

## Installation

Basic install:

```bash
pip install -r requirements.txt
```

If AbNumber/ANARCI installation fails, try:

```bash
conda install -c bioconda abnumber
```

or:

```bash
mamba install -c bioconda abnumber
```

Numbering is optional but recommended. If AbNumber or the ANARCI backend is not available, AbDev-Lite will still complete the v0.3 sequence-level analysis and will write numbering warnings to the report tables.

## Run

```bash
streamlit run app.py
```

Open the local Streamlit URL shown in the terminal, upload an input file, and click **Run Analysis**. The page provides both **Download Excel Results** and **Download HTML Report** buttons.

## Example Data

Example input data are provided at:

```text
data/example_input.xlsx
data/example_humanness_results.xlsx
```

Use `example_input.xlsx` alone for a quick compatibility smoke test. Use `example_humanness_results.xlsx` as the optional humanness/germline upload to test v0.6 import, merge, and candidate ranking behavior.

## Methodology

AbDev-Lite v0.6 performs:

- Input sequence cleaning
- Basic sequence QC
- Biopython-based physicochemical property calculation
- Rule-based liability motif scanning
- Optional AbNumber/ANARCI IMGT numbering for variable-region sequences
- IMGT-based computational FR/CDR region extraction
- Liability-to-region mapping
- Rule-based sequence-level risk scoring
- CDR-adjusted risk scoring while preserving the original risk score
- Optional external humanness/germline result import
- Conservative humanness risk classification
- Human germline summary
- Combined developability + humanness flag
- Rule-based candidate prioritization from developability, CDR/FR liability mapping, and optional imported humanness metrics
- Final priority score, A/B/C/D priority class, decision label, review reason, and next-step recommendation
- Antibody-level aggregation

Liability scanning currently covers sequence-level proxies such as deamidation motifs, isomerization motifs, oxidation-prone residues, N-glycosylation motifs, cysteine-count warnings, clipping motifs, acid-sensitive motifs, hydrophobic-segment proxies, and low-complexity repeats.

## Current Limitations

This MVP performs sequence-level antibody variable-region developability screening with optional IMGT-based computational CDR/FR mapping, optional imported humanness/germline result aggregation, and rule-based candidate ranking. IMGT-based CDR/FR mapping should not be treated as experimental validation. Human-likeness assessment is not equivalent to clinical immunogenicity prediction. Candidate ranking is rule-based and intended for triage only. Ranking does not predict experimental success and does not replace binding, expression, stability, or immunogenicity assays. It does not perform 3D structure prediction, structural paratope prediction, humanization design, humanization mutation recommendation, antigen-binding prediction, or experimental validation. Results should be interpreted as computational screening signals, not definitive developability conclusions.

Additional limitations:

- Hydrophobic patch detection is only a sequence-level proxy and does not represent true structural surface patch analysis.
- BsAb analysis is chain-level only and does not evaluate chain pairing, heterodimerization, Fc engineering, linker geometry, or full molecule architecture.
- Full-length sequences, if uploaded, are only analyzed for basic sequence-level properties in this MVP.
- Numbering is currently intended for variable-region sequences; full_length rows are marked as skipped by the numbering module.
- scFv and unusual sequences may fail direct numbering; the rest of the analysis will continue.
- Structural paratope prediction is not included in this version.
- BioPhi, OASis, Sapiens, and IgBLAST are not required dependencies.
- Automatic humanization design may be considered in a later version, but is not implemented in v0.6.

## Roadmap

Future extensions may include:

- Expreso-lite formulation recommendation
- BsAb architecture analysis
