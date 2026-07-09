# AbDev-Lite MVP v0.6 Validation

## Validation Date

2026-07-09

## Environment

- Local Windows workspace
- Python environment with dependencies from `requirements.txt`
- Streamlit app entry point: `app.py`
- AbNumber/ANARCI is recommended for IMGT numbering but optional
- BioPhi, IgBLAST, Sapiens, and OASis are not required
- No external web services, webpage automation, third-party sequence upload, structure prediction tools, or deep learning models are required

## Install Command

```bash
pip install -r requirements.txt
```

## Run Command

```bash
streamlit run app.py
```

## Test Input File

```text
data/example_input.xlsx
data/example_humanness_results.xlsx
```

## Expected Outputs

```text
outputs/abdev_lite_results.xlsx
outputs/abdev_lite_report.html
```

Expected Excel sheets:

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

Expected HTML sections:

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

## Verified Functions

- CSV/XLSX input parsing
- Sequence cleaning
- Basic sequence QC
- Biopython-based physicochemical property calculation
- Rule-based liability scanning
- Optional AbNumber/ANARCI IMGT numbering
- FR/CDR region extraction
- Liability-to-region mapping
- Chain-level risk scoring
- CDR-adjusted risk scoring while preserving original risk_score
- Optional humanness/germline result import
- Humanness risk classification
- Human germline summary
- Combined developability + humanness flag
- Candidate prioritization
- Rule-based final priority score
- A/B/C/D priority class
- Decision label
- Review reason
- Next-step recommendation
- Antibody-level summary aggregation
- Excel report generation
- Jinja2-based offline HTML report generation
- Streamlit preview of `Antibody_Summary`
- Streamlit preview of `Liability_Sites`
- Streamlit preview of `Region_Summary`
- Streamlit preview of `Liability_Region_Map`
- Streamlit optional humanness/germline upload component
- Streamlit preview of `Humanness_Results`
- Streamlit preview of `Candidate_Ranking`
- Streamlit Excel download button
- Streamlit HTML download button

## Known Limitations

- This MVP performs sequence-level antibody variable-region developability screening with optional IMGT-based computational CDR/FR assignment.
- If AbNumber/ANARCI fails, sequence-level analysis still completes and numbering-related tables record FAIL or SKIPPED states.
- It does not perform 3D structure prediction.
- It imports optional external humanness/germline assessment results, but does not run BioPhi, IgBLAST, Sapiens, or OASis directly.
- Human-likeness assessment is not equivalent to clinical immunogenicity prediction.
- Candidate ranking is rule-based and intended for triage only.
- Candidate ranking does not predict experimental success.
- Candidate ranking does not replace binding, expression, stability, or immunogenicity assays.
- It does not perform automatic humanization design.
- It does not provide humanization mutation suggestions.
- It does not predict antigen binding.
- It does not perform structural paratope prediction.
- It does not perform experimental validation.
- IMGT-based CDR/FR mapping should be reviewed before experimental decisions.
- Hydrophobic patch detection is only a sequence-level proxy and does not represent true structural surface patch analysis.
- BsAb analysis is chain-level only and does not evaluate chain pairing, heterodimerization, Fc engineering, linker geometry, or full molecule architecture.
- Full-length sequences, if uploaded, are only analyzed for basic sequence-level properties in this MVP.
- The humanness/germline file is optional and missing files must not block v0.4-compatible analysis.

## Manual Test Checklist

- [ ] `pip install -r requirements.txt` works
- [ ] `streamlit run app.py` works
- [ ] `data/example_input.xlsx` exists
- [ ] `data/example_humanness_results.xlsx` exists
- [ ] `example_input.xlsx` upload works
- [ ] App runs without humanness file
- [ ] App runs with `example_humanness_results.xlsx`
- [ ] Run Analysis button works
- [ ] AbNumber import works, or fallback mode works
- [ ] IMGT numbering runs on example VH/VL/VHH sequences when AbNumber/ANARCI is available
- [ ] Excel file is generated
- [ ] HTML report is generated
- [ ] Excel contains expected sheets
- [ ] Candidate_Ranking sheet is generated
- [ ] final_priority_score is generated
- [ ] final_priority_class is generated
- [ ] decision_label is generated
- [ ] review_reason is generated
- [ ] next_step_recommendation is generated
- [ ] Humanness_Results sheet is generated
- [ ] Humanness_Results sheet is generated with headers when no humanness file is uploaded
- [ ] Antibody_Summary contains humanness fields
- [ ] Region_Summary sheet is generated
- [ ] Numbering_Residues sheet is generated
- [ ] Liability_Region_Map sheet is generated
- [ ] Liability sites are mapped to CDR/FR when numbering succeeds
- [ ] CDR-adjusted risk score is generated
- [ ] HTML report shows Numbering Summary
- [ ] HTML report shows Region-Level Liability Summary
- [ ] HTML report shows Humanness Summary
- [ ] HTML report shows Candidate Prioritization section
- [ ] HTML report handles missing humanness file gracefully
- [ ] Antibody_Summary preview is visible
- [ ] Liability_Sites preview is visible
- [ ] Streamlit shows Region_Summary preview
- [ ] Streamlit shows Liability_Region_Map preview
- [ ] Streamlit shows optional humanness upload component
- [ ] Streamlit shows Humanness_Results preview
- [ ] Streamlit shows Candidate Ranking preview
- [ ] Combined developability + humanness flag is generated
- [ ] App runs with humanness file
- [ ] App runs without humanness file
- [ ] v0.5 humanness functions still work
- [ ] v0.4 numbering functions still work
- [ ] Download Excel button works
- [ ] Download HTML button works
- [ ] Report disclaimer is visible
- [ ] BsAb limitation note is visible when BsAb is present
- [ ] Existing v0.3 functions still work if numbering fails
- [ ] Existing v0.4 numbering and CDR/FR mapping still work
