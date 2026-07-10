# AbDev-Lite MVP v0.8 Validation

## Validation Date

2026-07-10

## Environment

- Local Windows workspace
- Streamlit app entry point: `app.py`
- Dependencies from `requirements.txt`
- AbNumber/ANARCI is recommended for IMGT numbering but optional
- No external web services, webpage automation, third-party sequence upload, structure prediction tools, docking, molecular dynamics, or deep learning structure models are required

## Run Command

```bash
streamlit run app.py
```

## Test Input Files

```text
data/example_input.xlsx
data/example_humanness_results.xlsx
data/example_structure_results.xlsx
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
- `Formulation_Features`
- `Expreso_Predictions`
- `Formulation_Recommendations`
- `Structure_Results`
- `Structural_Risk_Summary`

Expected HTML sections:

- Executive Summary
- Numbering Summary
- Region-Level Liability Summary
- Humanness Summary
- Structural Risk Summary
- Risk Distribution
- Candidate Prioritization
- Formulation Recommendation
- Antibody Summary Cards
- Liability Sites Table
- Liability Region Map Table
- Humanness Results Table
- Chain-Level Risk Table
- Method Summary
- Limitations and Disclaimer

## v0.8 Validation Checklist

- [ ] `streamlit run app.py` works
- [ ] App runs without structure file
- [ ] App runs with `data/example_input.xlsx`
- [ ] App runs with `data/example_humanness_results.xlsx`
- [ ] App runs with `data/example_structure_results.xlsx`
- [ ] Structure_Results sheet is generated
- [ ] Structural_Risk_Summary sheet is generated
- [ ] Structure_Results and Structural_Risk_Summary keep headers when no structure file is uploaded
- [ ] No-structure runs show `Not Available` structural status
- [ ] Antibody_Summary contains `structure_available`, `structural_risk_class`, `structural_risk_score`, `structural_review_reason`, and `structural_next_step_recommendation`
- [ ] Candidate_Ranking includes `structure_available`, `structural_risk_class`, and `structural_risk_score`
- [ ] Candidate_Ranking integrates structural risk penalties for Medium and High structural risk
- [ ] Candidate_Ranking review_reason reports high structural risk when present
- [ ] Candidate_Ranking review_reason reports no structure result when structure data are not available
- [ ] HTML report shows Structural Risk Summary
- [ ] HTML report includes the v0.8 structure disclaimer
- [ ] Streamlit shows optional structure prediction result upload component
- [ ] Streamlit shows Structure_Results preview table
- [ ] Streamlit shows Structural_Risk_Summary preview table
- [ ] Streamlit shows structural risk metrics: Low, Medium, High, Not Available
- [ ] v0.7 formulation functions still work
- [ ] v0.6 candidate ranking still works
- [ ] v0.5 humanness optional import still works
- [ ] v0.4 numbering still works
- [ ] Excel download remains available
- [ ] HTML download remains available
- [ ] No IgFold, ImmuneBuilder, AlphaFold, or ColabFold run is triggered
- [ ] No external webpage call or sequence upload is triggered
- [ ] No antigen-binding prediction is performed
- [ ] No docking is performed
- [ ] No molecular dynamics is performed

## Known Limitations

- v0.8 imports external structure summary metrics but does not generate or validate 3D structures.
- Structural risk interpretation is computational and should be reviewed.
- Imported model file names or paths are recorded as user-provided labels; they are not treated as experimental structures.
- Human-likeness assessment is not equivalent to clinical immunogenicity prediction.
- Candidate ranking is rule-based and intended for triage only.
- BsAb analysis remains chain-level and does not evaluate full molecular architecture.
