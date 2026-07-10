# AbDev-Lite MVP v0.9 Validation

## Validation Date

2026-07-10

## Environment

- Local Windows workspace
- Streamlit app entry point: `app.py`
- Dependencies from `requirements.txt`
- AbNumber/ANARCI is recommended for IMGT numbering but optional
- Browser automation is disabled by default
- No external webpage call, automatic sequence upload, CAPTCHA bypass, login bypass, credential storage, docking, molecular dynamics, or automatic structure model generation is required

## Run Command

```bash
streamlit run app.py
```

## Test Input Files

```text
data/example_input.xlsx
data/example_humanness_results.xlsx
data/example_structure_results.xlsx
data/example_external_tool_results.xlsx
```

## Expected Outputs

```text
outputs/abdev_lite_results.xlsx
outputs/abdev_lite_report.html
external_inputs/
```

Expected v0.9 Excel sheets:

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
- `Tool_Registry`
- `External_Tool_Run_Plan`
- `External_Tool_Results`
- `External_Tool_Summary`

Expected HTML sections include `External Tool Integration`.

## v0.9 Validation Checklist

- [ ] App runs without external tool selection
- [ ] Tool_Registry sheet is generated
- [ ] External_Tool_Run_Plan sheet is generated
- [ ] External input package can be generated
- [ ] External_Tool_Results sheet is generated
- [ ] External_Tool_Summary sheet is generated
- [ ] Example external tool result file can be imported
- [ ] Antibody_Summary includes external tool fields
- [ ] Candidate_Ranking includes external tool fields
- [ ] HTML report shows External Tool Integration section
- [ ] Streamlit shows Tool Registry and Run Plan
- [ ] Browser automation remains disabled by default
- [ ] v0.8 structure import still works
- [ ] v0.7 formulation module still works
- [ ] v0.6 candidate ranking still works
- [ ] v0.5 humanness import still works
- [ ] v0.4 numbering still works

## Safety Checklist

- [ ] No real antibody sequence input package is committed
- [ ] No real external tool result file is committed
- [ ] No browser cache, cookies, credentials, account password, or API key is committed
- [ ] No third-party webpage is opened automatically
- [ ] No sequence is uploaded automatically
- [ ] No CAPTCHA or login restriction is bypassed

## Known Limitations

- v0.9 creates traceable input packages and imports user-provided results; it does not run most external tools directly.
- Browser automation is a reserved interface only and remains disabled by default.
- Imported external results are computational evidence and require human review before decision-making.
- BsAb analysis remains chain-level and does not evaluate full molecular architecture.
