# AbDev-Lite v1.0 Validation

## Validation Date

2026-07-10

## Environment

- Local Windows workspace
- Streamlit entry point: `app.py`
- Dependencies from `requirements.txt`
- AbNumber/ANARCI is optional for IMGT numbering
- Browser automation remains disabled
- No external webpage call or automatic sequence upload is required

## Run Command

```bash
streamlit run app.py
```

## Example Inputs

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
- `Final_Assessment`
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

## v1.0 Validation Checklist

- [ ] App starts successfully
- [ ] App runs with only `example_input.xlsx`
- [ ] App runs with humanness file
- [ ] App runs with structure file
- [ ] App runs with external tool result file
- [ ] App runs with all optional files
- [ ] Excel report is generated
- [ ] HTML report is generated
- [ ] Final_Assessment sheet is generated
- [ ] Final Integrated Assessment section appears in HTML
- [ ] Streamlit dashboard displays final metrics
- [ ] Optional modules can be skipped safely
- [ ] Browser automation remains disabled
- [ ] No external upload occurs
- [ ] v0.9 external tool adapter still works
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

- v1.0 is an integrated computational screening and reporting release, not an experimental decision engine.
- `confidence_level` is evidence completeness only.
- Imported optional results are user-provided evidence and require review.
- Browser automation is a reserved interface only and remains disabled.
