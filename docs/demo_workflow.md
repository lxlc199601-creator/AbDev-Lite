# AbDev-Lite v1.0 Demo Workflow

## 1. Start the App

```bash
streamlit run app.py
```

## 2. Upload Main Input

Upload:

```text
data/example_input.xlsx
```

This is the only required file. The app can run without optional files.

## 3. Optional Humanness Upload

Upload:

```text
data/example_humanness_results.xlsx
```

If skipped, the module status shows `Not provided / optional module skipped`.

## 4. Optional Structure Upload

Upload:

```text
data/example_structure_results.xlsx
```

If skipped, structural risk remains `Not Available` where appropriate.

## 5. Optional External Tool Result Upload

Upload:

```text
data/example_external_tool_results.xlsx
```

External tool results are imported only from user-provided files. AbDev-Lite does not upload sequences or run browser automation.

## 6. Run Analysis

Click **Run Analysis**.

## 7. Review Final_Assessment

Start with the `Final_Assessment` table on the dashboard. Key fields:

- `final_priority_class`: A/B/C/D integrated priority class
- `go_no_go_suggestion`: rule-based next decision label
- `confidence_level`: evidence completeness level, not experimental confidence
- `major_review_flags`: main reasons for review or escalation
- `recommended_next_action`: planning-oriented next action

## 8. Download Excel

Click **Download Excel Results**. Confirm the workbook contains `Final_Assessment`.

## 9. Download HTML

Click **Download HTML Report**. Confirm the report contains **Final Integrated Assessment** and **Executive Summary**.

## 10. Interpret A/B/C/D and Go/No-Go

- A: strongest integrated computational profile; default suggestion is `Advance`
- B: generally suitable for progression with targeted review; default suggestion is `Advance with review`
- C: requires engineering review before progression; default suggestion is `Engineering review`
- D: low priority or substantial flags; default suggestion is `Deprioritize / redesign`

The suggestion escalates to at least `Engineering review` if high structural risk, high formulation risk, external high-risk flags, or combined high humanness plus CDR liability flags are present.

AbDev-Lite v1.0 provides rule-based and imported computational screening outputs for antibody variable-region developability assessment. Results are intended for candidate triage, reporting, and planning support only. They do not replace experimental binding, expression, stability, immunogenicity, structural, formulation, or CMC studies.
