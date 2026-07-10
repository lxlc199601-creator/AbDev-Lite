"""Streamlit app for AbDev-Lite."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from src.antibody_summary import summarize_antibodies
from src.expreso_adapter import MODEL_NOT_FOUND_WARNING, predict_with_expreso
from src.external_input_builder import build_external_input_package, empty_external_run_plan
from src.external_result_importer import (
    empty_external_results,
    load_external_tool_result,
    merge_external_results,
    standardize_external_results,
)
from src.external_tool_registry import get_available_tools
from src.external_tool_summary import build_external_tool_summary, empty_external_tool_summary
from src.final_assessment import build_executive_decision_summary, build_final_assessment
from src.formulation_features import build_formulation_features
from src.formulation_recommendation import build_formulation_recommendations
from src.humanness import load_humanness_file, merge_humanness_results
from src.input_parser import parse_input
from src.liability_scanner import scan_liabilities
from src.numbering import map_liabilities_to_regions, run_imgt_numbering
from src.physicochemical import calculate_properties_table
from src.prioritization import build_candidate_ranking
from src.report_generator import (
    BSAB_NOTICE,
    FULL_LENGTH_LIMITATION,
    FULL_LENGTH_NOTICE,
    HYDROPHOBIC_PATCH_NOTICE,
    NUMBERING_NOTICE,
    REPORT_DISCLAIMER,
    generate_reports,
)
from src.scoring import score_chains
from src.sequence_classifier import classify_sequences
from src.sequence_qc import run_qc
from src.structural_risk import build_structural_risk_summary
from src.structure_import import (
    annotate_structure_result_matches,
    complete_structural_risk_summary_for_antibodies,
    empty_structure_results_for_antibodies,
    load_structure_results,
    merge_structure_results,
)


EXAMPLE_ROWS = [
    {
        "antibody_id": "Ab001",
        "molecule_format": "mAb",
        "chain_id": "VH",
        "chain_type": "heavy",
        "region_type": "variable_region",
        "sequence_scope": "VH",
        "sequence": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR",
    },
    {
        "antibody_id": "Ab001",
        "molecule_format": "mAb",
        "chain_id": "VL",
        "chain_type": "light",
        "region_type": "variable_region",
        "sequence_scope": "VL",
        "sequence": "DIQMTQSPSSLSASVGDRVTITCRASQSVSSYLAWYQQKPGKAPKLLIYDASNRATGIPARFSGSGSGTDFTLTISSLQPEDFATYYCQQRSNWPLTFGQGTKVEIK",
    },
]


def _priority_counts(candidate_ranking_df: pd.DataFrame | None) -> dict[str, int]:
    priority_counts = {"A": 0, "B": 0, "C": 0, "D": 0}
    if (
        candidate_ranking_df is not None
        and not candidate_ranking_df.empty
        and "final_priority_class" in candidate_ranking_df.columns
    ):
        observed_counts = candidate_ranking_df["final_priority_class"].astype(str).value_counts().to_dict()
        for key in priority_counts:
            priority_counts[key] = int(observed_counts.get(key, 0))
    return priority_counts


def run_pipeline(
    uploaded_file,
    humanness_file=None,
    structure_file=None,
    selected_external_tools: list[str] | None = None,
    external_result_file=None,
    external_result_tool_id: str = "",
) -> tuple[dict[str, pd.DataFrame], Path, Path]:
    """Run the full AbDev-Lite analysis pipeline."""
    input_df = parse_input(uploaded_file, uploaded_file.name)
    classified_df = classify_sequences(input_df)
    qc_df = run_qc(classified_df)
    properties_df = calculate_properties_table(classified_df)
    liabilities_df = scan_liabilities(qc_df)
    numbering_df, region_summary_df = run_imgt_numbering(qc_df)
    liability_region_map_df = map_liabilities_to_regions(liabilities_df, numbering_df)
    scores_df = score_chains(qc_df, liabilities_df, liability_region_map_df)
    summary_df = summarize_antibodies(qc_df, liabilities_df, scores_df, liability_region_map_df)
    humanness_input_df = load_humanness_file(humanness_file)
    humanness_results_df, summary_df = merge_humanness_results(
        scores_df,
        region_summary_df,
        summary_df,
        humanness_input_df,
    )
    structure_input_df = load_structure_results(structure_file)
    structure_results_df = annotate_structure_result_matches(summary_df, pd.DataFrame(), structure_input_df)
    if structure_results_df.empty:
        structure_results_df = empty_structure_results_for_antibodies(summary_df)
    structural_risk_summary_df = build_structural_risk_summary(structure_results_df)
    structural_risk_summary_df = complete_structural_risk_summary_for_antibodies(
        summary_df,
        structural_risk_summary_df,
    )
    summary_df, _ = merge_structure_results(
        summary_df,
        pd.DataFrame(),
        structure_results_df,
    )
    tool_registry_df = get_available_tools()
    external_run_plan_df = build_external_input_package(
        qc_df,
        qc_df,
        summary_df,
        selected_external_tools or [],
        Path("external_inputs"),
    )
    external_results_df = empty_external_results()
    if external_result_file is not None and external_result_tool_id:
        imported_external_df = load_external_tool_result(external_result_file, external_result_tool_id)
        external_results_df = standardize_external_results(imported_external_df, external_result_tool_id)
    summary_df, _ = merge_external_results(summary_df, pd.DataFrame(), external_results_df)
    external_summary_df = build_external_tool_summary(external_results_df)
    if external_summary_df.empty:
        external_summary_df = empty_external_tool_summary(summary_df["antibody_id"].dropna().tolist() if "antibody_id" in summary_df.columns else [])
    candidate_ranking_df = build_candidate_ranking(
        summary_df,
        scores_df,
        liability_region_map_df,
        humanness_results_df,
    )
    formulation_features_df = build_formulation_features(
        summary_df,
        properties_df,
        scores_df,
        liability_region_map_df,
        humanness_results_df,
        candidate_ranking_df,
    )
    expreso_predictions_df = predict_with_expreso(
        formulation_features_df,
        Path("models") / "expreso_lite",
    )
    formulation_recommendations_df = build_formulation_recommendations(
        formulation_features_df,
        expreso_predictions_df,
    )
    final_assessment_df = build_final_assessment(
        summary_df,
        candidate_ranking_df,
        formulation_recommendations_df,
        structural_risk_summary_df,
        external_summary_df,
    )

    results = {
        "Input_Cleaned": classified_df,
        "Sequence_QC": qc_df,
        "Chain_Properties": properties_df,
        "Numbering_Residues": numbering_df,
        "Region_Summary": region_summary_df,
        "Liability_Sites": liabilities_df,
        "Liability_Region_Map": liability_region_map_df,
        "Chain_Risk_Scores": scores_df,
        "Humanness_Results": humanness_results_df,
        "Antibody_Summary": summary_df,
        "Final_Assessment": final_assessment_df,
        "Candidate_Ranking": candidate_ranking_df,
        "Formulation_Features": formulation_features_df,
        "Expreso_Predictions": expreso_predictions_df,
        "Formulation_Recommendations": formulation_recommendations_df,
        "Structure_Results": structure_results_df,
        "Structural_Risk_Summary": structural_risk_summary_df,
        "Tool_Registry": tool_registry_df,
        "External_Tool_Run_Plan": external_run_plan_df,
        "External_Tool_Results": external_results_df,
        "External_Tool_Summary": external_summary_df,
    }
    excel_path, html_path = generate_reports(results, uploaded_file.name)
    return results, excel_path, html_path


def _download_file(path: Path) -> bytes:
    return path.read_bytes()


def main() -> None:
    """Render the Streamlit UI."""
    st.set_page_config(page_title="AbDev-Lite", layout="wide")
    st.title("AbDev-Lite v1.0: Integrated Antibody Developability Screening Platform")
    st.write(
        "A local, lightweight platform for end-to-end antibody variable-region developability triage, "
        "candidate prioritization, final integrated assessment, and Excel/HTML report export."
    )
    st.caption(REPORT_DISCLAIMER)

    with st.expander("Input format", expanded=True):
        st.write("Upload a CSV or XLSX file with these required columns:")
        st.code(
            "antibody_id, molecule_format, chain_id, chain_type, region_type, sequence_scope, sequence",
            language="text",
        )
        st.dataframe(pd.DataFrame(EXAMPLE_ROWS), use_container_width=True)

    uploaded_file = st.file_uploader("Upload antibody variable-region input file", type=["csv", "xlsx", "xls"])
    humanness_file = st.file_uploader(
        "Optional humanness/germline result file",
        type=["csv", "xlsx", "xls"],
        help=(
            "Upload an optional external humanness/germline assessment file if available. "
            "If not provided, AbDev-Lite will run sequence-level developability analysis only."
        ),
    )
    structure_file = st.file_uploader(
        "Optional structure prediction result file",
        type=["csv", "xlsx", "xls"],
        help=(
            "Upload an optional external structure prediction summary file if available. "
            "If not provided, AbDev-Lite will run sequence-level, humanness, prioritization, and formulation analysis only."
        ),
    )
    tool_registry = get_available_tools()
    st.subheader("External Tool Adapter")
    st.info(
        "Browser automation is disabled by default. AbDev-Lite v1.0 prepares input files and imports external results, "
        "but does not automatically submit confidential sequences to third-party websites."
    )
    st.write("Tool Registry")
    st.dataframe(tool_registry, use_container_width=True)
    tool_options = tool_registry["tool_id"].astype(str).tolist() if not tool_registry.empty and "tool_id" in tool_registry.columns else []
    selected_external_tools = st.multiselect(
        "Select external tools to prepare input packages",
        options=tool_options,
    )
    if st.button("Generate External Tool Input Package", disabled=uploaded_file is None or not selected_external_tools):
        try:
            input_df = parse_input(uploaded_file, uploaded_file.name)
            classified_df = classify_sequences(input_df)
            qc_df = run_qc(classified_df)
            temp_run_plan = build_external_input_package(qc_df, qc_df, pd.DataFrame(), selected_external_tools, Path("external_inputs"))
            st.success("External tool input package generated in external_inputs.")
            st.write("External_Tool_Run_Plan")
            st.dataframe(temp_run_plan, use_container_width=True)
        except Exception as exc:
            st.error(f"External input package generation failed: {exc}")
    external_result_tool_id = ""
    if tool_options:
        external_result_tool_id = st.selectbox(
            "External result file tool_id",
            options=tool_options,
            index=0,
        )
    external_result_file = st.file_uploader(
        "Upload external tool result file",
        type=["csv", "xlsx", "xls"],
        help="Upload a user-generated external tool result file for import and evidence integration.",
    )

    st.subheader("Module Status")
    module_status_rows = [
        {"Module": "Numbering module status", "Status": "Available when AbNumber/ANARCI is installed; otherwise optional module skipped during run"},
        {"Module": "Humanness file status", "Status": "Provided" if humanness_file is not None else "Not provided / optional module skipped"},
        {"Module": "Structure file status", "Status": "Provided" if structure_file is not None else "Not provided / optional module skipped"},
        {"Module": "Expreso-lite model status", "Status": "Local model used if present; rule-based fallback if not provided"},
        {"Module": "External tool result status", "Status": "Provided" if external_result_file is not None else "Not provided / optional module skipped"},
        {"Module": "Browser automation status", "Status": "Disabled"},
    ]
    st.dataframe(pd.DataFrame(module_status_rows), use_container_width=True, hide_index=True)

    if st.button("Run Analysis", type="primary", disabled=uploaded_file is None):
        if uploaded_file is None:
            st.warning("Please upload a CSV or XLSX file first.")
            return
        try:
            with st.spinner("Running AbDev-Lite analysis and generating reports..."):
                results, excel_path, html_path = run_pipeline(
                    uploaded_file,
                    humanness_file,
                    structure_file,
                    selected_external_tools,
                    external_result_file,
                    external_result_tool_id,
                )
            st.success("Analysis complete. Excel and HTML reports were generated.")

            summary = results["Antibody_Summary"]
            qc = results["Sequence_QC"]
            liabilities = results["Liability_Sites"]
            region_summary = results["Region_Summary"]
            liability_region_map = results["Liability_Region_Map"]
            humanness_results = results["Humanness_Results"]
            candidate_ranking = results.get("Candidate_Ranking", pd.DataFrame())
            final_assessment = results.get("Final_Assessment", pd.DataFrame())
            executive_decision_summary = build_executive_decision_summary(final_assessment)
            formulation_features = results.get("Formulation_Features", pd.DataFrame())
            expreso_predictions = results.get("Expreso_Predictions", pd.DataFrame())
            formulation_recommendations = results.get("Formulation_Recommendations", pd.DataFrame())
            structure_results = results.get("Structure_Results", pd.DataFrame())
            structural_risk_summary = results.get("Structural_Risk_Summary", pd.DataFrame())
            tool_registry_results = results.get("Tool_Registry", pd.DataFrame())
            external_run_plan = results.get("External_Tool_Run_Plan", pd.DataFrame())
            external_results = results.get("External_Tool_Results", pd.DataFrame())
            external_summary = results.get("External_Tool_Summary", pd.DataFrame())
            if summary.empty or "max_chain_risk_class" not in summary.columns:
                risk_counts = pd.Series({"Low": 0, "Medium": 0, "High": 0})
            else:
                risk_counts = summary["max_chain_risk_class"].value_counts().reindex(["Low", "Medium", "High"], fill_value=0)

            st.subheader("Dashboard")
            dashboard_cols = st.columns(7)
            dashboard_cols[0].metric("Total candidates", executive_decision_summary["total_candidates"])
            dashboard_cols[1].metric("Advance", executive_decision_summary["advance_count"])
            dashboard_cols[2].metric("Advance with review", executive_decision_summary["advance_with_review_count"])
            dashboard_cols[3].metric("Engineering review", executive_decision_summary["engineering_review_count"])
            dashboard_cols[4].metric("Deprioritize / redesign", executive_decision_summary["deprioritize_count"])
            dashboard_cols[5].metric("High-risk candidates", executive_decision_summary["high_risk_candidate_count"])
            dashboard_cols[6].metric("Top ranked candidate", executive_decision_summary["top_ranked_candidate"] or "-")
            st.write(executive_decision_summary["overall_project_summary"])
            st.caption("Confidence level means evidence completeness level only; it is not experimental confidence or success probability.")

            st.subheader("Final_Assessment")
            final_preview_columns = [
                "antibody_id",
                "molecule_format",
                "final_priority_class",
                "final_priority_score",
                "go_no_go_suggestion",
                "confidence_level",
                "major_review_flags",
                "recommended_next_action",
            ]
            final_preview = final_assessment.copy()
            for column in final_preview_columns:
                if column not in final_preview.columns:
                    final_preview[column] = ""
            st.dataframe(final_preview[final_preview_columns].head(200), use_container_width=True)
            if len(final_assessment) > 200:
                st.caption("Showing the first 200 final assessment rows in the app preview. The Excel workbook contains all rows.")

            metric_cols = st.columns(6)
            metric_cols[0].metric("Total antibodies", summary["antibody_id"].nunique() if "antibody_id" in summary.columns else 0)
            metric_cols[1].metric("Total chains", len(qc))
            metric_cols[2].metric("Low antibodies", int(risk_counts["Low"]))
            metric_cols[3].metric("Medium antibodies", int(risk_counts["Medium"]))
            metric_cols[4].metric("High antibodies", int(risk_counts["High"]))
            metric_cols[5].metric("Liability sites", len(liabilities))

            numbering_counts = (
                region_summary["numbering_status"].astype(str).value_counts()
                if not region_summary.empty and "numbering_status" in region_summary.columns
                else pd.Series(dtype=int)
            )
            numbering_cols = st.columns(3)
            numbering_cols[0].metric("Numbered chains", int(numbering_counts.get("PASS", 0) + numbering_counts.get("WARNING", 0)))
            numbering_cols[1].metric("Numbering failed chains", int(numbering_counts.get("FAIL", 0)))
            numbering_cols[2].metric("Numbering skipped chains", int(numbering_counts.get("SKIPPED", 0)))

            region_cols = st.columns(3)
            if liability_region_map.empty or "cdr_or_fr" not in liability_region_map.columns:
                cdr_count = fr_count = unknown_count = 0
            else:
                region_categories = liability_region_map["cdr_or_fr"].astype(str)
                cdr_count = int(region_categories.eq("CDR").sum())
                fr_count = int(region_categories.eq("FR").sum())
                unknown_count = int(region_categories.isin(["Unknown", "Boundary"]).sum())
            region_cols[0].metric("CDR liabilities", cdr_count)
            region_cols[1].metric("FR liabilities", fr_count)
            region_cols[2].metric("Unknown-region liabilities", unknown_count)

            st.subheader("Humanness Summary")
            humanness_matched = (
                humanness_results[humanness_results["merge_status"].astype(str).str.startswith("matched")]
                if not humanness_results.empty and "merge_status" in humanness_results.columns
                else pd.DataFrame()
            )
            humanness_counts = (
                humanness_matched["humanness_risk_class"]
                .astype(str)
                .value_counts()
                .reindex(["Low", "Medium", "High", "Unknown"], fill_value=0)
                if not humanness_matched.empty and "humanness_risk_class" in humanness_matched.columns
                else pd.Series({"Low": 0, "Medium": 0, "High": 0, "Unknown": 0})
            )
            humanness_cols = st.columns(6)
            humanness_cols[0].metric("Humanness available", bool(not humanness_matched.empty))
            humanness_cols[1].metric("Chains with data", len(humanness_matched))
            humanness_cols[2].metric("Low", int(humanness_counts["Low"]))
            humanness_cols[3].metric("Medium", int(humanness_counts["Medium"]))
            humanness_cols[4].metric("High", int(humanness_counts["High"]))
            humanness_cols[5].metric("Unknown", int(humanness_counts["Unknown"]))
            if humanness_results.empty:
                st.info(
                    "No external humanness/germline assessment file was provided. "
                    "Humanness analysis was not performed in this run."
                )
            st.caption("Human-likeness is not equivalent to clinical immunogenicity prediction.")

            st.subheader("Structural Risk Summary")
            structural_counts = (
                structural_risk_summary["structural_risk_class"]
                .astype(str)
                .value_counts()
                .reindex(["Low", "Medium", "High", "Not Available"], fill_value=0)
                if not structural_risk_summary.empty and "structural_risk_class" in structural_risk_summary.columns
                else pd.Series({"Low": 0, "Medium": 0, "High": 0, "Not Available": 0})
            )
            structure_cols = st.columns(4)
            structure_cols[0].metric("Low", int(structural_counts["Low"]))
            structure_cols[1].metric("Medium", int(structural_counts["Medium"]))
            structure_cols[2].metric("High", int(structural_counts["High"]))
            structure_cols[3].metric("Not Available", int(structural_counts["Not Available"]))
            st.write("Structure_Results")
            st.dataframe(structure_results.head(200), use_container_width=True)
            if len(structure_results) > 200:
                st.caption("Showing the first 200 structure result rows in the app preview. The Excel workbook contains all rows.")
            st.write("Structural_Risk_Summary")
            st.dataframe(structural_risk_summary.head(200), use_container_width=True)
            st.caption(
                "Structural risk interpretation is based on imported computational metrics or user annotations only. "
                "This version does not generate or validate 3D structures."
            )

            st.subheader("External Tool Integration")
            st.write("Tool Registry")
            st.dataframe(tool_registry_results, use_container_width=True)
            st.write("External_Tool_Run_Plan")
            st.dataframe(external_run_plan.head(200), use_container_width=True)
            if external_run_plan.empty:
                st.info("No external input package was generated in this run.")
            st.write("External_Tool_Results")
            st.dataframe(external_results.head(200), use_container_width=True)
            if external_results.empty:
                st.info("No external tool result file was imported in this run.")
            st.write("External_Tool_Summary")
            st.dataframe(external_summary.head(200), use_container_width=True)
            st.caption(
                "External tool results are user-provided computational evidence. Review imported results before decision-making."
            )

            st.subheader("Candidate Ranking")
            priority_counts = _priority_counts(candidate_ranking)
            priority_cols = st.columns(4)
            priority_cols[0].metric("A candidates", int(priority_counts["A"]))
            priority_cols[1].metric("B candidates", int(priority_counts["B"]))
            priority_cols[2].metric("C candidates", int(priority_counts["C"]))
            priority_cols[3].metric("D candidates", int(priority_counts["D"]))
            ranking_preview_columns = [
                "antibody_id",
                "molecule_format",
                "final_priority_score",
                "final_priority_class",
                "external_high_risk_flags",
                "external_medium_risk_flags",
                "external_tool_results_available",
                "structural_risk_class",
                "decision_label",
                "review_reason",
                "next_step_recommendation",
            ]
            ranking_preview = candidate_ranking.sort_values(
                "final_priority_score",
                ascending=False,
            ) if not candidate_ranking.empty and "final_priority_score" in candidate_ranking.columns else candidate_ranking
            if ranking_preview.empty:
                st.info("No candidate ranking results were generated in this run.")
                st.dataframe(pd.DataFrame(columns=ranking_preview_columns), use_container_width=True)
            else:
                for column in ranking_preview_columns:
                    if column not in ranking_preview.columns:
                        ranking_preview[column] = ""
                st.dataframe(ranking_preview[ranking_preview_columns].head(200), use_container_width=True)
            if len(candidate_ranking) > 200:
                st.caption("Showing the first 200 candidate ranking rows in the app preview. The Excel workbook contains all rows.")
            st.caption(
                "Candidate prioritization is rule-based computational triage only. It does not predict experimental success."
            )

            st.subheader("Formulation Recommendation")
            expreso_available = (
                bool(expreso_predictions["model_available"].fillna(False).astype(bool).any())
                if not expreso_predictions.empty and "model_available" in expreso_predictions.columns
                else False
            )
            prediction_mode = (
                ", ".join(sorted(expreso_predictions["prediction_mode"].dropna().astype(str).unique().tolist()))
                if not expreso_predictions.empty and "prediction_mode" in expreso_predictions.columns
                else "rule_based_fallback"
            )
            st.write(f"Expreso-lite model available: {expreso_available}")
            st.write(f"Prediction mode: {prediction_mode}")
            if not expreso_available:
                st.info(MODEL_NOT_FOUND_WARNING)
            formulation_counts = (
                formulation_recommendations["formulation_risk_class"]
                .astype(str)
                .value_counts()
                .reindex(["Low", "Medium", "High"], fill_value=0)
                if not formulation_recommendations.empty and "formulation_risk_class" in formulation_recommendations.columns
                else pd.Series({"Low": 0, "Medium": 0, "High": 0})
            )
            formulation_metric_cols = st.columns(3)
            formulation_metric_cols[0].metric("Low formulation risk", int(formulation_counts["Low"]))
            formulation_metric_cols[1].metric("Medium formulation risk", int(formulation_counts["Medium"]))
            formulation_metric_cols[2].metric("High formulation risk", int(formulation_counts["High"]))

            st.write("Formulation_Features")
            st.dataframe(formulation_features.head(200), use_container_width=True)
            if len(formulation_features) > 200:
                st.caption("Showing the first 200 formulation feature rows in the app preview. The Excel workbook contains all rows.")
            st.write("Expreso_Predictions")
            st.dataframe(expreso_predictions.head(200), use_container_width=True)
            if len(expreso_predictions) > 200:
                st.caption("Showing the first 200 Expreso prediction rows in the app preview. The Excel workbook contains all rows.")
            st.write("Formulation_Recommendations")
            st.dataframe(formulation_recommendations.head(200), use_container_width=True)
            if len(formulation_recommendations) > 200:
                st.caption("Showing the first 200 formulation recommendation rows in the app preview. The Excel workbook contains all rows.")
            st.caption(
                "Formulation recommendations are early-stage computational screening signals and do not replace experimental formulation screening."
            )

            st.info(f"Reports saved to: {excel_path} and {html_path}")

            st.subheader("Antibody Summary")
            st.dataframe(summary, use_container_width=True)

            st.subheader("Humanness Results")
            st.dataframe(humanness_results.head(200), use_container_width=True)
            if len(humanness_results) > 200:
                st.caption("Showing the first 200 humanness result rows in the app preview. The Excel workbook contains all rows.")

            st.subheader("Liability Sites")
            st.dataframe(liabilities.head(200), use_container_width=True)
            if len(liabilities) > 200:
                st.caption("Showing the first 200 liability sites in the app preview. The Excel workbook contains all rows.")

            st.subheader("Region Summary")
            st.dataframe(region_summary.head(200), use_container_width=True)
            if len(region_summary) > 200:
                st.caption("Showing the first 200 region summary rows in the app preview. The Excel workbook contains all rows.")

            st.subheader("Liability Region Map")
            st.dataframe(liability_region_map.head(200), use_container_width=True)
            if len(liability_region_map) > 200:
                st.caption("Showing the first 200 mapped liability sites in the app preview. The Excel workbook contains all rows.")

            download_cols = st.columns(2)
            with download_cols[0]:
                st.download_button(
                    "Download Excel Results",
                    data=_download_file(excel_path),
                    file_name="abdev_lite_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )
            with download_cols[1]:
                st.download_button(
                    "Download HTML Report",
                    data=_download_file(html_path),
                    file_name="abdev_lite_report.html",
                    mime="text/html",
                )

            st.subheader("Method Notes")
            st.write(NUMBERING_NOTICE)
            st.write(FULL_LENGTH_NOTICE)
            st.write(HYDROPHOBIC_PATCH_NOTICE)

            st.subheader("Limitations")
            st.write(REPORT_DISCLAIMER)
            st.write(BSAB_NOTICE)
            st.write(FULL_LENGTH_LIMITATION)
        except ValueError as exc:
            st.error(str(exc))
        except Exception as exc:
            st.error(f"Analysis failed: {exc}")

    st.divider()
    st.subheader("v1.0 Scope")
    st.write(NUMBERING_NOTICE)
    st.write(FULL_LENGTH_NOTICE)
    st.write(HYDROPHOBIC_PATCH_NOTICE)
    st.write(BSAB_NOTICE)
    st.write(FULL_LENGTH_LIMITATION)
    st.write(REPORT_DISCLAIMER)


if __name__ == "__main__":
    main()
