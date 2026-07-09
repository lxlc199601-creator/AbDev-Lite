"""Streamlit app for AbDev-Lite."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from src.antibody_summary import summarize_antibodies
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


def run_pipeline(uploaded_file, humanness_file=None) -> tuple[dict[str, pd.DataFrame], Path, Path]:
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
    candidate_ranking_df = build_candidate_ranking(
        summary_df,
        scores_df,
        liability_region_map_df,
        humanness_results_df,
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
        "Candidate_Ranking": candidate_ranking_df,
    }
    excel_path, html_path = generate_reports(results, uploaded_file.name)
    return results, excel_path, html_path


def _download_file(path: Path) -> bytes:
    return path.read_bytes()


def main() -> None:
    """Render the Streamlit UI."""
    st.set_page_config(page_title="AbDev-Lite", layout="wide")
    st.title("AbDev-Lite: Antibody Variable-Region Developability Screening")
    st.write(
        "A local lightweight MVP for early antibody VH/VL/VHH/scFv variable-region sequence screening. "
        "It performs input checks, basic sequence QC, physicochemical calculations, sequence-level liability "
        "scanning, chain scoring, antibody-level summaries, and report export."
    )

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

    if st.button("Run Analysis", type="primary", disabled=uploaded_file is None):
        if uploaded_file is None:
            st.warning("Please upload a CSV or XLSX file first.")
            return
        try:
            with st.spinner("Running AbDev-Lite analysis and generating reports..."):
                results, excel_path, html_path = run_pipeline(uploaded_file, humanness_file)
            st.success("Analysis complete. Excel and HTML reports were generated.")

            summary = results["Antibody_Summary"]
            qc = results["Sequence_QC"]
            liabilities = results["Liability_Sites"]
            region_summary = results["Region_Summary"]
            liability_region_map = results["Liability_Region_Map"]
            humanness_results = results["Humanness_Results"]
            candidate_ranking = results.get("Candidate_Ranking", pd.DataFrame())
            if summary.empty or "max_chain_risk_class" not in summary.columns:
                risk_counts = pd.Series({"Low": 0, "Medium": 0, "High": 0})
            else:
                risk_counts = summary["max_chain_risk_class"].value_counts().reindex(["Low", "Medium", "High"], fill_value=0)

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
    st.subheader("MVP Scope")
    st.write(NUMBERING_NOTICE)
    st.write(FULL_LENGTH_NOTICE)
    st.write(HYDROPHOBIC_PATCH_NOTICE)
    st.write(BSAB_NOTICE)
    st.write(FULL_LENGTH_LIMITATION)
    st.write(REPORT_DISCLAIMER)


if __name__ == "__main__":
    main()
