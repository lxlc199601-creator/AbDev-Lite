"""Excel and HTML report generation."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from .utils import ensure_output_dir


REPORT_TITLE = "AbDev-Lite: Antibody Variable-Region Developability Screening Report"
REPORT_VERSION = "MVP v0.7"
ANALYSIS_SCOPE = (
    "This MVP focuses on sequence-level developability screening of antibody variable regions including VH, VL, "
    "VHH, and scFv-derived variable domains."
)
REPORT_DISCLAIMER = (
    "This MVP performs sequence-level antibody variable-region developability screening with optional IMGT-based "
    "computational numbering, CDR/FR region assignment, optional imported humanness/germline assessment, and "
    "early-stage formulation recommendation integration. It also provides rule-based candidate prioritization "
    "for triage support. It does not perform 3D structure "
    "prediction, humanization design, antigen-binding prediction, structural paratope prediction, or experimental "
    "validation. Human-likeness assessment is not equivalent to clinical immunogenicity prediction. Results should "
    "be interpreted as computational screening signals, not definitive developability conclusions."
)
FORMULATION_DISCLAIMER = (
    "Formulation recommendations are computational triage outputs based on sequence-level features and optional "
    "Expreso-lite model predictions. They are intended to support early preformulation planning only and do not "
    "replace experimental formulation screening, stability studies, or CMC evaluation."
)
PRIORITIZATION_DISCLAIMER = (
    "Candidate prioritization is based on rule-based computational scoring from sequence-level developability, "
    "CDR/FR mapping, and optional imported humanness metrics. It should be used for triage and reporting support "
    "only, not as a definitive experimental or clinical developability conclusion."
)
HYDROPHOBIC_PATCH_NOTICE = (
    "Hydrophobic patch detection is only a sequence-level proxy and does not represent true structural surface patch "
    "analysis."
)
BSAB_NOTICE = (
    "BsAb analysis is chain-level only and does not evaluate chain pairing, heterodimerization, Fc engineering, "
    "linker geometry, or full molecule architecture."
)
FULL_LENGTH_LIMITATION = (
    "Full-length sequences, if uploaded, are only analyzed for basic sequence-level properties in this MVP."
)

NUMBERING_NOTICE = (
    "This version adds IMGT-based computational numbering and CDR/FR region assignment for antibody variable regions. "
    "Region mapping depends on successful numbering and should be reviewed before experimental decisions. This MVP "
    "does not perform structural paratope prediction or binding-impact prediction."
)

FULL_LENGTH_NOTICE = (
    "Current version primarily supports variable-region developability screening. Full-length analysis is limited "
    "to basic sequence checks and physicochemical reference metrics."
)


SHEET_ORDER = [
    "Input_Cleaned",
    "Sequence_QC",
    "Chain_Properties",
    "Numbering_Residues",
    "Region_Summary",
    "Liability_Sites",
    "Liability_Region_Map",
    "Chain_Risk_Scores",
    "Humanness_Results",
    "Antibody_Summary",
    "Candidate_Ranking",
    "Formulation_Features",
    "Expreso_Predictions",
    "Formulation_Recommendations",
]


def _df(value: pd.DataFrame | None) -> pd.DataFrame:
    return value if isinstance(value, pd.DataFrame) else pd.DataFrame()


def _now_string(analysis_time: datetime | str | None = None) -> str:
    if isinstance(analysis_time, datetime):
        return analysis_time.strftime("%Y-%m-%d %H:%M:%S")
    if analysis_time:
        return str(analysis_time)
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _safe_records(df: pd.DataFrame, columns: list[str] | None = None, limit: int | None = None) -> list[dict[str, object]]:
    table = _df(df).copy()
    if columns is not None:
        for column in columns:
            if column not in table.columns:
                table[column] = ""
        table = table[columns]
    if limit is not None:
        table = table.head(limit)
    if table.empty:
        return []
    table = table.where(pd.notna(table), "")
    return table.to_dict("records")


def _series_value_counts(df: pd.DataFrame, column: str, labels: list[str]) -> dict[str, int]:
    counts_by_label = {label: 0 for label in labels}
    if df.empty or column not in df.columns:
        return counts_by_label
    counts = df[column].astype(str).value_counts()
    for label in labels:
        counts_by_label[label] = int(counts.get(label, 0))
    return counts_by_label


def _priority_counts(candidate_ranking_df: pd.DataFrame | None) -> dict[str, int]:
    priority_counts = {"A": 0, "B": 0, "C": 0, "D": 0}
    ranking = _df(candidate_ranking_df)
    if not ranking.empty and "final_priority_class" in ranking.columns:
        observed_counts = ranking["final_priority_class"].astype(str).value_counts().to_dict()
        for key in priority_counts:
            priority_counts[key] = int(observed_counts.get(key, 0))
    return priority_counts


def _summary_with_ranking(antibody_summary_df: pd.DataFrame, candidate_ranking_df: pd.DataFrame) -> pd.DataFrame:
    summary = _df(antibody_summary_df).copy()
    ranking = _df(candidate_ranking_df)
    if summary.empty or ranking.empty or "antibody_id" not in summary.columns or "antibody_id" not in ranking.columns:
        return summary
    ranking_columns = ["antibody_id", "final_priority_class", "final_priority_score", "decision_label"]
    for column in ranking_columns:
        if column not in ranking.columns:
            ranking[column] = ""
    summary = summary.drop(columns=[column for column in ranking_columns[1:] if column in summary.columns])
    return summary.merge(ranking[ranking_columns], on="antibody_id", how="left")


def _summary_with_formulation(summary_df: pd.DataFrame, formulation_recommendations_df: pd.DataFrame) -> pd.DataFrame:
    summary = _df(summary_df).copy()
    formulation = _df(formulation_recommendations_df)
    if summary.empty or formulation.empty or "antibody_id" not in summary.columns or "antibody_id" not in formulation.columns:
        return summary
    formulation_columns = [
        "antibody_id",
        "formulation_risk_class",
        "recommended_excipient_classes",
        "formulation_next_step_recommendation",
    ]
    for column in formulation_columns:
        if column not in formulation.columns:
            formulation[column] = ""
    summary = summary.drop(columns=[column for column in formulation_columns[1:] if column in summary.columns])
    return summary.merge(formulation[formulation_columns], on="antibody_id", how="left")


def _contains_bsab(*frames: pd.DataFrame) -> bool:
    for frame in frames:
        if not frame.empty and "molecule_format" in frame.columns:
            if frame["molecule_format"].astype(str).str.lower().eq("bsab").any():
                return True
    return False


def _build_results(
    input_cleaned_df: pd.DataFrame | None = None,
    qc_df: pd.DataFrame | None = None,
    properties_df: pd.DataFrame | None = None,
    numbering_df: pd.DataFrame | None = None,
    region_summary_df: pd.DataFrame | None = None,
    liability_df: pd.DataFrame | None = None,
    liability_region_map_df: pd.DataFrame | None = None,
    chain_scores_df: pd.DataFrame | None = None,
    humanness_df: pd.DataFrame | None = None,
    antibody_summary_df: pd.DataFrame | None = None,
    candidate_ranking_df: pd.DataFrame | None = None,
    formulation_features_df: pd.DataFrame | None = None,
    expreso_predictions_df: pd.DataFrame | None = None,
    formulation_recommendations_df: pd.DataFrame | None = None,
) -> dict[str, pd.DataFrame]:
    return {
        "Input_Cleaned": _df(input_cleaned_df),
        "Sequence_QC": _df(qc_df),
        "Chain_Properties": _df(properties_df),
        "Numbering_Residues": _df(numbering_df),
        "Region_Summary": _df(region_summary_df),
        "Liability_Sites": _df(liability_df),
        "Liability_Region_Map": _df(liability_region_map_df),
        "Chain_Risk_Scores": _df(chain_scores_df),
        "Humanness_Results": _df(humanness_df),
        "Antibody_Summary": _df(antibody_summary_df),
        "Candidate_Ranking": _df(candidate_ranking_df),
        "Formulation_Features": _df(formulation_features_df),
        "Expreso_Predictions": _df(expreso_predictions_df),
        "Formulation_Recommendations": _df(formulation_recommendations_df),
    }


def generate_excel_report(
    input_cleaned_df: pd.DataFrame | None = None,
    qc_df: pd.DataFrame | None = None,
    properties_df: pd.DataFrame | None = None,
    numbering_df: pd.DataFrame | None = None,
    region_summary_df: pd.DataFrame | None = None,
    liability_df: pd.DataFrame | None = None,
    liability_region_map_df: pd.DataFrame | None = None,
    chain_scores_df: pd.DataFrame | None = None,
    antibody_summary_df: pd.DataFrame | None = None,
    output_path: str | Path = "outputs/abdev_lite_results.xlsx",
    humanness_df: pd.DataFrame | None = None,
    candidate_ranking_df: pd.DataFrame | None = None,
    formulation_features_df: pd.DataFrame | None = None,
    expreso_predictions_df: pd.DataFrame | None = None,
    formulation_recommendations_df: pd.DataFrame | None = None,
) -> Path:
    """Write all result tables to a multi-sheet Excel workbook."""
    output_path = Path(output_path)
    ensure_output_dir(output_path.parent)
    results = _build_results(
        input_cleaned_df,
        qc_df,
        properties_df,
        numbering_df,
        region_summary_df,
        liability_df,
        liability_region_map_df,
        chain_scores_df,
        humanness_df,
        antibody_summary_df,
        candidate_ranking_df,
        formulation_features_df,
        expreso_predictions_df,
        formulation_recommendations_df,
    )
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for sheet in SHEET_ORDER:
            results[sheet].to_excel(writer, sheet_name=sheet, index=False)
    return output_path


def write_excel(results: dict[str, pd.DataFrame], output_dir: str | Path = "outputs") -> Path:
    """Write all result tables to a multi-sheet Excel workbook."""
    output_path = ensure_output_dir(output_dir) / "abdev_lite_results.xlsx"
    return generate_excel_report(
        results.get("Input_Cleaned"),
        results.get("Sequence_QC"),
        results.get("Chain_Properties"),
        results.get("Numbering_Residues"),
        results.get("Region_Summary"),
        results.get("Liability_Sites"),
        results.get("Liability_Region_Map"),
        results.get("Chain_Risk_Scores"),
        results.get("Antibody_Summary"),
        output_path,
        humanness_df=results.get("Humanness_Results"),
        candidate_ranking_df=results.get("Candidate_Ranking"),
        formulation_features_df=results.get("Formulation_Features"),
        expreso_predictions_df=results.get("Expreso_Predictions"),
        formulation_recommendations_df=results.get("Formulation_Recommendations"),
    )


def generate_html_report(
    input_cleaned_df: pd.DataFrame | None = None,
    qc_df: pd.DataFrame | None = None,
    properties_df: pd.DataFrame | None = None,
    numbering_df: pd.DataFrame | None = None,
    region_summary_df: pd.DataFrame | None = None,
    liability_df: pd.DataFrame | None = None,
    liability_region_map_df: pd.DataFrame | None = None,
    chain_scores_df: pd.DataFrame | None = None,
    antibody_summary_df: pd.DataFrame | None = None,
    output_path: str | Path = "outputs/abdev_lite_report.html",
    uploaded_filename: str = "",
    analysis_time: datetime | str | None = None,
    template_dir: str | Path = "templates",
    humanness_df: pd.DataFrame | None = None,
    candidate_ranking_df: pd.DataFrame | None = None,
    formulation_features_df: pd.DataFrame | None = None,
    expreso_predictions_df: pd.DataFrame | None = None,
    formulation_recommendations_df: pd.DataFrame | None = None,
) -> Path:
    """Render the HTML report with Jinja2."""
    output_path = Path(output_path)
    ensure_output_dir(output_path.parent)

    input_cleaned_df = _df(input_cleaned_df)
    qc_df = _df(qc_df)
    properties_df = _df(properties_df)
    numbering_df = _df(numbering_df)
    region_summary_df = _df(region_summary_df)
    liability_df = _df(liability_df)
    liability_region_map_df = _df(liability_region_map_df)
    chain_scores_df = _df(chain_scores_df)
    humanness_df = _df(humanness_df)
    antibody_summary_df = _df(antibody_summary_df)
    candidate_ranking_df = _df(candidate_ranking_df)
    formulation_features_df = _df(formulation_features_df)
    expreso_predictions_df = _df(expreso_predictions_df)
    formulation_recommendations_df = _df(formulation_recommendations_df)
    summary_display_df = _summary_with_ranking(antibody_summary_df, candidate_ranking_df)
    summary_display_df = _summary_with_formulation(summary_display_df, formulation_recommendations_df)

    risk_counts = _series_value_counts(antibody_summary_df, "max_chain_risk_class", ["Low", "Medium", "High"])
    total_antibodies = (
        int(antibody_summary_df["antibody_id"].nunique())
        if not antibody_summary_df.empty and "antibody_id" in antibody_summary_df.columns
        else 0
    )
    total_chains = int(len(qc_df))
    variable_count = (
        int(qc_df["region_type"].astype(str).str.lower().eq("variable_region").sum())
        if not qc_df.empty and "region_type" in qc_df.columns
        else 0
    )
    full_length_count = (
        int(qc_df["region_type"].astype(str).str.lower().eq("full_length").sum())
        if not qc_df.empty and "region_type" in qc_df.columns
        else 0
    )
    common_liabilities = (
        "; ".join(f"{risk_type} ({count})" for risk_type, count in liability_df["risk_type"].value_counts().head(5).items())
        if not liability_df.empty and "risk_type" in liability_df.columns
        else "None detected"
    )
    liability_display_count = min(len(liability_df), 200)
    region_map_display_count = min(len(liability_region_map_df), 200)
    numbering_status_counts = _series_value_counts(region_summary_df, "numbering_status", ["PASS", "WARNING", "FAIL", "SKIPPED"])
    numbering_success_count = numbering_status_counts["PASS"] + numbering_status_counts["WARNING"]
    numbering_failed_count = numbering_status_counts["FAIL"]
    numbering_skipped_count = numbering_status_counts["SKIPPED"]
    total_numbered_chains = numbering_success_count
    total_cdr_liabilities = (
        int(liability_region_map_df["cdr_or_fr"].astype(str).eq("CDR").sum())
        if not liability_region_map_df.empty and "cdr_or_fr" in liability_region_map_df.columns
        else 0
    )
    total_fr_liabilities = (
        int(liability_region_map_df["cdr_or_fr"].astype(str).eq("FR").sum())
        if not liability_region_map_df.empty and "cdr_or_fr" in liability_region_map_df.columns
        else 0
    )
    total_unknown_region_liabilities = (
        int(liability_region_map_df["cdr_or_fr"].astype(str).isin(["Unknown", "Boundary"]).sum())
        if not liability_region_map_df.empty and "cdr_or_fr" in liability_region_map_df.columns
        else 0
    )
    cdr_rows = (
        liability_region_map_df[liability_region_map_df["cdr_or_fr"].astype(str).eq("CDR")]
        if not liability_region_map_df.empty and "cdr_or_fr" in liability_region_map_df.columns
        else pd.DataFrame()
    )
    fr_rows = (
        liability_region_map_df[liability_region_map_df["cdr_or_fr"].astype(str).eq("FR")]
        if not liability_region_map_df.empty and "cdr_or_fr" in liability_region_map_df.columns
        else pd.DataFrame()
    )
    most_common_cdr = (
        "; ".join(f"{risk_type} ({count})" for risk_type, count in cdr_rows["risk_type"].value_counts().head(5).items())
        if not cdr_rows.empty and "risk_type" in cdr_rows.columns
        else "None detected"
    )
    most_common_fr = (
        "; ".join(f"{risk_type} ({count})" for risk_type, count in fr_rows["risk_type"].value_counts().head(5).items())
        if not fr_rows.empty and "risk_type" in fr_rows.columns
        else "None detected"
    )
    humanness_matched = (
        humanness_df[humanness_df["merge_status"].astype(str).str.startswith("matched")]
        if not humanness_df.empty and "merge_status" in humanness_df.columns
        else pd.DataFrame()
    )
    humanness_available = bool(not humanness_matched.empty)
    humanness_counts = _series_value_counts(
        humanness_matched,
        "humanness_risk_class",
        ["Low", "Medium", "High", "Unknown"],
    )
    priority_counts = _priority_counts(candidate_ranking_df)
    formulation_counts = _series_value_counts(
        formulation_recommendations_df,
        "formulation_risk_class",
        ["Low", "Medium", "High"],
    )
    expreso_model_available = (
        bool(expreso_predictions_df["model_available"].fillna(False).astype(bool).any())
        if not expreso_predictions_df.empty and "model_available" in expreso_predictions_df.columns
        else False
    )
    prediction_modes = (
        ", ".join(sorted(expreso_predictions_df["prediction_mode"].dropna().astype(str).unique().tolist()))
        if not expreso_predictions_df.empty and "prediction_mode" in expreso_predictions_df.columns
        else "rule_based_fallback"
    )
    common_germlines = (
        "; ".join(
            f"{germline} ({count})"
            for germline, count in humanness_matched["closest_human_germline"].replace("", pd.NA).dropna().value_counts().head(5).items()
        )
        if not humanness_matched.empty and "closest_human_germline" in humanness_matched.columns
        else "None available"
    )

    context = {
        "title": REPORT_TITLE,
        "project_name": "AbDev-Lite",
        "uploaded_filename": uploaded_filename or "-",
        "analysis_time": _now_string(analysis_time),
        "report_version": REPORT_VERSION,
        "analysis_scope": ANALYSIS_SCOPE,
        "executive_summary": {
            "total_antibodies": total_antibodies,
            "total_chains": total_chains,
            "variable_region_chain_count": variable_count,
            "full_length_chain_count": full_length_count,
            "low_risk_antibody_count": risk_counts["Low"],
            "medium_risk_antibody_count": risk_counts["Medium"],
            "high_risk_antibody_count": risk_counts["High"],
            "total_liability_sites": int(len(liability_df)),
            "most_common_liability_types": common_liabilities,
        },
        "risk_counts": risk_counts,
        "total_antibodies": max(total_antibodies, 1),
        "numbering_summary": {
            "total_numbered_chains": total_numbered_chains,
            "numbering_success_count": numbering_success_count,
            "numbering_failed_count": numbering_failed_count,
            "numbering_skipped_count": numbering_skipped_count,
            "numbering_scheme": "IMGT",
        },
        "region_liability_summary": {
            "total_cdr_liabilities": total_cdr_liabilities,
            "total_fr_liabilities": total_fr_liabilities,
            "total_unknown_region_liabilities": total_unknown_region_liabilities,
            "most_common_cdr_liability_types": most_common_cdr,
            "most_common_fr_liability_types": most_common_fr,
        },
        "humanness_summary": {
            "humanness_available": humanness_available,
            "number_of_chains_with_humanness_data": int(len(humanness_matched)),
            "low_humanness_risk_chain_count": humanness_counts["Low"],
            "medium_humanness_risk_chain_count": humanness_counts["Medium"],
            "high_humanness_risk_chain_count": humanness_counts["High"],
            "unknown_humanness_risk_chain_count": humanness_counts["Unknown"],
            "most_common_closest_human_germlines": common_germlines,
        },
        "priority_counts": priority_counts,
        "candidate_ranking_rows": _safe_records(
            candidate_ranking_df,
            [
                "antibody_id",
                "molecule_format",
                "final_priority_score",
                "final_priority_class",
                "decision_label",
                "review_reason",
                "next_step_recommendation",
            ],
        ),
        "formulation_summary": {
            "total_antibodies_with_formulation_recommendation": int(len(formulation_recommendations_df)),
            "low_formulation_risk_count": formulation_counts["Low"],
            "medium_formulation_risk_count": formulation_counts["Medium"],
            "high_formulation_risk_count": formulation_counts["High"],
            "expreso_model_available": expreso_model_available,
            "prediction_mode": prediction_modes,
        },
        "formulation_recommendation_rows": _safe_records(
            formulation_recommendations_df,
            [
                "antibody_id",
                "formulation_risk_class",
                "formulation_risk_score",
                "recommended_excipient_classes",
                "buffer_ph_direction",
                "surfactant_consideration",
                "sugar_polyol_consideration",
                "oxidation_control_consideration",
                "formulation_next_step_recommendation",
            ],
        ),
        "summary_rows": _safe_records(summary_display_df),
        "liability_rows": _safe_records(
            liability_df,
            [
                "antibody_id",
                "chain_id",
                "sequence_scope",
                "risk_type",
                "motif",
                "position_start",
                "position_end",
                "local_sequence_window",
                "risk_level",
                "explanation",
            ],
            200,
        ),
        "liability_total_count": int(len(liability_df)),
        "liability_display_count": liability_display_count,
        "liability_region_rows": _safe_records(
            liability_region_map_df,
            [
                "antibody_id",
                "chain_id",
                "sequence_scope",
                "risk_type",
                "motif",
                "position_start",
                "position_end",
                "mapped_regions",
                "cdr_or_fr",
                "imgt_positions",
                "risk_level",
            ],
            200,
        ),
        "liability_region_total_count": int(len(liability_region_map_df)),
        "liability_region_display_count": region_map_display_count,
        "chain_score_rows": _safe_records(
            chain_scores_df,
            [
                "antibody_id",
                "chain_id",
                "sequence_scope",
                "risk_score",
                "risk_class",
                "cdr_adjusted_risk_score",
                "cdr_adjusted_risk_class",
                "cdr_liability_count",
                "fr_liability_count",
                "main_chain_liabilities",
                "chain_recommendation",
            ],
        ),
        "humanness_rows": _safe_records(
            humanness_df,
            [
                "antibody_id",
                "chain_id",
                "sequence_scope",
                "humanness_tool",
                "humanness_score",
                "humanness_percentile",
                "closest_human_germline",
                "closest_species",
                "framework_identity",
                "humanness_risk_class",
                "humanness_interpretation",
            ],
        ),
        "has_bsab": _contains_bsab(input_cleaned_df, qc_df, antibody_summary_df),
        "disclaimer": REPORT_DISCLAIMER,
        "prioritization_disclaimer": PRIORITIZATION_DISCLAIMER,
        "formulation_disclaimer": FORMULATION_DISCLAIMER,
        "numbering_notice": NUMBERING_NOTICE,
        "full_length_notice": FULL_LENGTH_NOTICE,
        "hydrophobic_patch_notice": HYDROPHOBIC_PATCH_NOTICE,
        "bsab_notice": BSAB_NOTICE,
        "full_length_limitation": FULL_LENGTH_LIMITATION,
    }

    env = Environment(loader=FileSystemLoader(str(template_dir)), autoescape=select_autoescape())
    env.filters["dash"] = lambda value: "-" if value is None or str(value).strip() == "" or str(value).lower() == "nan" else value
    template = env.get_template("report_template.html")
    output_path.write_text(template.render(**context), encoding="utf-8")
    return output_path


def write_html_report(
    results: dict[str, pd.DataFrame],
    upload_filename: str,
    template_dir: str | Path = "templates",
    output_dir: str | Path = "outputs",
) -> Path:
    """Compatibility wrapper for older callers."""
    return generate_html_report(
        results.get("Input_Cleaned"),
        results.get("Sequence_QC"),
        results.get("Chain_Properties"),
        results.get("Numbering_Residues"),
        results.get("Region_Summary"),
        results.get("Liability_Sites"),
        results.get("Liability_Region_Map"),
        results.get("Chain_Risk_Scores"),
        results.get("Antibody_Summary"),
        ensure_output_dir(output_dir) / "abdev_lite_report.html",
        upload_filename,
        template_dir=template_dir,
        humanness_df=results.get("Humanness_Results"),
        candidate_ranking_df=results.get("Candidate_Ranking"),
        formulation_features_df=results.get("Formulation_Features"),
        expreso_predictions_df=results.get("Expreso_Predictions"),
        formulation_recommendations_df=results.get("Formulation_Recommendations"),
    )


def generate_reports(results: dict[str, pd.DataFrame], upload_filename: str) -> tuple[Path, Path]:
    """Generate both Excel and HTML reports."""
    excel_path = write_excel(results)
    html_path = write_html_report(results, upload_filename)
    return excel_path, html_path
