"""Final integrated antibody assessment for AbDev-Lite v1.0."""

from __future__ import annotations

import pandas as pd

from .utils import join_messages, normalize_text


FINAL_ASSESSMENT_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "final_priority_class",
    "final_priority_score",
    "decision_label",
    "developability_risk_class",
    "humanness_risk_class",
    "structural_risk_class",
    "formulation_risk_class",
    "external_tool_risk_flag",
    "key_strengths",
    "key_weaknesses",
    "major_review_flags",
    "recommended_next_action",
    "go_no_go_suggestion",
    "confidence_level",
    "final_interpretation",
]

GO_NO_GO_ORDER = {
    "Advance": 0,
    "Advance with review": 1,
    "Engineering review": 2,
    "Deprioritize / redesign": 3,
}


def _df(value: pd.DataFrame | None) -> pd.DataFrame:
    return value if isinstance(value, pd.DataFrame) else pd.DataFrame()


def _safe_float(value: object, default: float = 0.0) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return default
    if pd.isna(numeric):
        return default
    return numeric


def _safe_int(value: object, default: int = 0) -> int:
    return int(round(_safe_float(value, float(default))))


def _safe_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = normalize_text(value).lower()
    return text in {"true", "1", "yes", "y", "available", "matched"}


def _first_match(df: pd.DataFrame, antibody_id: object) -> pd.Series:
    if df.empty or "antibody_id" not in df.columns:
        return pd.Series(dtype=object)
    matches = df[df["antibody_id"] == antibody_id]
    if matches.empty:
        return pd.Series(dtype=object)
    return matches.iloc[0]


def _merge_unique(left: pd.DataFrame, right: pd.DataFrame, suffix: str) -> pd.DataFrame:
    if right.empty or "antibody_id" not in right.columns:
        return left
    right = right.drop_duplicates(subset=["antibody_id"], keep="first").copy()
    overlap = [column for column in right.columns if column != "antibody_id" and column in left.columns]
    right = right.rename(columns={column: f"{column}_{suffix}" for column in overlap})
    return left.merge(right, on="antibody_id", how="left")


def _base_go_no_go(priority_class: str) -> str:
    return {
        "A": "Advance",
        "B": "Advance with review",
        "C": "Engineering review",
        "D": "Deprioritize / redesign",
    }.get(priority_class, "Engineering review")


def _at_least_engineering_review(suggestion: str) -> str:
    if GO_NO_GO_ORDER.get(suggestion, 2) < GO_NO_GO_ORDER["Engineering review"]:
        return "Engineering review"
    return suggestion


def _risk_flag(row: pd.Series) -> str:
    high = _safe_int(row.get("external_high_risk_flags"))
    medium = _safe_int(row.get("external_medium_risk_flags"))
    available = _safe_bool(row.get("external_tool_results_available"))
    if high > 0:
        return "High"
    if medium > 0:
        return "Medium"
    if available:
        return "Low"
    return "Not Available"


def _evidence_classes(row: pd.Series) -> list[str]:
    evidence = []
    if normalize_text(row.get("developability_risk_class")):
        evidence.append("developability")
    if "total_cdr_liabilities" in row.index or "total_fr_liabilities" in row.index:
        evidence.append("CDR/FR")
    if normalize_text(row.get("humanness_risk_class")) not in {"", "Unknown", "Not Available"} or _safe_bool(row.get("humanness_available")):
        evidence.append("humanness")
    if normalize_text(row.get("structural_risk_class")) not in {"", "Not Available"} or _safe_bool(row.get("structure_available")):
        evidence.append("structure")
    if normalize_text(row.get("formulation_risk_class")):
        evidence.append("formulation")
    if _safe_bool(row.get("external_tool_results_available")):
        evidence.append("external tool")
    return evidence


def _confidence_level(row: pd.Series) -> str:
    evidence_count = len(_evidence_classes(row))
    if evidence_count >= 5:
        return "High"
    if evidence_count >= 3:
        return "Medium"
    return "Low"


def _review_flags(row: pd.Series) -> str:
    flags: list[str] = []
    if normalize_text(row.get("developability_risk_class")) == "High":
        flags.append("High sequence developability risk")
    if _safe_int(row.get("total_cdr_liabilities")) > 0:
        flags.append("CDR liability present")
    if normalize_text(row.get("humanness_risk_class")) == "High":
        flags.append("High humanness risk in imported data")
    if normalize_text(row.get("structural_risk_class")) == "High":
        flags.append("High structural risk in imported data")
    if normalize_text(row.get("formulation_risk_class")) == "High":
        flags.append("High formulation risk")
    if _safe_int(row.get("external_high_risk_flags")) > 0:
        flags.append("External high-risk flag")
    if not flags:
        flags.append("No major integrated review flag detected")
    return join_messages(flags)


def _strengths(row: pd.Series) -> str:
    existing = normalize_text(row.get("strengths"))
    if existing:
        return existing
    strengths: list[str] = []
    if normalize_text(row.get("developability_risk_class")) == "Low":
        strengths.append("Low sequence developability risk")
    if _safe_int(row.get("total_cdr_liabilities")) == 0:
        strengths.append("No CDR-mapped liability detected")
    if normalize_text(row.get("humanness_risk_class")) == "Low":
        strengths.append("Low humanness risk in imported data")
    if normalize_text(row.get("structural_risk_class")) == "Low":
        strengths.append("Low structural risk in imported data")
    if normalize_text(row.get("formulation_risk_class")) == "Low":
        strengths.append("Low formulation risk by screening rules")
    return join_messages(strengths) or "No specific strength highlighted by integrated rules"


def _weaknesses(row: pd.Series) -> str:
    existing = normalize_text(row.get("weaknesses"))
    if existing:
        return existing
    flags = _review_flags(row)
    if flags == "No major integrated review flag detected":
        return "No major integrated weakness detected"
    return flags


def _recommended_next_action(row: pd.Series, suggestion: str) -> str:
    existing = normalize_text(row.get("next_step_recommendation")) or normalize_text(row.get("formulation_next_step_recommendation"))
    if suggestion == "Advance":
        return existing or "Proceed to next-stage expression, binding, and developability validation planning"
    if suggestion == "Advance with review":
        return existing or "Proceed with targeted review of flagged computational signals"
    if suggestion == "Engineering review":
        return "Perform integrated engineering review before progression"
    return "Deprioritize, redesign, or hold pending additional review"


def _final_interpretation(row: pd.Series, suggestion: str, confidence_level: str) -> str:
    antibody_id = normalize_text(row.get("antibody_id")) or "This candidate"
    priority_class = normalize_text(row.get("final_priority_class")) or "unclassified"
    flags = _review_flags(row)
    return (
        f"{antibody_id} is assigned priority class {priority_class} with go/no-go suggestion "
        f"'{suggestion}'. Evidence completeness is {confidence_level}; this is a computational triage "
        f"assessment, not an experimental success probability. Key integrated review basis: {flags}."
    )


def build_final_assessment(
    antibody_summary_df,
    candidate_ranking_df,
    formulation_recommendations_df=None,
    structural_risk_summary_df=None,
    external_tool_summary_df=None,
) -> pd.DataFrame:
    """Build one final integrated assessment row per antibody."""
    summary = _df(antibody_summary_df).copy()
    ranking = _df(candidate_ranking_df).copy()
    formulation = _df(formulation_recommendations_df)
    structural = _df(structural_risk_summary_df)
    external = _df(external_tool_summary_df)

    if ranking.empty and summary.empty:
        return pd.DataFrame(columns=FINAL_ASSESSMENT_COLUMNS)
    base = ranking if not ranking.empty else summary
    if "antibody_id" not in base.columns:
        return pd.DataFrame(columns=FINAL_ASSESSMENT_COLUMNS)

    table = base.drop_duplicates(subset=["antibody_id"], keep="first").copy()
    if not summary.empty and "antibody_id" in summary.columns and base is not summary:
        table = _merge_unique(table, summary, "summary")
    table = _merge_unique(table, formulation, "formulation")
    table = _merge_unique(table, structural, "structure")
    table = _merge_unique(table, external, "external")

    records: list[dict[str, object]] = []
    for _, row in table.iterrows():
        priority_class = normalize_text(row.get("final_priority_class")) or "C"
        developability_risk = (
            normalize_text(row.get("cdr_adjusted_max_risk_class"))
            or normalize_text(row.get("max_chain_risk_class"))
            or "Low"
        )
        humanness_risk = normalize_text(row.get("max_humanness_risk_class")) or "Unknown"
        structural_risk = (
            normalize_text(row.get("structural_risk_class"))
            or normalize_text(row.get("structural_risk_class_structure"))
            or "Not Available"
        )
        formulation_risk = (
            normalize_text(row.get("formulation_risk_class"))
            or normalize_text(row.get("formulation_risk_class_formulation"))
            or "Not Available"
        )
        external_high = _safe_int(row.get("external_high_risk_flags", row.get("external_high_risk_flags_external", 0)))
        external_medium = _safe_int(row.get("external_medium_risk_flags", row.get("external_medium_risk_flags_external", 0)))
        external_available = _safe_bool(row.get("external_tool_results_available", row.get("external_tool_results_available_external", False)))

        enriched = row.copy()
        enriched["developability_risk_class"] = developability_risk
        enriched["humanness_risk_class"] = humanness_risk
        enriched["structural_risk_class"] = structural_risk
        enriched["formulation_risk_class"] = formulation_risk
        enriched["external_high_risk_flags"] = external_high
        enriched["external_medium_risk_flags"] = external_medium
        enriched["external_tool_results_available"] = external_available

        suggestion = _base_go_no_go(priority_class)
        if (
            structural_risk == "High"
            or formulation_risk == "High"
            or external_high > 0
            or (_safe_int(enriched.get("high_humanness_risk_chain_count")) > 0 and _safe_int(enriched.get("total_cdr_liabilities")) > 0)
        ):
            suggestion = _at_least_engineering_review(suggestion)
        confidence = _confidence_level(enriched)

        records.append(
            {
                "antibody_id": row.get("antibody_id", ""),
                "molecule_format": row.get("molecule_format", row.get("molecule_format_summary", "")),
                "final_priority_class": priority_class,
                "final_priority_score": row.get("final_priority_score", ""),
                "decision_label": row.get("decision_label", ""),
                "developability_risk_class": developability_risk,
                "humanness_risk_class": humanness_risk,
                "structural_risk_class": structural_risk,
                "formulation_risk_class": formulation_risk,
                "external_tool_risk_flag": _risk_flag(enriched),
                "key_strengths": _strengths(enriched),
                "key_weaknesses": _weaknesses(enriched),
                "major_review_flags": _review_flags(enriched),
                "recommended_next_action": _recommended_next_action(enriched, suggestion),
                "go_no_go_suggestion": suggestion,
                "confidence_level": confidence,
                "final_interpretation": _final_interpretation(enriched, suggestion, confidence),
            }
        )

    result = pd.DataFrame(records, columns=FINAL_ASSESSMENT_COLUMNS)
    if not result.empty:
        result = result.sort_values(
            ["final_priority_score", "antibody_id"],
            ascending=[False, True],
            kind="mergesort",
        ).reset_index(drop=True)
    return result


def build_executive_decision_summary(final_assessment_df) -> dict:
    """Return summary metrics for Streamlit and HTML reporting."""
    final_assessment = _df(final_assessment_df)
    if final_assessment.empty:
        return {
            "total_candidates": 0,
            "advance_count": 0,
            "advance_with_review_count": 0,
            "engineering_review_count": 0,
            "deprioritize_count": 0,
            "top_ranked_candidate": "",
            "high_risk_candidate_count": 0,
            "candidates_with_high_confidence": 0,
            "candidates_missing_optional_data": 0,
            "overall_project_summary": (
                "The current batch contains 0 candidates. 0 candidates are recommended for advancement or "
                "advancement with review. 0 candidates require engineering review or redesign based on "
                "integrated computational assessment."
            ),
        }

    suggestion_counts = final_assessment.get("go_no_go_suggestion", pd.Series(dtype=str)).astype(str).value_counts()
    top_ranked_candidate = ""
    if "final_priority_score" in final_assessment.columns:
        ranked = final_assessment.copy()
        ranked["_score"] = pd.to_numeric(ranked["final_priority_score"], errors="coerce").fillna(-1)
        ranked = ranked.sort_values(["_score", "antibody_id"], ascending=[False, True], kind="mergesort")
        if not ranked.empty:
            top_ranked_candidate = normalize_text(ranked.iloc[0].get("antibody_id"))
    elif "antibody_id" in final_assessment.columns and not final_assessment.empty:
        top_ranked_candidate = normalize_text(final_assessment.iloc[0].get("antibody_id"))

    high_risk_candidate_count = 0
    for column in ["developability_risk_class", "humanness_risk_class", "structural_risk_class", "formulation_risk_class", "external_tool_risk_flag"]:
        if column in final_assessment.columns:
            high_risk_candidate_count += final_assessment[column].astype(str).eq("High").astype(int)
    if not isinstance(high_risk_candidate_count, int):
        high_risk_candidate_count = int((high_risk_candidate_count > 0).sum())

    total = int(len(final_assessment))
    advance = int(suggestion_counts.get("Advance", 0))
    advance_review = int(suggestion_counts.get("Advance with review", 0))
    engineering = int(suggestion_counts.get("Engineering review", 0))
    deprioritize = int(suggestion_counts.get("Deprioritize / redesign", 0))
    missing_optional = (
        int(final_assessment["confidence_level"].astype(str).isin(["Low", "Medium"]).sum())
        if "confidence_level" in final_assessment.columns
        else 0
    )
    review_or_redesign = engineering + deprioritize

    return {
        "total_candidates": total,
        "advance_count": advance,
        "advance_with_review_count": advance_review,
        "engineering_review_count": engineering,
        "deprioritize_count": deprioritize,
        "top_ranked_candidate": top_ranked_candidate,
        "high_risk_candidate_count": high_risk_candidate_count,
        "candidates_with_high_confidence": (
            int(final_assessment["confidence_level"].astype(str).eq("High").sum())
            if "confidence_level" in final_assessment.columns
            else 0
        ),
        "candidates_missing_optional_data": missing_optional,
        "overall_project_summary": (
            f"The current batch contains {total} candidates. {advance + advance_review} candidates are "
            f"recommended for advancement or advancement with review. {review_or_redesign} candidates require "
            "engineering review or redesign based on integrated computational assessment."
        ),
    }
