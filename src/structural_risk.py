"""Structural risk summary from imported external structure metrics."""

from __future__ import annotations

import pandas as pd

from .utils import join_messages, normalize_text


STRUCTURAL_RISK_SUMMARY_COLUMNS = [
    "antibody_id",
    "structure_available",
    "structure_tools",
    "structure_model_files",
    "structure_status_summary",
    "mean_plddt_min",
    "mean_plddt_mean",
    "cdr1_plddt_min",
    "cdr2_plddt_min",
    "cdr3_plddt_min",
    "low_confidence_cdr_count",
    "high_hydrophobic_patch_flag",
    "high_aggregation_patch_flag",
    "high_charge_patch_flag",
    "structural_risk_score",
    "structural_risk_class",
    "structural_review_reason",
    "structural_next_step_recommendation",
]


def _numeric_min(rows: pd.DataFrame, column: str) -> float | pd.NA:
    if column not in rows.columns:
        return pd.NA
    values = pd.to_numeric(rows[column], errors="coerce").dropna()
    return round(float(values.min()), 2) if not values.empty else pd.NA


def _numeric_mean(rows: pd.DataFrame, column: str) -> float | pd.NA:
    if column not in rows.columns:
        return pd.NA
    values = pd.to_numeric(rows[column], errors="coerce").dropna()
    return round(float(values.mean()), 2) if not values.empty else pd.NA


def _any_below(rows: pd.DataFrame, column: str, cutoff: float) -> bool:
    if column not in rows.columns:
        return False
    values = pd.to_numeric(rows[column], errors="coerce").dropna()
    return bool((values < cutoff).any())


def _any_at_least(rows: pd.DataFrame, column: str, cutoff: float) -> bool:
    if column not in rows.columns:
        return False
    values = pd.to_numeric(rows[column], errors="coerce").dropna()
    return bool((values >= cutoff).any())


def _status_penalty(statuses: pd.Series) -> int:
    clean = statuses.fillna("").astype(str).str.upper()
    score = 0
    if clean.eq("FAIL").any():
        score += 3
    if clean.eq("WARNING").any():
        score += 1
    return score


def _risk_class(score: int, available: bool) -> str:
    if not available:
        return "Not Available"
    if score <= 2:
        return "Low"
    if score <= 5:
        return "Medium"
    return "High"


def _next_step(risk_class: str) -> str:
    return {
        "Low": "Imported structure metrics do not indicate a major structural review flag in this MVP assessment",
        "Medium": "Review imported structure confidence and patch annotations before progression",
        "High": "Structural review is recommended before candidate progression",
        "Not Available": "Import external structure prediction summary metrics if structural triage is needed",
    }.get(risk_class, "Review imported structure metrics before candidate progression")


def _review_reason(
    available: bool,
    low_confidence_cdr_count: int,
    cdr3_low: bool,
    hydrophobic_flag: bool,
    aggregation_flag: bool,
    charge_flag: bool,
) -> str:
    if not available:
        return "No external structure prediction result was provided"
    reasons: list[str] = []
    if low_confidence_cdr_count > 0:
        reasons.append("Low-confidence CDR loop prediction detected")
    if cdr3_low:
        reasons.append("Low-confidence CDR3 prediction detected")
    if hydrophobic_flag:
        reasons.append("Predicted hydrophobic surface patch risk detected")
    if aggregation_flag:
        reasons.append("Predicted aggregation-prone patch risk detected")
    if charge_flag:
        reasons.append("Predicted charge patch risk detected")
    return join_messages(reasons) or "No major structural risk flag detected from imported metrics"


def _score_rows(rows: pd.DataFrame) -> int:
    score = 0
    if _any_below(rows, "mean_plddt", 70):
        score += 2
    if _any_below(rows, "cdr1_plddt", 70):
        score += 1
    if _any_below(rows, "cdr2_plddt", 70):
        score += 1
    if _any_below(rows, "cdr3_plddt", 70):
        score += 2
    if _any_below(rows, "vh_vl_orientation_confidence", 0.6):
        score += 2
    if _any_at_least(rows, "predicted_surface_hydrophobic_patch_score", 0.7):
        score += 2
    if _any_at_least(rows, "predicted_aggregation_patch_score", 0.7):
        score += 2
    if _any_at_least(rows, "predicted_charge_patch_score", 0.7):
        score += 1
    if "structure_status" in rows.columns:
        score += _status_penalty(rows["structure_status"])
    return int(score)


def build_structural_risk_summary(structure_results_df: pd.DataFrame) -> pd.DataFrame:
    """Build antibody-level structural risk summary from imported result rows."""
    if structure_results_df is None or structure_results_df.empty or "antibody_id" not in structure_results_df.columns:
        return pd.DataFrame(columns=STRUCTURAL_RISK_SUMMARY_COLUMNS)

    table = structure_results_df.copy()
    records: list[dict[str, object]] = []
    if "merge_status" in table.columns:
        statuses = table["merge_status"].fillna("").astype(str)
        if statuses.str.startswith("matched").any():
            matched = table[statuses.str.startswith("matched")]
        elif statuses.eq("unmerged").all():
            matched = table
        else:
            matched = pd.DataFrame()
    else:
        matched = table
    if matched.empty:
        return pd.DataFrame(columns=STRUCTURAL_RISK_SUMMARY_COLUMNS)

    for antibody_id, rows in matched.groupby("antibody_id", dropna=False):
        tools = join_messages(rows.get("structure_tool", pd.Series(dtype=str)).map(normalize_text).tolist())
        model_files = join_messages(rows.get("structure_model_file", pd.Series(dtype=str)).map(normalize_text).tolist())
        statuses = rows.get("structure_status", pd.Series(dtype=str)).map(lambda value: normalize_text(value).upper())
        status_summary = join_messages(f"{status} ({count})" for status, count in statuses.value_counts().items() if status)
        cdr_low_flags = [
            _any_below(rows, "cdr1_plddt", 70),
            _any_below(rows, "cdr2_plddt", 70),
            _any_below(rows, "cdr3_plddt", 70),
        ]
        low_confidence_cdr_count = int(sum(cdr_low_flags))
        hydrophobic_flag = _any_at_least(rows, "predicted_surface_hydrophobic_patch_score", 0.7)
        aggregation_flag = _any_at_least(rows, "predicted_aggregation_patch_score", 0.7)
        charge_flag = _any_at_least(rows, "predicted_charge_patch_score", 0.7)
        score = _score_rows(rows)
        risk_class = _risk_class(score, True)

        records.append(
            {
                "antibody_id": antibody_id,
                "structure_available": True,
                "structure_tools": tools,
                "structure_model_files": model_files,
                "structure_status_summary": status_summary,
                "mean_plddt_min": _numeric_min(rows, "mean_plddt"),
                "mean_plddt_mean": _numeric_mean(rows, "mean_plddt"),
                "cdr1_plddt_min": _numeric_min(rows, "cdr1_plddt"),
                "cdr2_plddt_min": _numeric_min(rows, "cdr2_plddt"),
                "cdr3_plddt_min": _numeric_min(rows, "cdr3_plddt"),
                "low_confidence_cdr_count": low_confidence_cdr_count,
                "high_hydrophobic_patch_flag": hydrophobic_flag,
                "high_aggregation_patch_flag": aggregation_flag,
                "high_charge_patch_flag": charge_flag,
                "structural_risk_score": score,
                "structural_risk_class": risk_class,
                "structural_review_reason": _review_reason(
                    True,
                    low_confidence_cdr_count,
                    cdr_low_flags[2],
                    hydrophobic_flag,
                    aggregation_flag,
                    charge_flag,
                ),
                "structural_next_step_recommendation": _next_step(risk_class),
            }
        )
    return pd.DataFrame(records, columns=STRUCTURAL_RISK_SUMMARY_COLUMNS)
