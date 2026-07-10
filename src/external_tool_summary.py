"""Antibody-level summaries for imported external tool results."""

from __future__ import annotations

import pandas as pd


SUMMARY_COLUMNS = [
    "antibody_id",
    "external_tool_results_available",
    "tools_with_results",
    "humanness_tool_result_available",
    "developability_tool_result_available",
    "structure_tool_result_available",
    "external_high_risk_flags",
    "external_medium_risk_flags",
    "external_low_risk_flags",
    "external_tool_summary_text",
    "external_tool_review_recommendation",
]


def empty_external_tool_summary(antibody_ids: list[object] | None = None) -> pd.DataFrame:
    records = []
    for antibody_id in antibody_ids or []:
        records.append(
            {
                "antibody_id": antibody_id,
                "external_tool_results_available": False,
                "tools_with_results": "",
                "humanness_tool_result_available": False,
                "developability_tool_result_available": False,
                "structure_tool_result_available": False,
                "external_high_risk_flags": 0,
                "external_medium_risk_flags": 0,
                "external_low_risk_flags": 0,
                "external_tool_summary_text": "No external tool result was imported.",
                "external_tool_review_recommendation": "No external tool result was imported.",
            }
        )
    return pd.DataFrame(records, columns=SUMMARY_COLUMNS)


def build_external_tool_summary(external_results_df) -> pd.DataFrame:
    """Build one summary row per antibody from standardized external results."""
    if external_results_df is None or external_results_df.empty or "antibody_id" not in external_results_df.columns:
        return empty_external_tool_summary()

    results = external_results_df.copy()
    for column in ["tool_name", "tool_category", "result_risk_class"]:
        if column not in results.columns:
            results[column] = ""
    records: list[dict[str, object]] = []
    for antibody_id, rows in results.groupby("antibody_id", dropna=False):
        risk_classes = rows["result_risk_class"].astype(str).str.strip()
        tools = sorted(set(rows["tool_name"].replace("", pd.NA).dropna().astype(str).tolist()))
        categories = set(rows["tool_category"].astype(str).str.lower().tolist())
        high_count = int(risk_classes.eq("High").sum())
        medium_count = int(risk_classes.eq("Medium").sum())
        low_count = int(risk_classes.eq("Low").sum())
        if high_count > 0:
            recommendation = "External tool results indicate high-priority review."
        elif low_count > 0 and high_count == 0 and medium_count == 0:
            recommendation = "External tool results do not indicate major additional risk in imported results."
        elif medium_count > 0:
            recommendation = "External tool results indicate medium-priority review."
        else:
            recommendation = "External tool results were imported; review source metrics before decision-making."
        records.append(
            {
                "antibody_id": antibody_id,
                "external_tool_results_available": True,
                "tools_with_results": "; ".join(tools),
                "humanness_tool_result_available": "humanness" in categories,
                "developability_tool_result_available": "developability" in categories,
                "structure_tool_result_available": "structure" in categories,
                "external_high_risk_flags": high_count,
                "external_medium_risk_flags": medium_count,
                "external_low_risk_flags": low_count,
                "external_tool_summary_text": (
                    f"Imported external results from {len(tools)} tool(s): {', '.join(tools) if tools else 'unknown tool'}; "
                    f"High={high_count}, Medium={medium_count}, Low={low_count}."
                ),
                "external_tool_review_recommendation": recommendation,
            }
        )
    return pd.DataFrame(records, columns=SUMMARY_COLUMNS)
