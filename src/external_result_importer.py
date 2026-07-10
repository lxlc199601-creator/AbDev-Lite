"""Import and merge user-provided external tool result files."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd

from .external_tool_registry import load_tool_registry
from .external_tool_summary import build_external_tool_summary, empty_external_tool_summary
from .utils import join_messages, normalize_text


RESULT_COLUMNS = [
    "tool_id",
    "tool_name",
    "tool_category",
    "antibody_id",
    "chain_id",
    "sequence_scope",
    "result_status",
    "result_metric_name",
    "result_metric_value",
    "result_risk_class",
    "result_interpretation",
    "result_file_name",
    "imported_at",
    "merge_status",
    "merge_warning",
]


def empty_external_results() -> pd.DataFrame:
    return pd.DataFrame(columns=RESULT_COLUMNS)


def load_external_tool_result(uploaded_file, tool_id: str) -> pd.DataFrame:
    """Load a user-provided CSV or XLSX result file."""
    if uploaded_file is None:
        return empty_external_results()
    filename = getattr(uploaded_file, "name", str(uploaded_file))
    suffix = Path(filename).suffix.lower()
    if hasattr(uploaded_file, "seek"):
        uploaded_file.seek(0)
    if suffix in {".xlsx", ".xls"}:
        result = pd.read_excel(uploaded_file)
    elif suffix == ".csv":
        result = pd.read_csv(uploaded_file)
    else:
        raise ValueError("External tool result file must be CSV or XLSX.")
    result.attrs["result_file_name"] = Path(filename).name
    result.attrs["tool_id"] = tool_id
    return result


def _tool_lookup(tool_id: str) -> dict[str, object]:
    registry = load_tool_registry(Path("tool_specs") / "external_tools.yaml")
    if registry.empty or "tool_id" not in registry.columns:
        return {"tool_id": tool_id, "tool_name": tool_id, "tool_category": "custom"}
    rows = registry[registry["tool_id"].astype(str).eq(str(tool_id))]
    if rows.empty:
        return {"tool_id": tool_id, "tool_name": tool_id, "tool_category": "custom"}
    return rows.iloc[0].to_dict()


def standardize_external_results(result_df, tool_id: str) -> pd.DataFrame:
    """Normalize arbitrary imported external result tables to the v0.9 schema."""
    if result_df is None or result_df.empty:
        return empty_external_results()
    source = result_df.copy()
    tool = _tool_lookup(tool_id)
    file_name = source.attrs.get("result_file_name", "")
    imported_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    aliases = {
        "result_risk_class": ["result_risk_class", "risk_class", "risk", "external_risk_class"],
        "result_metric_name": ["result_metric_name", "metric_name", "metric", "score_name"],
        "result_metric_value": ["result_metric_value", "metric_value", "score", "value"],
        "result_interpretation": ["result_interpretation", "interpretation", "comment", "notes"],
        "result_status": ["result_status", "status"],
    }
    records: list[dict[str, object]] = []
    for _, row in source.iterrows():
        warnings: list[str] = []

        def pick(column: str, default: str = "") -> object:
            if column in source.columns:
                return row.get(column, default)
            for alias in aliases.get(column, []):
                if alias in source.columns:
                    return row.get(alias, default)
            warnings.append(f"Missing {column}")
            return default

        antibody_id = row.get("antibody_id", "")
        chain_id = row.get("chain_id", "")
        sequence_scope = row.get("sequence_scope", "")
        if not normalize_text(antibody_id):
            warnings.append("Missing antibody_id")
        risk_class = normalize_text(pick("result_risk_class"))
        if risk_class not in {"Low", "Medium", "High", "Unknown", ""}:
            risk_class = risk_class.capitalize()
        records.append(
            {
                "tool_id": row.get("tool_id", tool_id) or tool_id,
                "tool_name": row.get("tool_name", tool.get("tool_name", tool_id)) or tool.get("tool_name", tool_id),
                "tool_category": row.get("tool_category", tool.get("tool_category", "custom")) or tool.get("tool_category", "custom"),
                "antibody_id": antibody_id,
                "chain_id": chain_id,
                "sequence_scope": sequence_scope,
                "result_status": pick("result_status", "IMPORTED") or "IMPORTED",
                "result_metric_name": pick("result_metric_name"),
                "result_metric_value": pick("result_metric_value"),
                "result_risk_class": risk_class or "Unknown",
                "result_interpretation": pick("result_interpretation"),
                "result_file_name": row.get("result_file_name", file_name) or file_name,
                "imported_at": row.get("imported_at", imported_at) or imported_at,
                "merge_status": "pending_merge",
                "merge_warning": join_messages(warnings),
            }
        )
    return pd.DataFrame(records, columns=RESULT_COLUMNS)


def _merge_summary_fields(antibody_summary_df: pd.DataFrame, external_summary_df: pd.DataFrame) -> pd.DataFrame:
    summary = antibody_summary_df.copy() if isinstance(antibody_summary_df, pd.DataFrame) else pd.DataFrame()
    if summary.empty:
        return summary
    if "antibody_id" not in summary.columns:
        return summary
    default_summary = empty_external_tool_summary(summary["antibody_id"].dropna().tolist())
    if external_summary_df is not None and not external_summary_df.empty:
        default_summary = default_summary.drop(columns=[column for column in external_summary_df.columns if column != "antibody_id"], errors="ignore")
        default_summary = default_summary.merge(external_summary_df, on="antibody_id", how="left", suffixes=("", "_imported"))
        for column in default_summary.columns:
            if column.endswith("_imported"):
                base = column.removesuffix("_imported")
                default_summary[base] = default_summary[column].combine_first(default_summary.get(base))
        default_summary = default_summary[[column for column in empty_external_tool_summary().columns if column in default_summary.columns]]
    summary = summary.drop(columns=[column for column in default_summary.columns if column != "antibody_id" and column in summary.columns])
    summary = summary.merge(default_summary, on="antibody_id", how="left")
    summary["external_tool_results_available"] = summary["external_tool_results_available"].fillna(False).astype(bool)
    for column in ["external_high_risk_flags", "external_medium_risk_flags", "external_low_risk_flags"]:
        summary[column] = summary[column].fillna(0).astype(int)
    summary["tools_with_results"] = summary["tools_with_results"].fillna("")
    summary["external_tool_summary_text"] = summary["external_tool_summary_text"].fillna("No external tool result was imported.")
    summary["external_tool_review_recommendation"] = summary["external_tool_review_recommendation"].fillna("No external tool result was imported.")
    if "strengths" in summary.columns:
        mask = summary["external_tool_results_available"] & summary["external_high_risk_flags"].eq(0) & summary["external_medium_risk_flags"].eq(0)
        summary.loc[mask, "strengths"] = summary.loc[mask, "strengths"].apply(
            lambda value: join_messages([value, "Imported external tool results did not add major risk flags in this MVP assessment"])
        )
    if "weaknesses" in summary.columns:
        mask = summary["external_high_risk_flags"].gt(0)
        summary.loc[mask, "weaknesses"] = summary.loc[mask, "weaknesses"].apply(
            lambda value: join_messages([value, "Imported external tool results include high-risk flags requiring review"])
        )
    return summary


def _merge_ranking_fields(candidate_ranking_df: pd.DataFrame, external_summary_df: pd.DataFrame) -> pd.DataFrame:
    ranking = candidate_ranking_df.copy() if isinstance(candidate_ranking_df, pd.DataFrame) else pd.DataFrame()
    if ranking.empty or "antibody_id" not in ranking.columns:
        return ranking
    external = empty_external_tool_summary(ranking["antibody_id"].dropna().tolist())
    if external_summary_df is not None and not external_summary_df.empty:
        external = external.drop(columns=[column for column in external_summary_df.columns if column != "antibody_id"], errors="ignore")
        external = external.merge(external_summary_df, on="antibody_id", how="left")
    ranking = ranking.drop(columns=["external_high_risk_flags", "external_medium_risk_flags", "external_tool_results_available"], errors="ignore")
    ranking = ranking.merge(
        external[["antibody_id", "external_high_risk_flags", "external_medium_risk_flags", "external_tool_results_available"]],
        on="antibody_id",
        how="left",
    )
    ranking["external_high_risk_flags"] = ranking["external_high_risk_flags"].fillna(0).astype(int)
    ranking["external_medium_risk_flags"] = ranking["external_medium_risk_flags"].fillna(0).astype(int)
    ranking["external_tool_results_available"] = ranking["external_tool_results_available"].fillna(False).astype(bool)
    if "final_priority_score" in ranking.columns:
        ranking["final_priority_score"] = (
            pd.to_numeric(ranking["final_priority_score"], errors="coerce").fillna(0)
            - ranking["external_high_risk_flags"] * 10
            - ranking["external_medium_risk_flags"] * 5
        ).clip(lower=0, upper=100).round(2)
    if "review_reason" in ranking.columns:
        high_mask = ranking["external_high_risk_flags"].gt(0)
        no_result_mask = ~ranking["external_tool_results_available"]
        ranking.loc[high_mask, "review_reason"] = ranking.loc[high_mask, "review_reason"].apply(
            lambda value: join_messages([value, "High-risk flags imported from external tool results"])
        )
        ranking.loc[no_result_mask, "review_reason"] = ranking.loc[no_result_mask, "review_reason"].apply(
            lambda value: join_messages([value, "No external tool result imported"])
        )
    return ranking


def merge_external_results(
    antibody_summary_df,
    candidate_ranking_df,
    external_results_df,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge imported external result summary fields into antibody summary and ranking."""
    external_results = external_results_df if isinstance(external_results_df, pd.DataFrame) else empty_external_results()
    external_summary = build_external_tool_summary(external_results)
    summary = _merge_summary_fields(antibody_summary_df, external_summary)
    ranking = _merge_ranking_fields(candidate_ranking_df, external_summary)
    if not external_results.empty and "antibody_id" in external_results.columns:
        known_ids = set(summary["antibody_id"].astype(str)) if not summary.empty and "antibody_id" in summary.columns else set()
        external_results["merge_status"] = external_results["antibody_id"].astype(str).apply(
            lambda value: "matched_antibody" if value in known_ids else "unmatched_antibody"
        )
        external_results["merge_warning"] = external_results.apply(
            lambda row: join_messages(
                [
                    row.get("merge_warning", ""),
                    "No matching antibody_id in AbDev-Lite summary" if row.get("merge_status") == "unmatched_antibody" else "",
                ]
            ),
            axis=1,
        )
    return summary, ranking
