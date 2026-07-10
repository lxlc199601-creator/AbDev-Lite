"""Optional external structure prediction result import utilities."""

from __future__ import annotations

from pathlib import Path
from typing import BinaryIO

import pandas as pd

from .structural_risk import STRUCTURAL_RISK_SUMMARY_COLUMNS, build_structural_risk_summary
from .utils import join_messages, normalize_text


STRUCTURE_INPUT_COLUMNS = [
    "antibody_id",
    "chain_id",
    "sequence_scope",
    "structure_tool",
    "structure_model_file",
    "structure_status",
    "mean_plddt",
    "cdr1_plddt",
    "cdr2_plddt",
    "cdr3_plddt",
    "vh_vl_orientation_confidence",
    "predicted_surface_hydrophobic_patch_score",
    "predicted_aggregation_patch_score",
    "predicted_charge_patch_score",
    "structural_notes",
]

STRUCTURE_RESULT_COLUMNS = [
    *STRUCTURE_INPUT_COLUMNS,
    "merge_status",
    "merge_warning",
]

NUMERIC_COLUMNS = [
    "mean_plddt",
    "cdr1_plddt",
    "cdr2_plddt",
    "cdr3_plddt",
    "vh_vl_orientation_confidence",
    "predicted_surface_hydrophobic_patch_score",
    "predicted_aggregation_patch_score",
    "predicted_charge_patch_score",
]

STRUCTURAL_SUMMARY_COLUMNS = [
    "structure_available",
    "structural_risk_class",
    "structural_risk_score",
    "structural_review_reason",
    "structural_next_step_recommendation",
]


def _empty_results_df() -> pd.DataFrame:
    return pd.DataFrame(columns=STRUCTURE_RESULT_COLUMNS)


def _read_uploaded_table(uploaded_file: BinaryIO | str | Path) -> pd.DataFrame:
    name = getattr(uploaded_file, "name", "") or str(uploaded_file)
    suffix = Path(name).suffix.lower()
    try:
        if suffix == ".csv":
            return pd.read_csv(uploaded_file)
        if suffix in {".xlsx", ".xls"}:
            return pd.read_excel(uploaded_file)
    except pd.errors.EmptyDataError as exc:
        raise ValueError("The structure result file is empty.") from exc
    except Exception as exc:
        raise ValueError(f"Could not read the structure result file: {exc}") from exc
    raise ValueError("Unsupported structure result file format. Please upload CSV or XLSX.")


def load_structure_results(uploaded_file) -> pd.DataFrame:
    """Read an optional external structure prediction summary CSV or XLSX file."""
    if uploaded_file is None:
        return _empty_results_df()
    return validate_structure_results_df(_read_uploaded_table(uploaded_file))


def validate_structure_results_df(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize imported structure result columns."""
    if df is None or df.empty:
        return _empty_results_df()

    table = df.dropna(how="all").copy()
    if table.empty:
        return _empty_results_df()
    table.columns = [normalize_text(column).strip().lower().replace(" ", "_") for column in table.columns]
    if "antibody_id" not in table.columns:
        raise ValueError("Missing required structure result column: antibody_id.")

    for column in STRUCTURE_INPUT_COLUMNS:
        if column not in table.columns:
            table[column] = pd.NA
    table = table[STRUCTURE_INPUT_COLUMNS].copy()

    for column in [
        "antibody_id",
        "chain_id",
        "sequence_scope",
        "structure_tool",
        "structure_model_file",
        "structure_status",
        "structural_notes",
    ]:
        table[column] = table[column].map(normalize_text)
    table["structure_status"] = table["structure_status"].str.upper().replace("", "NOT_RUN")
    for column in NUMERIC_COLUMNS:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table["merge_status"] = "unmerged"
    table["merge_warning"] = ""
    return table[STRUCTURE_RESULT_COLUMNS]


def _build_chain_lookup(antibody_summary_df: pd.DataFrame, candidate_ranking_df: pd.DataFrame) -> tuple[set[str], set[tuple[str, str]], set[tuple[str, str]]]:
    antibody_ids: set[str] = set()
    by_chain: set[tuple[str, str]] = set()
    by_scope: set[tuple[str, str]] = set()
    for frame in [antibody_summary_df, candidate_ranking_df]:
        if frame is None or frame.empty:
            continue
        if "antibody_id" in frame.columns:
            antibody_ids.update(frame["antibody_id"].map(normalize_text).tolist())
        for _, row in frame.iterrows():
            antibody_id = normalize_text(row.get("antibody_id"))
            chain_id = normalize_text(row.get("chain_id")).lower()
            sequence_scope = normalize_text(row.get("sequence_scope")).lower()
            if antibody_id and chain_id:
                by_chain.add((antibody_id, chain_id))
            if antibody_id and sequence_scope:
                by_scope.add((antibody_id, sequence_scope))
            scopes = [item.strip().lower() for item in normalize_text(row.get("sequence_scopes")).split(",") if item.strip()]
            for scope in scopes:
                by_scope.add((antibody_id, scope))
    return antibody_ids, by_chain, by_scope


def _merge_status(row: pd.Series, antibody_ids: set[str], by_chain: set[tuple[str, str]], by_scope: set[tuple[str, str]]) -> tuple[str, str]:
    antibody_id = normalize_text(row.get("antibody_id"))
    chain_id = normalize_text(row.get("chain_id")).lower()
    sequence_scope = normalize_text(row.get("sequence_scope")).lower()
    if not antibody_id:
        return "unmatched", "Missing antibody_id in imported structure row"
    if chain_id and (antibody_id, chain_id) in by_chain:
        return "matched_by_antibody_id_chain_id", ""
    if sequence_scope and (antibody_id, sequence_scope) in by_scope:
        return "matched_by_antibody_id_sequence_scope", ""
    if antibody_id in antibody_ids:
        return "matched_by_antibody_id", ""
    warnings: list[str] = ["No matching analyzed antibody or chain found for imported structure row"]
    if not chain_id and not sequence_scope:
        warnings.append("Missing both chain_id and sequence_scope")
    return "unmatched", join_messages(warnings)


def _add_default_structural_columns(df: pd.DataFrame) -> pd.DataFrame:
    table = df.copy()
    defaults = {
        "structure_available": False,
        "structural_risk_class": "Not Available",
        "structural_risk_score": 0,
        "structural_review_reason": "No external structure prediction result was provided",
        "structural_next_step_recommendation": "Import external structure prediction summary metrics if structural triage is needed",
    }
    for column, default in defaults.items():
        if column not in table.columns:
            table[column] = default
        else:
            table[column] = table[column].fillna(default).replace("", default)
    return table


def _append_message(existing: object, message: str) -> str:
    return join_messages([normalize_text(existing), message])


def _merge_structural_summary(df: pd.DataFrame, structural_summary: pd.DataFrame) -> pd.DataFrame:
    table = _add_default_structural_columns(df if df is not None else pd.DataFrame())
    if table.empty or structural_summary is None or structural_summary.empty or "antibody_id" not in table.columns:
        return table
    merge_columns = ["antibody_id", *STRUCTURAL_SUMMARY_COLUMNS]
    for column in merge_columns:
        if column not in structural_summary.columns:
            structural_summary[column] = pd.NA
    table = table.drop(columns=[column for column in STRUCTURAL_SUMMARY_COLUMNS if column in table.columns])
    table = table.merge(structural_summary[merge_columns], on="antibody_id", how="left")
    table = _add_default_structural_columns(table)
    if "strengths" in table.columns:
        low_mask = table["structural_risk_class"].astype(str).eq("Low")
        table.loc[low_mask, "strengths"] = table.loc[low_mask, "strengths"].map(
            lambda value: _append_message(
                value,
                "Imported structure metrics suggest low structural risk in this MVP assessment",
            )
        )
    if "weaknesses" in table.columns:
        high_mask = table["structural_risk_class"].astype(str).eq("High")
        table.loc[high_mask, "weaknesses"] = table.loc[high_mask, "weaknesses"].map(
            lambda value: _append_message(
                value,
                "Imported structure metrics suggest high structural risk; structural review is recommended",
            )
        )
    return table


def annotate_structure_result_matches(
    antibody_summary_df: pd.DataFrame,
    candidate_ranking_df: pd.DataFrame,
    structure_results_df: pd.DataFrame,
) -> pd.DataFrame:
    """Add merge status and warnings to imported structure result rows."""
    if structure_results_df is None or structure_results_df.empty:
        return _empty_results_df()
    structure_results = validate_structure_results_df(structure_results_df)
    antibody_ids, by_chain, by_scope = _build_chain_lookup(antibody_summary_df, candidate_ranking_df)
    statuses = structure_results.apply(lambda row: _merge_status(row, antibody_ids, by_chain, by_scope), axis=1)
    structure_results["merge_status"] = statuses.map(lambda item: item[0])
    structure_results["merge_warning"] = statuses.map(lambda item: item[1])
    return structure_results[STRUCTURE_RESULT_COLUMNS]


def empty_structural_risk_summary_for_antibodies(antibody_summary_df: pd.DataFrame) -> pd.DataFrame:
    """Create Not Available structural rows for analyzed antibodies when no structure file is imported."""
    if antibody_summary_df is None or antibody_summary_df.empty or "antibody_id" not in antibody_summary_df.columns:
        return pd.DataFrame(columns=STRUCTURAL_RISK_SUMMARY_COLUMNS)
    records = []
    for antibody_id in antibody_summary_df["antibody_id"].map(normalize_text).drop_duplicates():
        records.append(
            {
                "antibody_id": antibody_id,
                "structure_available": False,
                "structure_tools": "",
                "structure_model_files": "",
                "structure_status_summary": "Not Available",
                "mean_plddt_min": pd.NA,
                "mean_plddt_mean": pd.NA,
                "cdr1_plddt_min": pd.NA,
                "cdr2_plddt_min": pd.NA,
                "cdr3_plddt_min": pd.NA,
                "low_confidence_cdr_count": 0,
                "high_hydrophobic_patch_flag": False,
                "high_aggregation_patch_flag": False,
                "high_charge_patch_flag": False,
                "structural_risk_score": 0,
                "structural_risk_class": "Not Available",
                "structural_review_reason": "No external structure prediction result was provided",
                "structural_next_step_recommendation": "Import external structure prediction summary metrics if structural triage is needed",
            }
        )
    return pd.DataFrame(records, columns=STRUCTURAL_RISK_SUMMARY_COLUMNS)


def empty_structure_results_for_antibodies(antibody_summary_df: pd.DataFrame) -> pd.DataFrame:
    """Create NOT_RUN structure result rows for analyzed antibodies when no structure file is imported."""
    if antibody_summary_df is None or antibody_summary_df.empty or "antibody_id" not in antibody_summary_df.columns:
        return _empty_results_df()
    records = []
    for antibody_id in antibody_summary_df["antibody_id"].map(normalize_text).drop_duplicates():
        record = {column: pd.NA for column in STRUCTURE_RESULT_COLUMNS}
        record.update(
            {
                "antibody_id": antibody_id,
                "structure_status": "NOT_RUN",
                "merge_status": "not_available",
                "merge_warning": "No external structure prediction result was provided",
            }
        )
        records.append(record)
    return pd.DataFrame(records, columns=STRUCTURE_RESULT_COLUMNS)


def complete_structural_risk_summary_for_antibodies(
    antibody_summary_df: pd.DataFrame,
    structural_risk_summary_df: pd.DataFrame,
) -> pd.DataFrame:
    """Ensure each analyzed antibody has a structural summary row."""
    defaults = empty_structural_risk_summary_for_antibodies(antibody_summary_df)
    if structural_risk_summary_df is None or structural_risk_summary_df.empty:
        return defaults
    if defaults.empty:
        return structural_risk_summary_df[STRUCTURAL_RISK_SUMMARY_COLUMNS]
    present = set(structural_risk_summary_df["antibody_id"].map(normalize_text).tolist())
    missing_defaults = defaults[~defaults["antibody_id"].map(normalize_text).isin(present)]
    return pd.concat(
        [structural_risk_summary_df[STRUCTURAL_RISK_SUMMARY_COLUMNS], missing_defaults],
        ignore_index=True,
    )


def merge_structure_results(
    antibody_summary_df: pd.DataFrame,
    candidate_ranking_df: pd.DataFrame,
    structure_results_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge optional imported structure metrics into antibody summary and candidate ranking tables."""
    summary = _add_default_structural_columns(antibody_summary_df if antibody_summary_df is not None else pd.DataFrame())
    ranking = _add_default_structural_columns(candidate_ranking_df if candidate_ranking_df is not None else pd.DataFrame())
    if structure_results_df is None or structure_results_df.empty:
        return summary, ranking

    structure_results = annotate_structure_result_matches(summary, ranking, structure_results_df)
    structural_summary = build_structural_risk_summary(structure_results)
    return _merge_structural_summary(summary, structural_summary), _merge_structural_summary(ranking, structural_summary)
