"""Optional humanness and germline assessment import utilities."""

from __future__ import annotations

from pathlib import Path
from typing import BinaryIO

import pandas as pd

from .utils import join_messages, normalize_text


HUMANNESS_INPUT_COLUMNS = [
    "antibody_id",
    "chain_id",
    "sequence_scope",
    "humanness_tool",
    "humanness_score",
    "humanness_percentile",
    "closest_human_germline",
    "closest_species",
    "v_gene",
    "j_gene",
    "identity_to_human_germline",
    "framework_identity",
    "cdr_identity",
    "humanness_notes",
]

HUMANNESS_RESULT_COLUMNS = [
    "antibody_id",
    "chain_id",
    "sequence_scope",
    "humanness_tool",
    "humanness_score",
    "humanness_percentile",
    "closest_human_germline",
    "closest_species",
    "v_gene",
    "j_gene",
    "identity_to_human_germline",
    "framework_identity",
    "cdr_identity",
    "humanness_risk_class",
    "humanness_interpretation",
    "humanness_notes",
    "merge_status",
    "merge_warning",
]

ANTIBODY_HUMANNESS_COLUMNS = [
    "humanness_available",
    "chain_humanness_risk_classes",
    "max_humanness_risk_class",
    "closest_human_germlines",
    "non_human_like_chain_count",
    "high_humanness_risk_chain_count",
    "humanness_summary",
    "combined_developability_humanness_flag",
]

RISK_RANK = {"Low": 0, "Medium": 1, "High": 2}
NON_HUMAN_SPECIES = {"mouse", "rabbit", "rat", "camelid"}

HUMANNESS_INTERPRETATION_SUFFIX = (
    "Human-likeness is not equivalent to clinical immunogenicity prediction."
)


def _empty_results_df() -> pd.DataFrame:
    return pd.DataFrame(columns=HUMANNESS_RESULT_COLUMNS)


def _to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _risk_from_percent(value: object, low_cutoff: float, medium_cutoff: float) -> str | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if pd.isna(numeric):
        return None
    if numeric >= low_cutoff:
        return "Low"
    if numeric >= medium_cutoff:
        return "Medium"
    return "High"


def _max_humanness_risk(classes: list[str]) -> str:
    clean = [risk for risk in classes if risk in RISK_RANK]
    if not clean:
        return "Unknown"
    return max(clean, key=lambda item: RISK_RANK[item])


def _classify_humanness(row: pd.Series) -> tuple[str, str]:
    candidates: list[str] = []
    evidence: list[str] = []

    percentile_risk = _risk_from_percent(row.get("humanness_percentile"), 70, 40)
    if percentile_risk:
        candidates.append(percentile_risk)
        evidence.append(f"humanness_percentile={row.get('humanness_percentile')}")

    framework_risk = _risk_from_percent(row.get("framework_identity"), 85, 70)
    if framework_risk:
        candidates.append(framework_risk)
        evidence.append(f"framework_identity={row.get('framework_identity')}")

    germline_identity_risk = _risk_from_percent(row.get("identity_to_human_germline"), 85, 70)
    if germline_identity_risk:
        candidates.append(germline_identity_risk)
        evidence.append(f"identity_to_human_germline={row.get('identity_to_human_germline')}")

    species = normalize_text(row.get("closest_species")).lower()
    if species == "human":
        candidates.append("Low")
        evidence.append("closest_species=human")
    elif species in NON_HUMAN_SPECIES:
        candidates.append("Medium")
        evidence.append(f"closest_species={species}")
    elif species == "unknown":
        evidence.append("closest_species=unknown")

    risk_class = _max_humanness_risk(candidates)
    if risk_class == "Unknown":
        interpretation = "No sufficient imported humanness/germline metrics were available for classification."
    else:
        interpretation = (
            f"Imported humanness/germline metrics indicate {risk_class.lower()} human-likeness risk "
            f"based on conservative review of: {join_messages(evidence)}."
        )
    return risk_class, f"{interpretation} {HUMANNESS_INTERPRETATION_SUFFIX}"


def _read_uploaded_table(uploaded_file: BinaryIO | str | Path) -> pd.DataFrame:
    name = getattr(uploaded_file, "name", "") or str(uploaded_file)
    suffix = Path(name).suffix.lower()
    try:
        if suffix == ".csv":
            return pd.read_csv(uploaded_file)
        if suffix in {".xlsx", ".xls"}:
            return pd.read_excel(uploaded_file)
    except pd.errors.EmptyDataError as exc:
        raise ValueError("The humanness/germline file is empty.") from exc
    except Exception as exc:
        raise ValueError(f"Could not read the humanness/germline file: {exc}") from exc
    raise ValueError("Unsupported humanness/germline file format. Please upload CSV or XLSX.")


def load_humanness_file(uploaded_file) -> pd.DataFrame:
    """Read a user-provided optional humanness/germline CSV or XLSX file."""
    if uploaded_file is None:
        return pd.DataFrame(columns=HUMANNESS_INPUT_COLUMNS)
    return validate_humanness_df(_read_uploaded_table(uploaded_file))


def validate_humanness_df(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize optional humanness/germline result columns and risk labels."""
    if df is None or df.empty:
        return pd.DataFrame(columns=HUMANNESS_RESULT_COLUMNS)

    table = df.dropna(how="all").copy()
    table.columns = [normalize_text(column).strip() for column in table.columns]
    for column in HUMANNESS_INPUT_COLUMNS:
        if column not in table.columns:
            table[column] = pd.NA

    table = table[HUMANNESS_INPUT_COLUMNS].copy()
    for column in [
        "antibody_id",
        "chain_id",
        "sequence_scope",
        "humanness_tool",
        "closest_human_germline",
        "closest_species",
        "v_gene",
        "j_gene",
        "humanness_notes",
    ]:
        table[column] = table[column].map(normalize_text)

    for column in [
        "humanness_score",
        "humanness_percentile",
        "identity_to_human_germline",
        "framework_identity",
        "cdr_identity",
    ]:
        table[column] = _to_numeric(table[column])

    risk_and_interpretation = table.apply(_classify_humanness, axis=1)
    table["humanness_risk_class"] = risk_and_interpretation.map(lambda item: item[0])
    table["humanness_interpretation"] = risk_and_interpretation.map(lambda item: item[1])
    table["merge_status"] = "unmerged"
    table["merge_warning"] = ""

    return table[HUMANNESS_RESULT_COLUMNS]


def _chain_key(row: pd.Series) -> tuple[str, str, str]:
    return (
        normalize_text(row.get("antibody_id")),
        normalize_text(row.get("chain_id")).lower(),
        normalize_text(row.get("sequence_scope")).lower(),
    )


def _build_chain_lookup(chain_scores_df: pd.DataFrame) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    by_chain: set[tuple[str, str]] = set()
    by_scope: set[tuple[str, str]] = set()
    if chain_scores_df is None or chain_scores_df.empty:
        return by_chain, by_scope
    for _, row in chain_scores_df.iterrows():
        antibody_id = normalize_text(row.get("antibody_id"))
        chain_id = normalize_text(row.get("chain_id")).lower()
        sequence_scope = normalize_text(row.get("sequence_scope")).lower()
        if antibody_id and chain_id:
            by_chain.add((antibody_id, chain_id))
        if antibody_id and sequence_scope:
            by_scope.add((antibody_id, sequence_scope))
    return by_chain, by_scope


def _merge_status(row: pd.Series, by_chain: set[tuple[str, str]], by_scope: set[tuple[str, str]]) -> tuple[str, str]:
    antibody_id = normalize_text(row.get("antibody_id"))
    chain_id = normalize_text(row.get("chain_id")).lower()
    sequence_scope = normalize_text(row.get("sequence_scope")).lower()
    warnings: list[str] = []
    if not antibody_id:
        return "unmatched", "Missing antibody_id in imported humanness row"
    if chain_id and (antibody_id, chain_id) in by_chain:
        return "matched_by_antibody_id_chain_id", ""
    if sequence_scope and (antibody_id, sequence_scope) in by_scope:
        return "matched_by_antibody_id_sequence_scope", ""
    if not chain_id and not sequence_scope:
        warnings.append("Missing both chain_id and sequence_scope")
    else:
        warnings.append("No matching analyzed chain found for imported humanness row")
    return "unmatched", join_messages(warnings)


def _records_for_chain(row: pd.Series, humanness_df: pd.DataFrame) -> pd.DataFrame:
    if humanness_df.empty:
        return pd.DataFrame()
    antibody_id, chain_id, sequence_scope = _chain_key(row)
    chain_matches = humanness_df[
        (humanness_df["antibody_id"].astype(str) == antibody_id)
        & (humanness_df["chain_id"].astype(str).str.lower() == chain_id)
    ]
    if not chain_matches.empty:
        return chain_matches
    return humanness_df[
        (humanness_df["antibody_id"].astype(str) == antibody_id)
        & (humanness_df["sequence_scope"].astype(str).str.lower() == sequence_scope)
    ]


def _combined_flag(max_humanness_risk: str, developability_score: object) -> str:
    try:
        score = float(developability_score)
    except (TypeError, ValueError):
        score = 0.0
    high_dev = score > 8
    if max_humanness_risk == "High" and high_dev:
        return "High priority engineering review"
    if max_humanness_risk == "High":
        return "Humanness review recommended"
    if high_dev:
        return "Developability review recommended"
    return "No major sequence-level humanness/developability flag in MVP assessment"


def _append_message(existing: object, message: str) -> str:
    return join_messages([normalize_text(existing), message])


def summarize_humanness_by_antibody(merged_chain_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate imported humanness results at antibody level."""
    if merged_chain_df is None or merged_chain_df.empty or "merge_status" not in merged_chain_df.columns:
        return pd.DataFrame(columns=["antibody_id", *ANTIBODY_HUMANNESS_COLUMNS])

    matched = merged_chain_df[merged_chain_df["merge_status"].astype(str).str.startswith("matched")].copy()
    if matched.empty:
        return pd.DataFrame(columns=["antibody_id", *ANTIBODY_HUMANNESS_COLUMNS])

    records: list[dict[str, object]] = []
    for antibody_id, rows in matched.groupby("antibody_id", dropna=False):
        risk_classes = rows["humanness_risk_class"].fillna("Unknown").astype(str).tolist()
        max_risk = _max_humanness_risk(risk_classes)
        germlines = [
            value
            for value in rows["closest_human_germline"].map(normalize_text).tolist()
            if value
        ]
        non_human_like_count = int(rows["humanness_risk_class"].astype(str).isin(["Medium", "High"]).sum())
        high_count = int(rows["humanness_risk_class"].astype(str).eq("High").sum())
        records.append(
            {
                "antibody_id": antibody_id,
                "humanness_available": True,
                "chain_humanness_risk_classes": join_messages(
                    f"{normalize_text(row.get('chain_id')) or normalize_text(row.get('sequence_scope'))}: {row.get('humanness_risk_class')}"
                    for _, row in rows.iterrows()
                ),
                "max_humanness_risk_class": max_risk,
                "closest_human_germlines": join_messages(germlines) or "",
                "non_human_like_chain_count": non_human_like_count,
                "high_humanness_risk_chain_count": high_count,
                "humanness_summary": (
                    f"Imported humanness/germline assessment available for {len(rows)} chain(s); "
                    f"maximum conservative human-likeness risk class: {max_risk}. "
                    f"{HUMANNESS_INTERPRETATION_SUFFIX}"
                ),
            }
        )
    return pd.DataFrame(records)


def _add_default_humanness_columns(antibody_summary_df: pd.DataFrame) -> pd.DataFrame:
    summary = antibody_summary_df.copy()
    defaults = {
        "humanness_available": False,
        "chain_humanness_risk_classes": "",
        "max_humanness_risk_class": "Unknown",
        "closest_human_germlines": "",
        "non_human_like_chain_count": 0,
        "high_humanness_risk_chain_count": 0,
        "humanness_summary": "No external humanness/germline assessment file was provided or matched.",
        "combined_developability_humanness_flag": "",
    }
    for column, default in defaults.items():
        if column not in summary.columns:
            summary[column] = default
        else:
            summary[column] = summary[column].fillna(default)
    return summary


def merge_humanness_results(
    chain_scores_df: pd.DataFrame,
    region_summary_df: pd.DataFrame,
    antibody_summary_df: pd.DataFrame,
    humanness_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge optional imported humanness results and update antibody summaries."""
    del region_summary_df
    summary = _add_default_humanness_columns(antibody_summary_df if antibody_summary_df is not None else pd.DataFrame())

    if humanness_df is None or humanness_df.empty:
        if not summary.empty:
            summary["combined_developability_humanness_flag"] = summary.apply(
                lambda row: _combined_flag(row.get("max_humanness_risk_class"), row.get("cdr_adjusted_total_risk_score")),
                axis=1,
            )
        return _empty_results_df(), summary

    humanness = validate_humanness_df(humanness_df)
    by_chain, by_scope = _build_chain_lookup(chain_scores_df)
    statuses = humanness.apply(lambda row: _merge_status(row, by_chain, by_scope), axis=1)
    humanness["merge_status"] = statuses.map(lambda item: item[0])
    humanness["merge_warning"] = statuses.map(lambda item: item[1])

    antibody_humanness = summarize_humanness_by_antibody(humanness)
    if not antibody_humanness.empty and not summary.empty:
        summary = summary.drop(columns=[column for column in ANTIBODY_HUMANNESS_COLUMNS if column in summary.columns])
        summary = summary.merge(antibody_humanness, on="antibody_id", how="left")
        summary = _add_default_humanness_columns(summary)

    if not summary.empty:
        summary["combined_developability_humanness_flag"] = summary.apply(
            lambda row: _combined_flag(row.get("max_humanness_risk_class"), row.get("cdr_adjusted_total_risk_score")),
            axis=1,
        )
        low_mask = summary["max_humanness_risk_class"].astype(str).eq("Low")
        high_mask = summary["max_humanness_risk_class"].astype(str).eq("High")
        summary.loc[low_mask, "strengths"] = summary.loc[low_mask, "strengths"].map(
            lambda value: _append_message(
                value,
                "Human-likeness metrics suggest low humanness risk based on imported assessment",
            )
        )
        summary.loc[high_mask, "weaknesses"] = summary.loc[high_mask, "weaknesses"].map(
            lambda value: _append_message(
                value,
                "Imported humanness metrics suggest high human-likeness risk; humanization review may be needed",
            )
        )
        recommendation_mask = high_mask & summary.get("total_cdr_liabilities", pd.Series(0, index=summary.index)).fillna(0).gt(0)
        summary.loc[recommendation_mask, "recommendation"] = summary.loc[recommendation_mask, "recommendation"].map(
            lambda value: _append_message(
                value,
                "Prioritize humanization/developability review before progression. Avoid automatic CDR mutation decisions without binding and functional data.",
            )
        )

    return humanness[HUMANNESS_RESULT_COLUMNS], summary
