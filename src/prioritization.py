"""Candidate-level prioritization and decision matrix."""

from __future__ import annotations

import pandas as pd

from .utils import join_messages, normalize_text


RANKING_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "number_of_chains",
    "total_risk_score",
    "cdr_adjusted_total_risk_score",
    "max_chain_risk_class",
    "cdr_adjusted_max_risk_class",
    "total_cdr_liabilities",
    "total_fr_liabilities",
    "max_humanness_risk_class",
    "high_humanness_risk_chain_count",
    "non_human_like_chain_count",
    "combined_developability_humanness_flag",
    "final_priority_score",
    "final_priority_class",
    "decision_label",
    "review_reason",
    "next_step_recommendation",
]

HUMANNESS_RISK_RANK = {"Low": 0, "Medium": 1, "High": 2, "Unknown": -1}


def _to_float(value: object, default: float = 0.0) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return default
    if pd.isna(numeric):
        return default
    return numeric


def _to_int(value: object, default: int = 0) -> int:
    return int(round(_to_float(value, float(default))))


def _priority_class(score: float) -> str:
    if score >= 80:
        return "A"
    if score >= 60:
        return "B"
    if score >= 40:
        return "C"
    return "D"


def _decision_label(priority_class: str) -> str:
    return {
        "A": "Advance to next-stage screening",
        "B": "Advance with targeted review",
        "C": "Engineering review before progression",
        "D": "Deprioritize or redesign before progression",
    }[priority_class]


def _next_step(priority_class: str) -> str:
    return {
        "A": "Proceed to expression and functional validation",
        "B": "Proceed with targeted sequence liability and humanness review",
        "C": "Perform engineering review before expression scale-up",
        "D": "Consider redesign, deprioritization, or additional humanization/developability review",
    }[priority_class]


def _risk_penalty(risk_class: object) -> int:
    risk = normalize_text(risk_class)
    if risk == "Medium":
        return 10
    if risk == "High":
        return 25
    return 0


def _max_humanness_risk(classes: pd.Series) -> str:
    clean = [normalize_text(value) for value in classes.tolist()]
    valid = [value for value in clean if value in HUMANNESS_RISK_RANK and value != "Unknown"]
    if not valid:
        return "Unknown"
    return max(valid, key=lambda value: HUMANNESS_RISK_RANK[value])


def _matched_humanness(humanness_results_df: pd.DataFrame | None) -> pd.DataFrame:
    if humanness_results_df is None or humanness_results_df.empty:
        return pd.DataFrame()
    if "merge_status" not in humanness_results_df.columns:
        return humanness_results_df.copy()
    return humanness_results_df[
        humanness_results_df["merge_status"].astype(str).str.startswith("matched")
    ].copy()


def _aggregate_humanness(humanness_results_df: pd.DataFrame | None) -> pd.DataFrame:
    humanness = _matched_humanness(humanness_results_df)
    if humanness.empty or "antibody_id" not in humanness.columns:
        return pd.DataFrame(
            columns=[
                "antibody_id",
                "max_humanness_risk_class",
                "high_humanness_risk_chain_count",
                "non_human_like_chain_count",
            ]
        )

    records: list[dict[str, object]] = []
    for antibody_id, rows in humanness.groupby("antibody_id", dropna=False):
        risk_classes = rows.get("humanness_risk_class", pd.Series(dtype=str)).astype(str)
        records.append(
            {
                "antibody_id": antibody_id,
                "max_humanness_risk_class": _max_humanness_risk(risk_classes),
                "high_humanness_risk_chain_count": int(risk_classes.eq("High").sum()),
                "non_human_like_chain_count": int(risk_classes.isin(["Medium", "High"]).sum()),
            }
        )
    return pd.DataFrame(records)


def _liability_rows_for_antibody(
    antibody_id: object,
    liability_region_map_df: pd.DataFrame | None,
) -> pd.DataFrame:
    if (
        liability_region_map_df is None
        or liability_region_map_df.empty
        or "antibody_id" not in liability_region_map_df.columns
    ):
        return pd.DataFrame()
    return liability_region_map_df[liability_region_map_df["antibody_id"] == antibody_id]


def _chain_rows_for_antibody(antibody_id: object, chain_scores_df: pd.DataFrame | None) -> pd.DataFrame:
    if chain_scores_df is None or chain_scores_df.empty or "antibody_id" not in chain_scores_df.columns:
        return pd.DataFrame()
    return chain_scores_df[chain_scores_df["antibody_id"] == antibody_id]


def _review_reason(row: pd.Series, region_rows: pd.DataFrame, chain_rows: pd.DataFrame) -> str:
    reasons: list[str] = []
    if _to_int(row.get("total_cdr_liabilities")) > 0:
        reasons.append("CDR-mapped liability detected")
    if normalize_text(row.get("max_humanness_risk_class")) == "High":
        reasons.append("High humanness risk based on imported assessment")
    if _to_float(row.get("cdr_adjusted_total_risk_score")) > 8:
        reasons.append("High CDR-adjusted developability burden")

    if not region_rows.empty and {"risk_type", "cdr_or_fr"}.issubset(region_rows.columns):
        cdr_risk_types = region_rows[region_rows["cdr_or_fr"].astype(str).eq("CDR")]["risk_type"].astype(str)
        if cdr_risk_types.eq("N_glycosylation").any():
            reasons.append("CDR N-glycosylation motif requires review")
        hydrophobic_count = int(region_rows["risk_type"].astype(str).eq("hydrophobic_patch_proxy").sum())
    else:
        hydrophobic_count = 0
    if not chain_rows.empty and "main_chain_liabilities" in chain_rows.columns:
        hydrophobic_count = max(
            hydrophobic_count,
            int(chain_rows["main_chain_liabilities"].astype(str).str.contains("hydrophobic_patch_proxy", na=False).sum()),
        )
    if hydrophobic_count >= 2:
        reasons.append("Sequence-level hydrophobic patch proxy detected")

    return (
        join_messages(reasons)
        or "No major sequence-level developability or humanness flag detected in MVP assessment"
    )


def _base_ranking_table(
    antibody_summary_df: pd.DataFrame,
    humanness_results_df: pd.DataFrame | None,
) -> pd.DataFrame:
    ranking = antibody_summary_df.copy()
    for column in [
        "antibody_id",
        "molecule_format",
        "number_of_chains",
        "total_risk_score",
        "cdr_adjusted_total_risk_score",
        "max_chain_risk_class",
        "cdr_adjusted_max_risk_class",
        "total_cdr_liabilities",
        "total_fr_liabilities",
        "max_humanness_risk_class",
        "high_humanness_risk_chain_count",
        "non_human_like_chain_count",
        "combined_developability_humanness_flag",
    ]:
        if column not in ranking.columns:
            ranking[column] = pd.NA

    if ranking["max_humanness_risk_class"].isna().all() or ranking["max_humanness_risk_class"].astype(str).eq("").all():
        aggregated = _aggregate_humanness(humanness_results_df)
        if not aggregated.empty:
            ranking = ranking.drop(
                columns=[
                    column
                    for column in [
                        "max_humanness_risk_class",
                        "high_humanness_risk_chain_count",
                        "non_human_like_chain_count",
                    ]
                    if column in ranking.columns
                ]
            ).merge(aggregated, on="antibody_id", how="left")

    ranking["max_humanness_risk_class"] = ranking["max_humanness_risk_class"].fillna("Unknown").replace("", "Unknown")
    ranking["high_humanness_risk_chain_count"] = ranking["high_humanness_risk_chain_count"].fillna(0)
    ranking["non_human_like_chain_count"] = ranking["non_human_like_chain_count"].fillna(0)
    ranking["combined_developability_humanness_flag"] = (
        ranking["combined_developability_humanness_flag"]
        .fillna("")
        .replace("", "No major sequence-level humanness/developability flag in MVP assessment")
    )
    return ranking


def build_candidate_ranking(
    antibody_summary_df: pd.DataFrame,
    chain_scores_df: pd.DataFrame,
    liability_region_map_df: pd.DataFrame | None = None,
    humanness_results_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build rule-based antibody candidate ranking from pipeline result tables."""
    if antibody_summary_df is None or antibody_summary_df.empty:
        return pd.DataFrame(columns=RANKING_COLUMNS)

    ranking = _base_ranking_table(antibody_summary_df, humanness_results_df)
    records: list[dict[str, object]] = []

    for _, row in ranking.iterrows():
        score = 100.0
        score -= _to_float(row.get("cdr_adjusted_total_risk_score")) * 3
        score -= _to_int(row.get("total_cdr_liabilities")) * 5
        score -= _to_int(row.get("total_fr_liabilities")) * 1
        score -= _risk_penalty(row.get("max_chain_risk_class"))
        score -= _risk_penalty(row.get("cdr_adjusted_max_risk_class"))
        score -= _risk_penalty(row.get("max_humanness_risk_class"))
        score -= _to_int(row.get("high_humanness_risk_chain_count")) * 10
        score -= _to_int(row.get("non_human_like_chain_count")) * 5
        score = max(0.0, min(100.0, score))
        priority_class = _priority_class(score)
        region_rows = _liability_rows_for_antibody(row.get("antibody_id"), liability_region_map_df)
        chain_rows = _chain_rows_for_antibody(row.get("antibody_id"), chain_scores_df)

        records.append(
            {
                "antibody_id": row.get("antibody_id", ""),
                "molecule_format": row.get("molecule_format", ""),
                "number_of_chains": _to_int(row.get("number_of_chains")),
                "total_risk_score": round(_to_float(row.get("total_risk_score")), 2),
                "cdr_adjusted_total_risk_score": round(_to_float(row.get("cdr_adjusted_total_risk_score")), 2),
                "max_chain_risk_class": normalize_text(row.get("max_chain_risk_class")) or "Low",
                "cdr_adjusted_max_risk_class": normalize_text(row.get("cdr_adjusted_max_risk_class")) or "Low",
                "total_cdr_liabilities": _to_int(row.get("total_cdr_liabilities")),
                "total_fr_liabilities": _to_int(row.get("total_fr_liabilities")),
                "max_humanness_risk_class": normalize_text(row.get("max_humanness_risk_class")) or "Unknown",
                "high_humanness_risk_chain_count": _to_int(row.get("high_humanness_risk_chain_count")),
                "non_human_like_chain_count": _to_int(row.get("non_human_like_chain_count")),
                "combined_developability_humanness_flag": row.get("combined_developability_humanness_flag", ""),
                "final_priority_score": round(score, 2),
                "final_priority_class": priority_class,
                "decision_label": _decision_label(priority_class),
                "review_reason": _review_reason(row, region_rows, chain_rows),
                "next_step_recommendation": _next_step(priority_class),
            }
        )

    result = pd.DataFrame(records, columns=RANKING_COLUMNS)
    if not result.empty:
        result = result.sort_values(
            ["final_priority_score", "antibody_id"],
            ascending=[False, True],
            kind="mergesort",
        ).reset_index(drop=True)
    return result
