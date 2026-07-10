"""Formulation-related feature extraction for AbDev-Lite v0.7."""

from __future__ import annotations

import pandas as pd


FEATURE_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "number_of_chains",
    "average_sequence_length",
    "total_sequence_length",
    "average_theoretical_pI",
    "max_theoretical_pI",
    "min_theoretical_pI",
    "average_gravy",
    "max_gravy",
    "average_instability_index",
    "total_cysteine_count",
    "total_methionine_count",
    "total_tryptophan_count",
    "total_asparagine_count",
    "total_aspartic_acid_count",
    "total_lysine_count",
    "total_arginine_count",
    "hydrophobic_residue_fraction_mean",
    "aromatic_residue_fraction_mean",
    "positive_residue_fraction_mean",
    "negative_residue_fraction_mean",
    "total_risk_score",
    "cdr_adjusted_total_risk_score",
    "total_cdr_liabilities",
    "total_fr_liabilities",
    "total_nglycosylation_motifs",
    "total_deamidation_motifs",
    "total_isomerization_motifs",
    "total_oxidation_motifs",
    "total_hydrophobic_patch_proxy",
    "max_humanness_risk_class",
    "final_priority_class",
    "final_priority_score",
]


def _df(value: pd.DataFrame | None) -> pd.DataFrame:
    return value if isinstance(value, pd.DataFrame) else pd.DataFrame()


def _numeric(series: pd.Series | None) -> pd.Series:
    if series is None:
        return pd.Series(dtype=float)
    return pd.to_numeric(series, errors="coerce")


def _first_text(rows: pd.DataFrame, column: str, default: str = "") -> str:
    if rows.empty or column not in rows.columns:
        return default
    values = rows[column].dropna().astype(str)
    values = values[values.str.len() > 0]
    return values.iloc[0] if not values.empty else default


def _safe_mean(rows: pd.DataFrame, column: str) -> float:
    values = _numeric(rows[column]) if column in rows.columns else pd.Series(dtype=float)
    return round(float(values.mean()), 4) if not values.dropna().empty else 0.0


def _safe_max(rows: pd.DataFrame, column: str) -> float:
    values = _numeric(rows[column]) if column in rows.columns else pd.Series(dtype=float)
    return round(float(values.max()), 4) if not values.dropna().empty else 0.0


def _safe_min(rows: pd.DataFrame, column: str) -> float:
    values = _numeric(rows[column]) if column in rows.columns else pd.Series(dtype=float)
    return round(float(values.min()), 4) if not values.dropna().empty else 0.0


def _safe_sum(rows: pd.DataFrame, column: str) -> float:
    values = _numeric(rows[column]) if column in rows.columns else pd.Series(dtype=float)
    return round(float(values.sum()), 4) if not values.dropna().empty else 0.0


def _risk_type_count(rows: pd.DataFrame, risk_type: str) -> int:
    if rows.empty or "risk_type" not in rows.columns:
        return 0
    return int(rows["risk_type"].astype(str).eq(risk_type).sum())


def _aggregate_humanness(humanness_results_df: pd.DataFrame) -> pd.DataFrame:
    if humanness_results_df.empty or "antibody_id" not in humanness_results_df.columns:
        return pd.DataFrame(columns=["antibody_id", "max_humanness_risk_class"])
    risk_rank = {"Low": 0, "Medium": 1, "High": 2, "Unknown": -1}
    records: list[dict[str, object]] = []
    rows = humanness_results_df
    if "merge_status" in rows.columns:
        rows = rows[rows["merge_status"].astype(str).str.startswith("matched")]
    for antibody_id, group in rows.groupby("antibody_id", dropna=False):
        classes = group.get("humanness_risk_class", pd.Series(dtype=str)).fillna("Unknown").astype(str)
        valid = [value for value in classes if value in risk_rank and value != "Unknown"]
        max_class = max(valid, key=lambda value: risk_rank[value]) if valid else "Unknown"
        records.append({"antibody_id": antibody_id, "max_humanness_risk_class": max_class})
    return pd.DataFrame(records)


def build_formulation_features(
    antibody_summary_df: pd.DataFrame,
    chain_properties_df: pd.DataFrame,
    chain_scores_df: pd.DataFrame,
    liability_region_map_df: pd.DataFrame | None = None,
    humanness_results_df: pd.DataFrame | None = None,
    candidate_ranking_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build antibody-level formulation screening features.

    The function is deliberately permissive: missing tables or columns produce
    defaults instead of blocking the pipeline.
    """
    summary = _df(antibody_summary_df)
    properties = _df(chain_properties_df)
    scores = _df(chain_scores_df)
    liability_map = _df(liability_region_map_df)
    humanness = _df(humanness_results_df)
    ranking = _df(candidate_ranking_df)

    antibody_ids: list[object] = []
    for frame in [summary, properties, scores, liability_map, humanness, ranking]:
        if not frame.empty and "antibody_id" in frame.columns:
            antibody_ids.extend(frame["antibody_id"].dropna().tolist())
    antibody_ids = sorted(set(antibody_ids), key=lambda value: str(value))
    if not antibody_ids:
        return pd.DataFrame(columns=FEATURE_COLUMNS)

    humanness_agg = _aggregate_humanness(humanness)
    records: list[dict[str, object]] = []

    for antibody_id in antibody_ids:
        summary_rows = summary[summary["antibody_id"] == antibody_id] if "antibody_id" in summary.columns else pd.DataFrame()
        property_rows = properties[properties["antibody_id"] == antibody_id] if "antibody_id" in properties.columns else pd.DataFrame()
        score_rows = scores[scores["antibody_id"] == antibody_id] if "antibody_id" in scores.columns else pd.DataFrame()
        liability_rows = (
            liability_map[liability_map["antibody_id"] == antibody_id]
            if "antibody_id" in liability_map.columns
            else pd.DataFrame()
        )
        humanness_rows = (
            humanness_agg[humanness_agg["antibody_id"] == antibody_id]
            if "antibody_id" in humanness_agg.columns
            else pd.DataFrame()
        )
        ranking_rows = ranking[ranking["antibody_id"] == antibody_id] if "antibody_id" in ranking.columns else pd.DataFrame()

        number_of_chains = int(_safe_sum(summary_rows, "number_of_chains")) if not summary_rows.empty else len(property_rows)
        if number_of_chains == 0:
            number_of_chains = len(score_rows)

        max_humanness = _first_text(summary_rows, "max_humanness_risk_class", "")
        if not max_humanness:
            max_humanness = _first_text(humanness_rows, "max_humanness_risk_class", "Unknown")

        total_cdr = int(_safe_sum(summary_rows, "total_cdr_liabilities")) if not summary_rows.empty else 0
        total_fr = int(_safe_sum(summary_rows, "total_fr_liabilities")) if not summary_rows.empty else 0
        if total_cdr == 0 and not liability_rows.empty and "cdr_or_fr" in liability_rows.columns:
            total_cdr = int(liability_rows["cdr_or_fr"].astype(str).eq("CDR").sum())
        if total_fr == 0 and not liability_rows.empty and "cdr_or_fr" in liability_rows.columns:
            total_fr = int(liability_rows["cdr_or_fr"].astype(str).eq("FR").sum())

        total_risk = _safe_sum(summary_rows, "total_risk_score")
        if total_risk == 0:
            total_risk = _safe_sum(score_rows, "risk_score")
        cdr_adjusted_total = _safe_sum(summary_rows, "cdr_adjusted_total_risk_score")
        if cdr_adjusted_total == 0:
            cdr_adjusted_total = _safe_sum(score_rows, "cdr_adjusted_risk_score")

        records.append(
            {
                "antibody_id": antibody_id,
                "molecule_format": _first_text(summary_rows, "molecule_format")
                or _first_text(property_rows, "molecule_format"),
                "number_of_chains": int(number_of_chains),
                "average_sequence_length": _safe_mean(property_rows, "sequence_length"),
                "total_sequence_length": int(_safe_sum(property_rows, "sequence_length")),
                "average_theoretical_pI": _safe_mean(property_rows, "theoretical_pI"),
                "max_theoretical_pI": _safe_max(property_rows, "theoretical_pI"),
                "min_theoretical_pI": _safe_min(property_rows, "theoretical_pI"),
                "average_gravy": _safe_mean(property_rows, "gravy"),
                "max_gravy": _safe_max(property_rows, "gravy"),
                "average_instability_index": _safe_mean(property_rows, "instability_index"),
                "total_cysteine_count": int(_safe_sum(property_rows, "cysteine_count")),
                "total_methionine_count": int(_safe_sum(property_rows, "methionine_count")),
                "total_tryptophan_count": int(_safe_sum(property_rows, "tryptophan_count")),
                "total_asparagine_count": int(_safe_sum(property_rows, "asparagine_count")),
                "total_aspartic_acid_count": int(_safe_sum(property_rows, "aspartic_acid_count")),
                "total_lysine_count": int(_safe_sum(property_rows, "lysine_count")),
                "total_arginine_count": int(_safe_sum(property_rows, "arginine_count")),
                "hydrophobic_residue_fraction_mean": _safe_mean(property_rows, "hydrophobic_residue_fraction"),
                "aromatic_residue_fraction_mean": _safe_mean(property_rows, "aromatic_residue_fraction"),
                "positive_residue_fraction_mean": _safe_mean(property_rows, "positive_residue_fraction"),
                "negative_residue_fraction_mean": _safe_mean(property_rows, "negative_residue_fraction"),
                "total_risk_score": round(float(total_risk), 2),
                "cdr_adjusted_total_risk_score": round(float(cdr_adjusted_total), 2),
                "total_cdr_liabilities": total_cdr,
                "total_fr_liabilities": total_fr,
                "total_nglycosylation_motifs": _risk_type_count(liability_rows, "N_glycosylation"),
                "total_deamidation_motifs": _risk_type_count(liability_rows, "deamidation"),
                "total_isomerization_motifs": _risk_type_count(liability_rows, "isomerization"),
                "total_oxidation_motifs": _risk_type_count(liability_rows, "oxidation"),
                "total_hydrophobic_patch_proxy": _risk_type_count(liability_rows, "hydrophobic_patch_proxy"),
                "max_humanness_risk_class": max_humanness or "Unknown",
                "final_priority_class": _first_text(ranking_rows, "final_priority_class", ""),
                "final_priority_score": _safe_max(ranking_rows, "final_priority_score"),
            }
        )

    return pd.DataFrame(records, columns=FEATURE_COLUMNS)
