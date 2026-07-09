"""Chain-level sequence risk scoring."""

from __future__ import annotations

import pandas as pd

from .utils import join_messages, normalize_text, risk_class_from_score


RISK_WEIGHTS = {
    "deamidation": 1.0,
    "isomerization": 1.0,
    "N_glycosylation": 2.0,
    "clipping_risk": 0.5,
    "acid_sensitive": 0.5,
    "hydrophobic_patch_proxy": 1.5,
    "low_complexity": 1.5,
}

CDR_EXTRA_WEIGHTS = {
    "deamidation": 1.0,
    "isomerization": 1.0,
    "N_glycosylation": 2.0,
    "hydrophobic_patch_proxy": 1.0,
}


SCORE_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "chain_id",
    "chain_type",
    "region_type",
    "sequence_scope",
    "risk_score",
    "risk_class",
    "cdr_adjusted_risk_score",
    "cdr_adjusted_risk_class",
    "cdr_liability_count",
    "fr_liability_count",
    "unknown_region_liability_count",
    "main_chain_liabilities",
    "chain_recommendation",
]


def _chain_key(row: pd.Series) -> tuple[object, object, object]:
    return (row.get("antibody_id"), row.get("chain_id"), row.get("sequence_scope"))


def _liability_score(site: pd.Series) -> float:
    risk_type = normalize_text(site.get("risk_type"))
    motif = normalize_text(site.get("motif"))
    if risk_type == "oxidation":
        return 0.5 if motif == "M" else 1.0 if motif == "W" else 0.0
    if risk_type == "cysteine_risk":
        if motif == "potential_unpaired_cysteine":
            return 3.0
        if motif == "high_cysteine_count_warning":
            return 2.0
        return 0.0
    return RISK_WEIGHTS.get(risk_type, 0.0)


def _cdr_extra_score(site: pd.Series) -> float:
    if normalize_text(site.get("cdr_or_fr")) != "CDR":
        return 0.0
    risk_type = normalize_text(site.get("risk_type"))
    motif = normalize_text(site.get("motif"))
    if risk_type == "oxidation" and motif == "W":
        return 1.0
    return CDR_EXTRA_WEIGHTS.get(risk_type, 0.0)


def _recommendation(risk_class: str) -> str:
    return {
        "Low": "Suitable for next-stage screening based on sequence-level MVP assessment",
        "Medium": "Proceed with caution; review highlighted liability motifs",
        "High": "Engineering review recommended before progression",
    }[risk_class]


def score_chains(
    qc_df: pd.DataFrame,
    liabilities_df: pd.DataFrame,
    liability_region_map_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Calculate chain-level MVP risk scores from QC and liability tables."""
    records: list[dict[str, object]] = []
    if qc_df.empty:
        return pd.DataFrame(records, columns=SCORE_COLUMNS)

    liability_lookup: dict[tuple[object, object, object], pd.DataFrame] = {}
    if not liabilities_df.empty:
        for key, group in liabilities_df.groupby(["antibody_id", "chain_id", "sequence_scope"], dropna=False):
            liability_lookup[key] = group

    region_lookup: dict[tuple[object, object, object], pd.DataFrame] = {}
    region_map = liability_region_map_df if isinstance(liability_region_map_df, pd.DataFrame) else pd.DataFrame()
    if not region_map.empty:
        for key, group in region_map.groupby(["antibody_id", "chain_id", "sequence_scope"], dropna=False):
            region_lookup[key] = group

    for _, row in qc_df.iterrows():
        sequence = normalize_text(row.get("cleaned_sequence"))
        score = 0.0
        cdr_extra_score = 0.0
        main_liabilities: list[str] = []
        cdr_liability_count = 0
        fr_liability_count = 0
        unknown_region_liability_count = 0

        key = _chain_key(row)
        chain_liabilities = liability_lookup.get(key, pd.DataFrame())
        if not chain_liabilities.empty:
            for _, site in chain_liabilities.iterrows():
                score += _liability_score(site)
            main_liabilities = chain_liabilities["risk_type"].value_counts().head(5).index.astype(str).tolist()

        mapped_liabilities = region_lookup.get(key, pd.DataFrame())
        if not mapped_liabilities.empty:
            region_categories = mapped_liabilities["cdr_or_fr"].astype(str)
            cdr_liability_count = int(region_categories.eq("CDR").sum())
            fr_liability_count = int(region_categories.eq("FR").sum())
            mapping_statuses = mapped_liabilities["region_mapping_status"].astype(str)
            unknown_region_liability_count = int(
                (region_categories.isin(["Unknown", "Boundary"]) | mapping_statuses.isin(["FAIL", "WARNING"])).sum()
            )
            for _, site in mapped_liabilities.iterrows():
                cdr_extra_score += _cdr_extra_score(site)

        if normalize_text(row.get("variable_region_length_warning")):
            score += 2.0
            main_liabilities.append("length_warning")
        if normalize_text(row.get("illegal_characters")):
            score += 5.0
            main_liabilities.append("illegal_characters")
        if not sequence:
            score += 10.0
            main_liabilities.append("empty_sequence")

        risk_class = risk_class_from_score(score)
        cdr_adjusted_score = score + cdr_extra_score
        cdr_adjusted_risk_class = risk_class_from_score(cdr_adjusted_score)

        records.append(
            {
                "antibody_id": row.get("antibody_id", ""),
                "molecule_format": row.get("molecule_format", ""),
                "chain_id": row.get("chain_id", ""),
                "chain_type": row.get("chain_type", ""),
                "region_type": row.get("region_type", ""),
                "sequence_scope": row.get("sequence_scope", ""),
                "risk_score": round(score, 2),
                "risk_class": risk_class,
                "cdr_adjusted_risk_score": round(cdr_adjusted_score, 2),
                "cdr_adjusted_risk_class": cdr_adjusted_risk_class,
                "cdr_liability_count": cdr_liability_count,
                "fr_liability_count": fr_liability_count,
                "unknown_region_liability_count": unknown_region_liability_count,
                "main_chain_liabilities": join_messages(main_liabilities) or "None detected",
                "chain_recommendation": _recommendation(cdr_adjusted_risk_class),
            }
        )

    return pd.DataFrame(records, columns=SCORE_COLUMNS)
