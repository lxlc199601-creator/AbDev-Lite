"""Antibody-level result summarization."""

from __future__ import annotations

import pandas as pd

from .utils import join_messages, max_risk_class, normalize_text


def _top_liabilities(df: pd.DataFrame) -> str:
    if df.empty:
        return "None detected"
    counts = df["risk_type"].value_counts().head(5)
    return "; ".join(f"{risk_type} ({count})" for risk_type, count in counts.items())


def _region_warning(region_rows: pd.DataFrame) -> str:
    if region_rows.empty:
        return ""
    warnings = region_rows.get("region_mapping_warning", pd.Series(dtype=str)).astype(str)
    statuses = region_rows.get("region_mapping_status", pd.Series(dtype=str)).astype(str)
    if statuses.isin(["FAIL", "WARNING"]).any() or warnings.str.len().gt(0).any():
        return "Some liability sites could not be confidently mapped to IMGT CDR/FR regions"
    return ""


SUMMARY_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "number_of_chains",
    "sequence_scopes",
    "variable_region_chain_count",
    "full_length_chain_count",
    "total_risk_score",
    "average_chain_risk_score",
    "max_chain_risk_class",
    "low_chain_count",
    "medium_chain_count",
    "high_chain_count",
    "major_liabilities",
    "total_cdr_liabilities",
    "total_fr_liabilities",
    "total_unknown_region_liabilities",
    "major_cdr_liabilities",
    "major_fr_liabilities",
    "cdr_region_warning",
    "cdr_adjusted_total_risk_score",
    "cdr_adjusted_max_risk_class",
    "strengths",
    "weaknesses",
    "recommendation",
]

BSAB_NOTICE = (
    "BsAb format detected. This MVP evaluates sequence-level liabilities of each variable region chain independently. "
    "It does not yet evaluate chain pairing, heterodimerization, Fc engineering, linker geometry, or full molecule architecture."
)


def _strengths(qc_rows: pd.DataFrame, liability_rows: pd.DataFrame, total_score: float, total_cdr_liabilities: int = 0) -> str:
    strengths: list[str] = []
    variable_liabilities = liability_rows[
        liability_rows["region_type"].astype(str).str.lower().eq("variable_region")
    ] if not liability_rows.empty else pd.DataFrame()
    variable_risk_types = set(variable_liabilities["risk_type"]) if not variable_liabilities.empty else set()
    risk_types = set(liability_rows["risk_type"]) if not liability_rows.empty else set()
    if not qc_rows.empty and (qc_rows["qc_status"] == "PASS").all():
        strengths.append("Input sequences passed basic QC")
    if total_score <= 3:
        strengths.append("Low sequence-level liability burden")
    if "N_glycosylation" not in variable_risk_types:
        strengths.append("No obvious N-glycosylation motif detected in variable region")
    odd_cys = (
        liability_rows["motif"].astype(str).eq("potential_unpaired_cysteine").any()
        if not liability_rows.empty
        else False
    )
    if not odd_cys:
        strengths.append("No obvious unpaired cysteine proxy detected")
    if "hydrophobic_patch_proxy" not in risk_types:
        strengths.append("No long hydrophobic segment detected by sequence-level proxy")
    if total_cdr_liabilities == 0:
        strengths.append("No obvious liability motif mapped to CDR regions by IMGT-based computational assignment")
    return join_messages(strengths) or "No specific strengths assigned by MVP rules"


def _weaknesses(
    qc_rows: pd.DataFrame,
    liability_rows: pd.DataFrame,
    molecule_format: str,
    total_cdr_liabilities: int = 0,
) -> str:
    weaknesses: list[str] = []
    risk_types = set(liability_rows["risk_type"]) if not liability_rows.empty else set()
    labels = {
        "N_glycosylation": "N-glycosylation motif detected",
        "deamidation": "Deamidation motif detected",
        "isomerization": "Isomerization motif detected",
        "oxidation": "Oxidation-prone residue detected",
        "hydrophobic_patch_proxy": "Long hydrophobic segment detected by sequence-level proxy",
    }
    for risk_type, label in labels.items():
        if risk_type in risk_types:
            weaknesses.append(label)
    odd_cys = (
        liability_rows["motif"].astype(str).eq("potential_unpaired_cysteine").any()
        if not liability_rows.empty
        else False
    )
    if odd_cys:
        weaknesses.append("Odd cysteine count proxy detected")
    if not qc_rows.empty and qc_rows["variable_region_length_warning"].astype(str).str.len().gt(0).any():
        weaknesses.append("Length warning detected")
    if not qc_rows.empty and (
        qc_rows["illegal_characters"].astype(str).str.len().gt(0).any() or (qc_rows["qc_status"] == "FAIL").any()
    ):
        weaknesses.append("Illegal characters or QC failure detected")
    if normalize_text(molecule_format).lower() == "bsab":
        weaknesses.append(BSAB_NOTICE)
    if total_cdr_liabilities > 0:
        weaknesses.append(
            "Liability motifs detected in CDR regions; these may require careful review because CDRs are directly involved in antigen recognition"
        )
    return join_messages(weaknesses) or "No major weaknesses assigned by MVP rules"


def summarize_antibodies(
    qc_df: pd.DataFrame,
    liabilities_df: pd.DataFrame,
    chain_scores_df: pd.DataFrame,
    liability_region_map_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Summarize chain results at antibody level."""
    records: list[dict[str, object]] = []
    if qc_df.empty:
        return pd.DataFrame(records, columns=SUMMARY_COLUMNS)

    for antibody_id, chains in qc_df.groupby("antibody_id", dropna=False):
        molecule_format = normalize_text(chains["molecule_format"].iloc[0])
        qc_rows = chains
        liability_rows = liabilities_df[liabilities_df["antibody_id"] == antibody_id] if not liabilities_df.empty else pd.DataFrame()
        score_rows = chain_scores_df[chain_scores_df["antibody_id"] == antibody_id]
        region_rows = (
            liability_region_map_df[liability_region_map_df["antibody_id"] == antibody_id]
            if isinstance(liability_region_map_df, pd.DataFrame)
            and not liability_region_map_df.empty
            and "antibody_id" in liability_region_map_df.columns
            else pd.DataFrame()
        )

        total_score = float(score_rows["risk_score"].sum()) if not score_rows.empty else 0.0
        average_score = float(score_rows["risk_score"].mean()) if not score_rows.empty else 0.0
        max_class = max_risk_class(score_rows["risk_class"].tolist()) if not score_rows.empty else "Low"
        adjusted_total_score = (
            float(score_rows["cdr_adjusted_risk_score"].sum())
            if not score_rows.empty and "cdr_adjusted_risk_score" in score_rows.columns
            else total_score
        )
        adjusted_max_class = (
            max_risk_class(score_rows["cdr_adjusted_risk_class"].tolist())
            if not score_rows.empty and "cdr_adjusted_risk_class" in score_rows.columns
            else max_class
        )
        total_cdr_liabilities = int(region_rows["cdr_or_fr"].astype(str).eq("CDR").sum()) if not region_rows.empty else 0
        total_fr_liabilities = int(region_rows["cdr_or_fr"].astype(str).eq("FR").sum()) if not region_rows.empty else 0
        total_unknown_region_liabilities = (
            int(region_rows["cdr_or_fr"].astype(str).isin(["Unknown", "Boundary"]).sum()) if not region_rows.empty else 0
        )
        cdr_rows = region_rows[region_rows["cdr_or_fr"].astype(str).eq("CDR")] if not region_rows.empty else pd.DataFrame()
        fr_rows = region_rows[region_rows["cdr_or_fr"].astype(str).eq("FR")] if not region_rows.empty else pd.DataFrame()
        region_types = chains["region_type"].astype(str).str.lower()
        variable_count = region_types.eq("variable_region").sum()
        full_count = region_types.eq("full_length").sum()

        recommendation = "Recommend for next-stage screening"
        if adjusted_total_score > 8:
            recommendation = "High sequence-level developability risk; recommend engineering review before progression"
        elif adjusted_total_score > 3:
            recommendation = "Proceed with caution; recommend targeted liability review"
        recommendation_parts = [recommendation]
        if adjusted_total_score >= total_score + 2:
            recommendation_parts.append("CDR-mapped liabilities increase review priority.")
        if not cdr_rows.empty:
            cdr3_major = cdr_rows[
                cdr_rows["mapped_regions"].astype(str).str.contains("CDR3", na=False)
                & cdr_rows["risk_type"].astype(str).isin(["N_glycosylation", "deamidation", "isomerization"])
            ]
            if not cdr3_major.empty:
                recommendation_parts.append(
                    "CDR3 liability detected; targeted engineering review is recommended if binding activity permits mutation."
                )

        records.append(
            {
                "antibody_id": antibody_id,
                "molecule_format": molecule_format,
                "number_of_chains": len(chains),
                "sequence_scopes": ", ".join(chains["sequence_scope"].astype(str).tolist()),
                "variable_region_chain_count": int(variable_count),
                "full_length_chain_count": int(full_count),
                "total_risk_score": round(total_score, 2),
                "average_chain_risk_score": round(average_score, 2),
                "max_chain_risk_class": max_class,
                "low_chain_count": int((score_rows["risk_class"] == "Low").sum()) if not score_rows.empty else 0,
                "medium_chain_count": int((score_rows["risk_class"] == "Medium").sum()) if not score_rows.empty else 0,
                "high_chain_count": int((score_rows["risk_class"] == "High").sum()) if not score_rows.empty else 0,
                "major_liabilities": _top_liabilities(liability_rows),
                "total_cdr_liabilities": total_cdr_liabilities,
                "total_fr_liabilities": total_fr_liabilities,
                "total_unknown_region_liabilities": total_unknown_region_liabilities,
                "major_cdr_liabilities": _top_liabilities(cdr_rows),
                "major_fr_liabilities": _top_liabilities(fr_rows),
                "cdr_region_warning": _region_warning(region_rows),
                "cdr_adjusted_total_risk_score": round(adjusted_total_score, 2),
                "cdr_adjusted_max_risk_class": adjusted_max_class,
                "strengths": _strengths(qc_rows, liability_rows, total_score, total_cdr_liabilities),
                "weaknesses": _weaknesses(qc_rows, liability_rows, molecule_format, total_cdr_liabilities),
                "recommendation": join_messages(recommendation_parts),
            }
        )

    return pd.DataFrame(records, columns=SUMMARY_COLUMNS)
