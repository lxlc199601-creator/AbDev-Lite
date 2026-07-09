"""Sequence-level antibody variable-region liability scanning."""

from __future__ import annotations

import re

import pandas as pd

from .utils import HYDROPHOBIC, local_window, normalize_text


LIABILITY_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "chain_id",
    "chain_type",
    "region_type",
    "sequence_scope",
    "risk_type",
    "motif",
    "position_start",
    "position_end",
    "local_sequence_window",
    "risk_level",
    "explanation",
]

MOTIF_RULES = [
    ("deamidation", ["NG", "NS", "NT", "NH"], "Medium", "Potential Asn deamidation motif."),
    ("isomerization", ["DG", "DS", "DT"], "Medium", "Potential Asp isomerization motif."),
    ("clipping_risk", ["KK", "KR", "RR", "RK"], "Low", "Basic residue pair may indicate clipping sensitivity."),
    ("acid_sensitive", ["DP"], "Low", "DP motif can be acid-sensitive in some contexts."),
]


def _base_record(row: pd.Series) -> dict[str, object]:
    return {
        "antibody_id": row.get("antibody_id", ""),
        "molecule_format": row.get("molecule_format", ""),
        "chain_id": row.get("chain_id", ""),
        "chain_type": row.get("chain_type", ""),
        "region_type": row.get("region_type", ""),
        "sequence_scope": row.get("sequence_scope", ""),
    }


def _add_site(
    records: list[dict[str, object]],
    row: pd.Series,
    risk_type: str,
    motif: str,
    start: int,
    end: int,
    risk_level: str,
    explanation: str,
) -> None:
    sequence = normalize_text(row.get("cleaned_sequence"))
    record = _base_record(row)
    record.update(
        {
            "risk_type": risk_type,
            "motif": motif,
            "position_start": start,
            "position_end": end,
            "local_sequence_window": local_window(sequence, start, end),
            "risk_level": risk_level,
            "explanation": explanation,
        }
    )
    records.append(record)


def _scan_literal_motifs(records: list[dict[str, object]], row: pd.Series, sequence: str) -> None:
    for risk_type, motifs, risk_level, explanation in MOTIF_RULES:
        for motif in motifs:
            for match in re.finditer(f"(?={motif})", sequence):
                start = match.start() + 1
                _add_site(records, row, risk_type, motif, start, start + len(motif) - 1, risk_level, explanation)


def _scan_regex(records: list[dict[str, object]], row: pd.Series, sequence: str) -> None:
    for match in re.finditer(r"N[^P][ST]", sequence):
        _add_site(
            records,
            row,
            "N_glycosylation",
            match.group(0),
            match.start() + 1,
            match.end(),
            "High",
            "N-X-S/T motif detected where X is not proline.",
        )

    hydro = "".join(sorted(HYDROPHOBIC))
    for match in re.finditer(f"[{hydro}]{{4,}}", sequence):
        _add_site(
            records,
            row,
            "hydrophobic_patch_proxy",
            match.group(0),
            match.start() + 1,
            match.end(),
            "Medium",
            "Four or more consecutive hydrophobic residues detected as a sequence-level proxy.",
        )

    for match in re.finditer(r"(.)\1{4,}", sequence):
        _add_site(
            records,
            row,
            "low_complexity",
            match.group(0),
            match.start() + 1,
            match.end(),
            "Medium",
            "Same amino acid repeated five or more times.",
        )

    for match in re.finditer(r"([A-Z]{2})\1{3,}", sequence):
        _add_site(
            records,
            row,
            "low_complexity",
            match.group(0),
            match.start() + 1,
            match.end(),
            "Medium",
            "Dipeptide repeat detected four or more times.",
        )


def _scan_residue_risks(records: list[dict[str, object]], row: pd.Series, sequence: str) -> None:
    for index, aa in enumerate(sequence, start=1):
        if aa == "M":
            _add_site(records, row, "oxidation", "M", index, index, "Medium", "Methionine oxidation liability.")
        elif aa == "W":
            _add_site(records, row, "oxidation", "W", index, index, "High", "Tryptophan oxidation liability.")

    cys_count = sequence.count("C")
    region_type = normalize_text(row.get("region_type")).lower()
    if cys_count % 2 == 1 and cys_count > 0:
        _add_site(
            records,
            row,
            "cysteine_risk",
            "potential_unpaired_cysteine",
            1,
            max(len(sequence), 1),
            "High",
            "Odd cysteine count is a proxy for potential unpaired cysteine.",
        )
    if region_type == "variable_region" and cys_count > 4:
        _add_site(
            records,
            row,
            "cysteine_risk",
            "high_cysteine_count_warning",
            1,
            max(len(sequence), 1),
            "Medium",
            "More than four cysteines detected in a variable-region sequence.",
        )


def scan_liabilities(df_cleaned_or_qc: pd.DataFrame) -> pd.DataFrame:
    """Scan all cleaned sequences for MVP liability motifs."""
    records: list[dict[str, object]] = []
    for _, row in df_cleaned_or_qc.iterrows():
        sequence = normalize_text(row.get("cleaned_sequence"))
        if not sequence:
            continue
        _scan_literal_motifs(records, row, sequence)
        _scan_regex(records, row, sequence)
        _scan_residue_risks(records, row, sequence)

    return pd.DataFrame(records, columns=LIABILITY_COLUMNS)
