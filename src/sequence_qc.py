"""Basic sequence quality-control checks."""

from __future__ import annotations

import re

import pandas as pd

from .utils import STOP_GAP_CHARS, join_messages, normalize_text


HOMOPOLYMER_PATTERN = re.compile(r"(.)\1{4,}")
DIPEPTIDE_REPEAT_PATTERN = re.compile(r"([A-Z]{2})\1{3,}")


def _length_warning(row: pd.Series, length: int) -> str:
    region_type = normalize_text(row.get("region_type")).lower()
    scope = normalize_text(row.get("sequence_scope")).lower()
    if scope == "scfv":
        if length < 180 or length > 350:
            return "scFv length is outside the MVP suggested range of 180-350 aa"
        return ""
    if region_type == "variable_region" or region_type == "":
        if length < 70:
            return "variable_region sequence is shorter than 70 aa"
        if length > 180:
            return "variable_region sequence is longer than 180 aa"
    return ""


def qc_row(row: pd.Series) -> dict[str, object]:
    """Run sequence QC for a single row."""
    sequence = normalize_text(row.get("cleaned_sequence"))
    original = normalize_text(row.get("original_sequence")).upper()
    illegal = normalize_text(row.get("illegal_characters"))
    length = len(sequence)
    cys_count = sequence.count("C")
    warnings: list[str] = []
    fail = False

    if not sequence:
        warnings.append("empty cleaned sequence")
        fail = True
    if illegal:
        warnings.append(f"illegal character(s) detected: {illegal}")
    if any(char in original for char in STOP_GAP_CHARS):
        warnings.append("stop-like or gap-like character detected in original sequence")

    length_warning = _length_warning(row, length)
    if length_warning:
        warnings.append(length_warning)

    region_type = normalize_text(row.get("region_type")).lower()
    if region_type == "variable_region" and length and (length < 80 or length > 160):
        warnings.append("sequence is outside typical VH/VL/VHH variable-region range of 80-160 aa")

    if cys_count == 0 and length:
        warnings.append("no cysteine detected; verify variable-region input")
    elif cys_count % 2 == 1:
        warnings.append("odd cysteine count may indicate unpaired cysteine")
    elif cys_count > 4 and region_type == "variable_region":
        warnings.append("high cysteine count for a variable-region sequence")

    low_complexity_warning = ""
    if HOMOPOLYMER_PATTERN.search(sequence) or DIPEPTIDE_REPEAT_PATTERN.search(sequence):
        low_complexity_warning = "low-complexity repeat detected"
        warnings.append(low_complexity_warning)

    qc_status = "FAIL" if fail else "WARNING" if warnings else "PASS"
    return {
        "qc_status": qc_status,
        "qc_warnings": join_messages(warnings),
        "sequence_length": length,
        "illegal_characters": illegal,
        "cys_count": cys_count,
        "low_complexity_warning": low_complexity_warning,
        "variable_region_length_warning": length_warning,
    }


def run_qc(classified_df: pd.DataFrame) -> pd.DataFrame:
    """Run QC checks for all rows."""
    qc = classified_df.apply(qc_row, axis=1, result_type="expand")
    base = classified_df[
        [
            "antibody_id",
            "molecule_format",
            "chain_id",
            "chain_type",
            "region_type",
            "sequence_scope",
            "cleaned_sequence",
        ]
    ].reset_index(drop=True)
    return pd.concat([base, qc.reset_index(drop=True)], axis=1)
