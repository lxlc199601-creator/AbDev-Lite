"""Lightweight sequence and chain classification rules."""

from __future__ import annotations

import pandas as pd

from .utils import join_messages, normalize_text


HEAVY_SCOPES = {"VH", "VH_ARM1", "VH_ARM2", "FULL_HEAVY"}
LIGHT_SCOPES = {"VL", "VL_ARM1", "VL_ARM2", "FULL_LIGHT"}


def classify_row(row: pd.Series) -> dict[str, str]:
    """Infer a simple sequence family from user-provided metadata."""
    region_type = normalize_text(row.get("region_type")).lower()
    sequence_scope = normalize_text(row.get("sequence_scope"))
    scope_upper = sequence_scope.upper()
    chain_type = normalize_text(row.get("chain_type")).lower()

    warnings: list[str] = []
    inferred_region_type = "full_length" if region_type == "full_length" else "variable_region"

    if scope_upper in HEAVY_SCOPES:
        family = "heavy_variable" if inferred_region_type == "variable_region" else "heavy_full_length"
    elif scope_upper in LIGHT_SCOPES:
        family = "light_variable" if inferred_region_type == "variable_region" else "light_full_length"
    elif scope_upper == "VHH":
        family = "vhh"
    elif scope_upper == "SCFV":
        family = "scfv"
    else:
        family = "unknown"
        warnings.append("sequence_scope is not recognized by MVP rules")

    if region_type == "full_length":
        warnings.append(
            "Current version focuses on variable-region screening; full-length analysis is basic reference only"
        )
    elif region_type and region_type != "variable_region":
        warnings.append("region_type is not standard; treated as variable_region for MVP analysis")

    if family.startswith("heavy") and chain_type not in {"heavy", "unknown", ""}:
        warnings.append("chain_type may not match heavy sequence_scope")
    if family.startswith("light") and chain_type not in {"light", "unknown", ""}:
        warnings.append("chain_type may not match light sequence_scope")
    if family == "vhh" and chain_type not in {"vhh", "heavy", "unknown", ""}:
        warnings.append("chain_type may not match VHH sequence_scope")
    if family == "scfv" and chain_type not in {"scfv", "unknown", ""}:
        warnings.append("chain_type may not match scFv sequence_scope")

    return {
        "inferred_region_type": inferred_region_type,
        "inferred_sequence_scope": sequence_scope or "unknown",
        "inferred_chain_family": family,
        "classification_warning": join_messages(warnings),
    }


def classify_sequences(input_df: pd.DataFrame) -> pd.DataFrame:
    """Classify all input rows with lightweight rules."""
    classified = input_df.apply(classify_row, axis=1, result_type="expand")
    return pd.concat([input_df.reset_index(drop=True), classified.reset_index(drop=True)], axis=1)
