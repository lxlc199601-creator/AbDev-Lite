"""Biopython-based physicochemical property calculation."""

from __future__ import annotations

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from .utils import HYDROPHOBIC, NEGATIVE, POLAR, POSITIVE, VALID_AA, normalize_text


def _fraction(sequence: str, residues: set[str]) -> float | None:
    if not sequence:
        return None
    return sum(1 for aa in sequence if aa in residues) / len(sequence)


def calculate_properties(row: pd.Series) -> dict[str, object]:
    """Calculate sequence properties for one cleaned sequence."""
    sequence = normalize_text(row.get("cleaned_sequence"))
    invalid = sorted({aa for aa in sequence if aa not in VALID_AA})
    if not sequence or invalid:
        return {
            "sequence_length": len(sequence),
            "property_warning": "properties not calculated because sequence is empty or invalid",
        }

    try:
        analysis = ProteinAnalysis(sequence)
        return {
            "sequence_length": len(sequence),
            "molecular_weight": round(analysis.molecular_weight(), 3),
            "theoretical_pI": round(analysis.isoelectric_point(), 3),
            "aromaticity": round(analysis.aromaticity(), 4),
            "instability_index": round(analysis.instability_index(), 3),
            "gravy": round(analysis.gravy(), 4),
            "hydrophobic_residue_fraction": round(_fraction(sequence, HYDROPHOBIC), 4),
            "aromatic_residue_fraction": round(_fraction(sequence, set("FWY")), 4),
            "positive_residue_fraction": round(_fraction(sequence, POSITIVE), 4),
            "negative_residue_fraction": round(_fraction(sequence, NEGATIVE), 4),
            "polar_residue_fraction": round(_fraction(sequence, POLAR), 4),
            "cysteine_count": sequence.count("C"),
            "methionine_count": sequence.count("M"),
            "tryptophan_count": sequence.count("W"),
            "asparagine_count": sequence.count("N"),
            "aspartic_acid_count": sequence.count("D"),
            "lysine_count": sequence.count("K"),
            "arginine_count": sequence.count("R"),
            "property_warning": "",
        }
    except Exception as exc:
        return {
            "sequence_length": len(sequence),
            "property_warning": f"Biopython property calculation failed: {exc}",
        }


def calculate_properties_table(classified_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate properties for every chain."""
    properties = classified_df.apply(calculate_properties, axis=1, result_type="expand")
    base = classified_df[
        ["antibody_id", "molecule_format", "chain_id", "chain_type", "region_type", "sequence_scope"]
    ].reset_index(drop=True)
    return pd.concat([base, properties.reset_index(drop=True)], axis=1)
