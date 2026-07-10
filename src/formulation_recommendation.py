"""Rule-based formulation recommendation summaries."""

from __future__ import annotations

import pandas as pd

from .expreso_adapter import EXCIPIENT_CLASSES


FORMULATION_DISCLAIMER = (
    "Formulation recommendations are computational triage outputs based on sequence-level features and optional "
    "Expreso-lite model predictions. They are intended to support early preformulation planning only and do not "
    "replace experimental formulation screening, stability studies, or CMC evaluation."
)

RECOMMENDATION_COLUMNS = [
    "antibody_id",
    "formulation_risk_class",
    "formulation_risk_score",
    "recommended_excipient_classes",
    "top_excipient_class_1",
    "top_excipient_class_2",
    "top_excipient_class_3",
    "buffer_ph_direction",
    "surfactant_consideration",
    "sugar_polyol_consideration",
    "oxidation_control_consideration",
    "ionic_strength_consideration",
    "formulation_review_reason",
    "formulation_next_step_recommendation",
    "formulation_disclaimer",
]

DISPLAY_NAMES = {
    "recommended_buffer_ph_modifier": "buffer / pH modifier",
    "sugar_stabilizer": "sugar stabilizer",
    "polyol_stabilizer": "polyol stabilizer",
    "surfactant": "surfactant",
    "amino_acid_stabilizer": "amino acid stabilizer",
    "antioxidant": "antioxidant",
    "chelating_agent": "chelating agent",
    "salt_ionic_strength_modifier": "salt / ionic strength modifier",
    "preservative": "preservative",
}


def _safe_float(value: object, default: float = 0.0) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return default
    if pd.isna(numeric):
        return default
    return numeric


def _safe_text(value: object, default: str = "") -> str:
    if value is None or pd.isna(value):
        return default
    return str(value)


def _risk_score(row: pd.Series) -> int:
    score = 0
    if _safe_float(row.get("cdr_adjusted_total_risk_score")) > 8:
        score += 3
    if _safe_float(row.get("total_hydrophobic_patch_proxy")) > 0:
        score += 2
    if _safe_float(row.get("total_oxidation_motifs")) >= 3:
        score += 2
    if _safe_float(row.get("total_nglycosylation_motifs")) > 0:
        score += 1
    if _safe_float(row.get("average_instability_index")) > 40:
        score += 2
    if _safe_float(row.get("average_gravy")) > 0:
        score += 2
    pi = _safe_float(row.get("average_theoretical_pI"), 7.0)
    if pi > 9 or pi < 5:
        score += 1
    priority_class = _safe_text(row.get("final_priority_class")).upper()
    if priority_class == "C":
        score += 1
    elif priority_class == "D":
        score += 2
    return score


def _risk_class(score: int) -> str:
    if score <= 2:
        return "Low"
    if score <= 5:
        return "Medium"
    return "High"


def _next_step(risk_class: str) -> str:
    return {
        "Low": "Proceed with standard preformulation screening",
        "Medium": "Perform targeted buffer, stabilizer, and surfactant screening",
        "High": "Prioritize preformulation risk review before scale-up or advanced developability studies",
    }[risk_class]


def _top_classes(prediction_row: pd.Series | None) -> list[str]:
    scored: list[tuple[str, float]] = []
    if prediction_row is not None:
        for class_name in EXCIPIENT_CLASSES:
            scored.append((class_name, _safe_float(prediction_row.get(f"{class_name}_probability"))))
    if not scored or max(score for _, score in scored) == 0:
        scored = [
            ("recommended_buffer_ph_modifier", 0.7),
            ("sugar_stabilizer", 0.6),
            ("surfactant", 0.6),
        ]
    scored.sort(key=lambda item: item[1], reverse=True)
    ordered = [DISPLAY_NAMES[name] for name, _ in scored]
    defaults = ["buffer / pH modifier", "sugar stabilizer", "surfactant"]
    for item in defaults:
        if item not in ordered:
            ordered.append(item)
    return ordered[:3]


def _review_reasons(row: pd.Series, risk_class: str) -> str:
    reasons: list[str] = []
    if _safe_float(row.get("total_oxidation_motifs")) >= 3:
        reasons.append("Elevated oxidation-prone residue count")
    if _safe_float(row.get("hydrophobic_residue_fraction_mean")) > 0.42 or _safe_float(row.get("total_hydrophobic_patch_proxy")) > 0:
        reasons.append("Hydrophobicity or hydrophobic patch proxy suggests surfactant review")
    if _safe_float(row.get("cdr_adjusted_total_risk_score")) > 8:
        reasons.append("High CDR-adjusted developability burden")
    pi = _safe_float(row.get("average_theoretical_pI"), 7.0)
    if pi > 9 or pi < 5:
        reasons.append("Theoretical pI is outside the broad neutral range")
    if _safe_float(row.get("average_instability_index")) > 40:
        reasons.append("Instability index suggests stabilizer screening")
    if _safe_float(row.get("total_deamidation_motifs")) + _safe_float(row.get("total_isomerization_motifs")) > 0:
        reasons.append("Deamidation/isomerization motif burden suggests pH sensitivity review")
    cys = _safe_float(row.get("total_cysteine_count"))
    chains = max(1.0, _safe_float(row.get("number_of_chains"), 1.0))
    if cys % 2 == 1 or cys / chains > 4:
        reasons.append("Cysteine count suggests redox-related review")
    if not reasons:
        reasons.append(f"{risk_class} formulation risk by sequence-level screening rules")
    return "; ".join(reasons)


def build_formulation_recommendations(
    formulation_features_df: pd.DataFrame,
    expreso_predictions_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build early-stage formulation recommendation rows."""
    features = formulation_features_df if isinstance(formulation_features_df, pd.DataFrame) else pd.DataFrame()
    predictions = expreso_predictions_df if isinstance(expreso_predictions_df, pd.DataFrame) else pd.DataFrame()
    if features.empty or "antibody_id" not in features.columns:
        return pd.DataFrame(columns=RECOMMENDATION_COLUMNS)

    records: list[dict[str, object]] = []
    for _, row in features.iterrows():
        antibody_id = row.get("antibody_id", "")
        prediction_row = None
        if not predictions.empty and "antibody_id" in predictions.columns:
            matches = predictions[predictions["antibody_id"] == antibody_id]
            if not matches.empty:
                prediction_row = matches.iloc[0]
        score = _risk_score(row)
        risk_class = _risk_class(score)
        top_classes = _top_classes(prediction_row)
        pi = _safe_float(row.get("average_theoretical_pI"), 7.0)
        buffer_direction = "Include broad buffer and pH screening"
        if pi > 9:
            buffer_direction = "Prioritize pH screening below the high theoretical pI range"
        elif pi < 5:
            buffer_direction = "Prioritize pH screening above the low theoretical pI range"
        if _safe_float(row.get("total_deamidation_motifs")) + _safe_float(row.get("total_isomerization_motifs")) > 0:
            buffer_direction += "; review pH-sensitive motifs"

        surfactant = "Consider standard surfactant screening"
        if _safe_float(row.get("hydrophobic_residue_fraction_mean")) > 0.42 or _safe_float(row.get("total_hydrophobic_patch_proxy")) > 0:
            surfactant = "Prioritize surfactant screening due to hydrophobicity signal"

        sugar_polyol = "Consider sugar stabilizer as a default early screen"
        if _safe_float(row.get("average_instability_index")) > 40:
            sugar_polyol = "Prioritize sugar/polyol stabilizer screening due to instability signal"

        oxidation = "Routine oxidation control review"
        if _safe_float(row.get("total_oxidation_motifs")) >= 3:
            oxidation = "Consider antioxidant and oxidation-control review"

        ionic = "Optional ionic strength screen"
        if pi > 9 or pi < 5:
            ionic = "Include ionic strength screen because theoretical pI is shifted"

        records.append(
            {
                "antibody_id": antibody_id,
                "formulation_risk_class": risk_class,
                "formulation_risk_score": score,
                "recommended_excipient_classes": "; ".join(top_classes),
                "top_excipient_class_1": top_classes[0],
                "top_excipient_class_2": top_classes[1],
                "top_excipient_class_3": top_classes[2],
                "buffer_ph_direction": buffer_direction,
                "surfactant_consideration": surfactant,
                "sugar_polyol_consideration": sugar_polyol,
                "oxidation_control_consideration": oxidation,
                "ionic_strength_consideration": ionic,
                "formulation_review_reason": _review_reasons(row, risk_class),
                "formulation_next_step_recommendation": _next_step(risk_class),
                "formulation_disclaimer": FORMULATION_DISCLAIMER,
            }
        )

    return pd.DataFrame(records, columns=RECOMMENDATION_COLUMNS)
