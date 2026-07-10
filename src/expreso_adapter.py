"""Optional Expreso-lite model adapter."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


EXCIPIENT_CLASSES = [
    "recommended_buffer_ph_modifier",
    "sugar_stabilizer",
    "polyol_stabilizer",
    "surfactant",
    "amino_acid_stabilizer",
    "antioxidant",
    "chelating_agent",
    "salt_ionic_strength_modifier",
    "preservative",
]

PREDICTION_COLUMNS = [
    "antibody_id",
    "model_available",
    "model_status",
    "prediction_mode",
    *[f"{name}_probability" for name in EXCIPIENT_CLASSES],
    "expreso_prediction_warning",
]

MODEL_NOT_FOUND_WARNING = (
    "Expreso-lite model files were not found. Rule-based fallback recommendation was used."
)


def _safe_float(value: object, default: float = 0.0) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return default
    if pd.isna(numeric):
        return default
    return numeric


def _rough_scores(row: pd.Series) -> dict[str, float]:
    oxidation = _safe_float(row.get("total_oxidation_motifs"))
    hydrophobic = _safe_float(row.get("hydrophobic_residue_fraction_mean"))
    hydrophobic_patch = _safe_float(row.get("total_hydrophobic_patch_proxy"))
    instability = _safe_float(row.get("average_instability_index"))
    gravy = _safe_float(row.get("average_gravy"))
    pi = _safe_float(row.get("average_theoretical_pI"), 7.0)
    deamidation = _safe_float(row.get("total_deamidation_motifs"))
    isomerization = _safe_float(row.get("total_isomerization_motifs"))
    risk = _safe_float(row.get("cdr_adjusted_total_risk_score"))

    return {
        "recommended_buffer_ph_modifier": min(0.95, 0.55 + (0.2 if pi > 9 or pi < 5 else 0.0) + (0.1 if deamidation + isomerization > 0 else 0.0)),
        "sugar_stabilizer": min(0.9, 0.55 + (0.2 if instability > 40 else 0.0) + (0.1 if risk > 8 else 0.0)),
        "polyol_stabilizer": min(0.85, 0.35 + (0.2 if instability > 40 else 0.0) + (0.1 if gravy > 0 else 0.0)),
        "surfactant": min(0.95, 0.55 + (0.2 if hydrophobic > 0.42 else 0.0) + (0.2 if hydrophobic_patch > 0 else 0.0)),
        "amino_acid_stabilizer": min(0.8, 0.3 + (0.15 if pi > 9 or pi < 5 else 0.0) + (0.15 if risk > 8 else 0.0)),
        "antioxidant": min(0.9, 0.25 + (0.15 * min(oxidation, 4))),
        "chelating_agent": min(0.65, 0.2 + (0.1 if oxidation >= 3 else 0.0)),
        "salt_ionic_strength_modifier": min(0.8, 0.3 + (0.15 if pi > 9 or pi < 5 else 0.0)),
        "preservative": 0.15,
    }


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_expreso_models(model_dir: Path) -> dict:
    """Load optional Expreso-lite model metadata and joblib estimators."""
    model_dir = Path(model_dir)
    result: dict[str, Any] = {
        "model_available": False,
        "model_status": "model_not_found",
        "feature_columns": [],
        "manifest": {},
        "models": {},
        "warnings": [],
    }
    feature_path = model_dir / "feature_columns.json"
    if not model_dir.exists() or not feature_path.exists():
        result["warnings"].append(MODEL_NOT_FOUND_WARNING)
        return result

    try:
        feature_columns = _load_json(feature_path)
        if not isinstance(feature_columns, list):
            raise ValueError("feature_columns.json must contain a JSON list")
        result["feature_columns"] = [str(column) for column in feature_columns]
    except Exception as exc:
        result["model_status"] = "feature_columns_load_failed"
        result["warnings"].append(f"feature_columns.json could not be loaded: {exc}")
        return result

    manifest_path = model_dir / "model_manifest.json"
    if manifest_path.exists():
        try:
            manifest = _load_json(manifest_path)
            result["manifest"] = manifest if isinstance(manifest, dict) else {}
        except Exception as exc:
            result["warnings"].append(f"model_manifest.json could not be loaded: {exc}")

    joblib_files = sorted(model_dir.glob("*.joblib"))
    if not joblib_files:
        result["warnings"].append(MODEL_NOT_FOUND_WARNING)
        return result

    try:
        import joblib  # type: ignore
    except Exception as exc:
        result["model_status"] = "joblib_unavailable"
        result["warnings"].append(f"joblib could not be imported: {exc}")
        return result

    for model_path in joblib_files:
        try:
            loaded = joblib.load(model_path)
            if isinstance(loaded, dict):
                for key, value in loaded.items():
                    result["models"][str(key)] = value
            else:
                result["models"][model_path.stem] = loaded
        except Exception as exc:
            result["warnings"].append(f"{model_path.name} failed to load: {exc}")

    if result["models"]:
        result["model_available"] = True
        result["model_status"] = "model_loaded"
    else:
        result["model_status"] = "model_load_failed"
    return result


def _model_for_class(models: dict[str, Any], class_name: str) -> Any | None:
    normalized = class_name.lower()
    candidates = [
        class_name,
        f"{class_name}_probability",
        class_name.replace("_", "-"),
        class_name.replace("_", ""),
    ]
    for candidate in candidates:
        if candidate in models:
            return models[candidate]
    for key, model in models.items():
        if normalized in str(key).lower():
            return model
    return None


def _predict_probability(model: Any, features: pd.DataFrame) -> list[float]:
    if hasattr(model, "predict_proba"):
        values = model.predict_proba(features)
        probabilities = values[:, 1] if getattr(values, "ndim", 1) > 1 and values.shape[1] > 1 else values.ravel()
        return [float(value) for value in probabilities]
    if hasattr(model, "predict"):
        values = model.predict(features)
        return [float(value) for value in values]
    raise TypeError("model does not provide predict_proba or predict")


def predict_with_expreso(
    formulation_features_df: pd.DataFrame,
    model_dir: Path,
) -> pd.DataFrame:
    """Predict excipient-class probabilities or fall back to rule-based scores."""
    features = formulation_features_df if isinstance(formulation_features_df, pd.DataFrame) else pd.DataFrame()
    if features.empty or "antibody_id" not in features.columns:
        return pd.DataFrame(columns=PREDICTION_COLUMNS)

    loaded = load_expreso_models(Path(model_dir))
    warnings = list(loaded.get("warnings", []))

    if not loaded.get("model_available"):
        records = []
        for _, row in features.iterrows():
            scores = _rough_scores(row)
            record = {
                "antibody_id": row.get("antibody_id", ""),
                "model_available": False,
                "model_status": loaded.get("model_status", "model_not_found"),
                "prediction_mode": "rule_based_fallback",
                "expreso_prediction_warning": MODEL_NOT_FOUND_WARNING,
            }
            for class_name in EXCIPIENT_CLASSES:
                record[f"{class_name}_probability"] = round(scores[class_name], 4)
            records.append(record)
        return pd.DataFrame(records, columns=PREDICTION_COLUMNS)

    feature_columns = loaded.get("feature_columns", [])
    aligned = features.copy()
    missing_columns = [column for column in feature_columns if column not in aligned.columns]
    for column in missing_columns:
        aligned[column] = 0
    if missing_columns:
        warnings.append(f"Missing feature columns were filled with 0: {', '.join(missing_columns)}")
    model_input = aligned[feature_columns].apply(pd.to_numeric, errors="coerce").fillna(0)

    predictions_by_class: dict[str, list[float]] = {}
    models = loaded.get("models", {})
    for class_name in EXCIPIENT_CLASSES:
        model = _model_for_class(models, class_name)
        if model is None:
            warnings.append(f"No Expreso-lite model found for {class_name}; fallback score used.")
            predictions_by_class[class_name] = [
                _rough_scores(row)[class_name] for _, row in features.iterrows()
            ]
            continue
        try:
            predictions_by_class[class_name] = _predict_probability(model, model_input)
        except Exception as exc:
            warnings.append(f"Prediction failed for {class_name}: {exc}; fallback score used.")
            predictions_by_class[class_name] = [
                _rough_scores(row)[class_name] for _, row in features.iterrows()
            ]

    warning_text = "; ".join(dict.fromkeys(warnings))
    records = []
    for index, row in features.reset_index(drop=True).iterrows():
        record = {
            "antibody_id": row.get("antibody_id", ""),
            "model_available": True,
            "model_status": loaded.get("model_status", "model_loaded"),
            "prediction_mode": "expreso_lite_model",
            "expreso_prediction_warning": warning_text,
        }
        for class_name in EXCIPIENT_CLASSES:
            values = predictions_by_class.get(class_name, [])
            value = values[index] if index < len(values) else 0.0
            record[f"{class_name}_probability"] = round(float(value), 4)
        records.append(record)
    return pd.DataFrame(records, columns=PREDICTION_COLUMNS)
