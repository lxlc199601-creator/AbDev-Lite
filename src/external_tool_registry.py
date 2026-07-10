"""External tool registry for manual-assisted adapter workflows."""

from __future__ import annotations

from pathlib import Path
import re

import pandas as pd


REGISTRY_COLUMNS = [
    "tool_id",
    "tool_name",
    "tool_category",
    "input_type",
    "output_type",
    "recommended_use",
    "requires_web_upload",
    "supports_batch",
    "supports_api",
    "supports_cli",
    "default_mode",
    "automation_allowed",
    "notes",
    "expected_result_fields",
]

TOOL_REGISTRY_COLUMNS = [
    "tool_id",
    "tool_name",
    "tool_category",
    "default_mode",
    "automation_allowed",
    "supports_batch",
    "supports_api",
    "supports_cli",
    "requires_web_upload",
    "recommended_use",
    "notes",
]

VALID_CATEGORIES = {"humanness", "germline", "developability", "structure", "formulation", "custom"}
VALID_MODES = {"manual_export_import", "local_cli", "api", "browser_automation_disabled"}

DEFAULT_TOOL_SPECS = [
    {
        "tool_id": "BioPhi_OASis",
        "tool_name": "BioPhi / OASis",
        "tool_category": "humanness",
        "input_type": "FASTA",
        "output_type": "CSV/XLSX",
        "recommended_use": "Human-likeness and OASis-style humanness review.",
        "requires_web_upload": True,
        "supports_batch": True,
        "supports_api": False,
        "supports_cli": False,
        "default_mode": "manual_export_import",
        "automation_allowed": False,
        "notes": "Prepare FASTA locally, run manually according to tool terms, then import the result.",
        "expected_result_fields": "antibody_id, chain_id, sequence_scope, result_risk_class, result_metric_name, result_metric_value",
    },
    {
        "tool_id": "IgBLAST",
        "tool_name": "IgBLAST",
        "tool_category": "germline",
        "input_type": "FASTA",
        "output_type": "CSV/TSV/XLSX",
        "recommended_use": "Germline assignment and V(D)J annotation support.",
        "requires_web_upload": False,
        "supports_batch": True,
        "supports_api": False,
        "supports_cli": True,
        "default_mode": "local_cli",
        "automation_allowed": False,
        "notes": "Local CLI use is preferred when available; no browser automation is enabled by default.",
        "expected_result_fields": "antibody_id, chain_id, sequence_scope, result_metric_name, result_metric_value",
    },
    {
        "tool_id": "TAP",
        "tool_name": "TAP",
        "tool_category": "developability",
        "input_type": "FASTA",
        "output_type": "CSV/XLSX",
        "recommended_use": "External developability flag review.",
        "requires_web_upload": True,
        "supports_batch": True,
        "supports_api": False,
        "supports_cli": False,
        "default_mode": "manual_export_import",
        "automation_allowed": False,
        "notes": "User-managed export/import workflow only in v0.9.",
        "expected_result_fields": "antibody_id, chain_id, sequence_scope, result_risk_class, result_interpretation",
    },
    {
        "tool_id": "IgFold",
        "tool_name": "IgFold",
        "tool_category": "structure",
        "input_type": "FASTA",
        "output_type": "CSV/XLSX/PDB",
        "recommended_use": "Optional structure prediction evidence generated outside AbDev-Lite.",
        "requires_web_upload": False,
        "supports_batch": True,
        "supports_api": False,
        "supports_cli": True,
        "default_mode": "local_cli",
        "automation_allowed": False,
        "notes": "Import only summary-level user-provided structure evidence in this framework.",
        "expected_result_fields": "antibody_id, chain_id, sequence_scope, result_metric_name, result_metric_value, result_risk_class",
    },
    {
        "tool_id": "ImmuneBuilder",
        "tool_name": "ImmuneBuilder",
        "tool_category": "structure",
        "input_type": "FASTA",
        "output_type": "CSV/XLSX/PDB",
        "recommended_use": "Optional antibody structure model evidence generated outside AbDev-Lite.",
        "requires_web_upload": False,
        "supports_batch": True,
        "supports_api": False,
        "supports_cli": True,
        "default_mode": "local_cli",
        "automation_allowed": False,
        "notes": "Browser automation and credential handling are not implemented.",
        "expected_result_fields": "antibody_id, chain_id, sequence_scope, result_metric_name, result_metric_value, result_risk_class",
    },
    {
        "tool_id": "Custom_Manual_Tool",
        "tool_name": "Custom Manual Tool",
        "tool_category": "custom",
        "input_type": "CSV/FASTA",
        "output_type": "CSV/XLSX",
        "recommended_use": "User-defined external computational evidence.",
        "requires_web_upload": False,
        "supports_batch": True,
        "supports_api": False,
        "supports_cli": False,
        "default_mode": "manual_export_import",
        "automation_allowed": False,
        "notes": "Generic manual export/import path for internal or user-defined tools.",
        "expected_result_fields": "antibody_id, chain_id, sequence_scope, result_risk_class, result_interpretation",
    },
]


def _to_bool(value: object, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "yes", "1"}:
        return True
    if text in {"false", "no", "0"}:
        return False
    return default


def _parse_simple_yaml(path: Path) -> list[dict[str, object]]:
    """Parse the small registry YAML format without requiring PyYAML."""
    records: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip()
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped == "tools:":
            continue
        item_match = re.match(r"^-\s+([^:]+):\s*(.*)$", stripped)
        if item_match:
            if current is not None:
                records.append(current)
            current = {item_match.group(1).strip(): item_match.group(2).strip()}
            continue
        if current is not None and ":" in stripped:
            key, value = stripped.split(":", 1)
            current[key.strip()] = value.strip()
    if current is not None:
        records.append(current)
    return records


def validate_tool_spec(tool_spec: dict) -> dict:
    """Normalize and validate a single external tool specification."""
    spec = {column: tool_spec.get(column, "") for column in REGISTRY_COLUMNS}
    spec["tool_id"] = str(spec["tool_id"]).strip()
    spec["tool_name"] = str(spec["tool_name"]).strip() or spec["tool_id"]
    spec["tool_category"] = str(spec["tool_category"]).strip() or "custom"
    if spec["tool_category"] not in VALID_CATEGORIES:
        spec["tool_category"] = "custom"
    spec["default_mode"] = str(spec["default_mode"]).strip() or "manual_export_import"
    if spec["default_mode"] not in VALID_MODES:
        spec["default_mode"] = "manual_export_import"
    for column in ["automation_allowed", "supports_batch", "supports_api", "supports_cli", "requires_web_upload"]:
        spec[column] = _to_bool(spec[column], False)
    if not spec["tool_id"]:
        raise ValueError("External tool spec is missing tool_id.")
    return spec


def load_tool_registry(tool_specs_path: Path) -> pd.DataFrame:
    """Load external tool specifications, falling back to built-in defaults."""
    specs = DEFAULT_TOOL_SPECS
    path = Path(tool_specs_path)
    if path.exists():
        try:
            parsed = _parse_simple_yaml(path)
            if parsed:
                specs = parsed
        except Exception:
            specs = DEFAULT_TOOL_SPECS
    records = [validate_tool_spec(spec) for spec in specs]
    return pd.DataFrame(records, columns=REGISTRY_COLUMNS)


def get_available_tools() -> pd.DataFrame:
    """Return the registry display table used by reports and the app."""
    registry = load_tool_registry(Path("tool_specs") / "external_tools.yaml")
    for column in TOOL_REGISTRY_COLUMNS:
        if column not in registry.columns:
            registry[column] = ""
    return registry[TOOL_REGISTRY_COLUMNS].copy()
