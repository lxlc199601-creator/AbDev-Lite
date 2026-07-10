"""Build auditable input packages for external antibody analysis tools."""

from __future__ import annotations

from datetime import datetime
import json
from pathlib import Path
import re

import pandas as pd

from .external_tool_registry import load_tool_registry
from .utils import ensure_output_dir, normalize_text


RUN_PLAN_COLUMNS = [
    "run_id",
    "tool_id",
    "tool_name",
    "antibody_id",
    "chain_id",
    "sequence_scope",
    "input_file",
    "input_format",
    "run_mode",
    "run_status",
    "created_at",
    "user_action_required",
    "run_instruction",
    "warning",
]

RUN_INSTRUCTION = (
    "Exported input file generated. Please run this file manually with the selected external tool and upload the "
    "result file back to AbDev-Lite."
)


def _safe_name(value: object) -> str:
    text = normalize_text(value) or "unknown"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "unknown"


def _sequence_rows(input_cleaned_df: pd.DataFrame, qc_df: pd.DataFrame | None = None) -> pd.DataFrame:
    rows = input_cleaned_df.copy() if isinstance(input_cleaned_df, pd.DataFrame) else pd.DataFrame()
    if rows.empty:
        return pd.DataFrame(columns=["antibody_id", "chain_id", "sequence_scope", "sequence", "qc_status"])
    for column in ["antibody_id", "chain_id", "sequence_scope", "sequence"]:
        if column == "sequence" and column not in rows.columns and "cleaned_sequence" in rows.columns:
            rows[column] = rows["cleaned_sequence"]
        elif column not in rows.columns:
            rows[column] = ""
    rows["qc_status"] = ""
    if isinstance(qc_df, pd.DataFrame) and not qc_df.empty:
        keys = ["antibody_id", "chain_id", "sequence_scope"]
        available = [column for column in keys + ["qc_status"] if column in qc_df.columns]
        if set(keys).issubset(available) and "qc_status" in available:
            rows = rows.drop(columns=["qc_status"], errors="ignore").merge(qc_df[available], on=keys, how="left")
    return rows


def _write_fasta(rows: pd.DataFrame, output_path: Path) -> None:
    lines: list[str] = []
    for _, row in rows.iterrows():
        header = "|".join(
            [
                _safe_name(row.get("antibody_id")),
                _safe_name(row.get("chain_id")),
                _safe_name(row.get("sequence_scope")),
            ]
        )
        sequence = normalize_text(row.get("sequence"))
        lines.extend([f">{header}", sequence])
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _warnings(row: pd.Series) -> str:
    warnings: list[str] = []
    if not normalize_text(row.get("sequence")):
        warnings.append("Sequence is empty")
    if normalize_text(row.get("qc_status")).upper() == "FAIL":
        warnings.append("QC status is FAIL")
    return "; ".join(warnings)


def empty_external_run_plan() -> pd.DataFrame:
    return pd.DataFrame(columns=RUN_PLAN_COLUMNS)


def build_external_input_package(
    input_cleaned_df,
    region_summary_df,
    antibody_summary_df,
    selected_tools: list[str],
    output_dir: Path,
) -> pd.DataFrame:
    """Create one local input package directory per selected external tool."""
    if not selected_tools:
        return empty_external_run_plan()

    registry = load_tool_registry(Path("tool_specs") / "external_tools.yaml")
    registry_by_id = {str(row["tool_id"]): row for _, row in registry.iterrows()}
    package_root = ensure_output_dir(output_dir)
    created_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sequence_rows = _sequence_rows(input_cleaned_df, region_summary_df)
    records: list[dict[str, object]] = []

    for tool_id in selected_tools:
        if tool_id not in registry_by_id:
            continue
        tool = registry_by_id[tool_id]
        safe_tool_id = _safe_name(tool_id)
        tool_dir = ensure_output_dir(package_root / safe_tool_id)
        csv_path = tool_dir / f"{safe_tool_id}_input.csv"
        fasta_path = tool_dir / f"{safe_tool_id}_input.fasta"
        metadata_path = tool_dir / "metadata.json"

        export_columns = ["antibody_id", "chain_id", "sequence_scope", "sequence", "qc_status"]
        for column in export_columns:
            if column not in sequence_rows.columns:
                sequence_rows[column] = ""
        sequence_rows[export_columns].to_csv(csv_path, index=False)
        _write_fasta(sequence_rows, fasta_path)
        metadata = {
            "tool_id": tool_id,
            "tool_name": tool.get("tool_name", ""),
            "created_at": created_at,
            "record_count": int(len(sequence_rows)),
            "source": "AbDev-Lite v0.9 external input package",
            "browser_automation_enabled": False,
        }
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

        for _, row in sequence_rows.iterrows():
            for input_file, input_format in [(fasta_path, "FASTA"), (csv_path, "CSV")]:
                records.append(
                    {
                        "run_id": f"{safe_tool_id}_{_safe_name(row.get('antibody_id'))}_{_safe_name(row.get('chain_id'))}_{input_format}",
                        "tool_id": tool_id,
                        "tool_name": tool.get("tool_name", ""),
                        "antibody_id": row.get("antibody_id", ""),
                        "chain_id": row.get("chain_id", ""),
                        "sequence_scope": row.get("sequence_scope", ""),
                        "input_file": str(input_file),
                        "input_format": input_format,
                        "run_mode": tool.get("default_mode", "manual_export_import"),
                        "run_status": "PENDING_MANUAL_RUN",
                        "created_at": created_at,
                        "user_action_required": True,
                        "run_instruction": RUN_INSTRUCTION,
                        "warning": _warnings(row),
                    }
                )

    return pd.DataFrame(records, columns=RUN_PLAN_COLUMNS)
