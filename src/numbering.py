"""Optional IMGT numbering and CDR/FR region mapping."""

from __future__ import annotations

import importlib
import re
from collections.abc import Iterable
from typing import Any

import pandas as pd

from .utils import join_messages, normalize_text


NUMBERING_SCHEME = "IMGT"

NUMBERING_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "chain_id",
    "chain_type",
    "region_type",
    "sequence_scope",
    "original_position",
    "residue",
    "numbering_scheme",
    "imgt_position",
    "insertion_code",
    "region_label",
    "numbering_status",
    "numbering_warning",
]

REGION_SUMMARY_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "chain_id",
    "chain_type",
    "region_type",
    "sequence_scope",
    "numbering_scheme",
    "numbering_status",
    "numbering_warning",
    "inferred_chain_type",
    "FR1_sequence",
    "CDR1_sequence",
    "FR2_sequence",
    "CDR2_sequence",
    "FR3_sequence",
    "CDR3_sequence",
    "FR4_sequence",
    "CDR1_length",
    "CDR2_length",
    "CDR3_length",
    "total_numbered_residues",
]

LIABILITY_REGION_MAP_COLUMNS = [
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
    "mapped_region_start",
    "mapped_region_end",
    "mapped_regions",
    "cdr_or_fr",
    "imgt_positions",
    "region_mapping_status",
    "region_mapping_warning",
]

REGION_ORDER = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
CDR_REGIONS = {"CDR1", "CDR2", "CDR3"}
FR_REGIONS = {"FR1", "FR2", "FR3", "FR4"}


def _base_chain_record(row: pd.Series) -> dict[str, object]:
    return {
        "antibody_id": row.get("antibody_id", ""),
        "molecule_format": row.get("molecule_format", ""),
        "chain_id": row.get("chain_id", ""),
        "chain_type": row.get("chain_type", ""),
        "region_type": row.get("region_type", ""),
        "sequence_scope": row.get("sequence_scope", ""),
    }


def _empty_region_sequences() -> dict[str, object]:
    return {
        "FR1_sequence": "",
        "CDR1_sequence": "",
        "FR2_sequence": "",
        "CDR2_sequence": "",
        "FR3_sequence": "",
        "CDR3_sequence": "",
        "FR4_sequence": "",
        "CDR1_length": 0,
        "CDR2_length": 0,
        "CDR3_length": 0,
    }


def _summary_record(
    row: pd.Series,
    status: str,
    warning: str,
    inferred_chain_type: str = "unknown",
    region_sequences: dict[str, object] | None = None,
    total_numbered_residues: int = 0,
) -> dict[str, object]:
    record = _base_chain_record(row)
    sequences = _empty_region_sequences()
    if region_sequences:
        sequences.update(region_sequences)
    record.update(
        {
            "numbering_scheme": NUMBERING_SCHEME,
            "numbering_status": status,
            "numbering_warning": warning,
            "inferred_chain_type": inferred_chain_type,
            **sequences,
            "total_numbered_residues": total_numbered_residues,
        }
    )
    return record


def _load_abnumber_chain() -> tuple[Any | None, str]:
    try:
        module = importlib.import_module("abnumber")
        chain_class = getattr(module, "Chain", None)
        if chain_class is None:
            return None, "abnumber is installed but Chain API was not found"
        return chain_class, ""
    except Exception as exc:
        return None, f"AbNumber/ANARCI is not available: {exc}"


def _infer_chain_type(row: pd.Series, chain_obj: Any | None = None) -> str:
    scope = normalize_text(row.get("sequence_scope")).upper()
    provided = normalize_text(row.get("chain_type")).lower()
    api_value = normalize_text(getattr(chain_obj, "chain_type", "") if chain_obj is not None else "").upper()
    if provided == "vhh" or scope == "VHH":
        return "vhh"
    if api_value in {"H", "HEAVY"} or provided == "heavy" or scope.startswith("VH"):
        return "heavy"
    if api_value in {"K", "L", "LIGHT", "KAPPA", "LAMBDA"} or provided == "light" or scope.startswith("VL"):
        return "light"
    return "unknown"


def _is_scfv_scope(row: pd.Series) -> bool:
    scope = normalize_text(row.get("sequence_scope")).lower()
    chain_type = normalize_text(row.get("chain_type")).lower()
    return scope == "scfv" or chain_type == "scfv"


def _iter_numbered_residues(chain_obj: Any) -> list[tuple[Any, str]]:
    positions = getattr(chain_obj, "positions", None)
    if isinstance(positions, dict):
        return [(position, normalize_text(residue)) for position, residue in positions.items() if normalize_text(residue)]

    residues: list[tuple[Any, str]] = []
    try:
        iterator: Iterable[Any] = iter(chain_obj)
    except TypeError:
        iterator = []
    for item in iterator:
        if isinstance(item, tuple) and len(item) >= 2:
            residues.append((item[0], normalize_text(item[1])))
        else:
            position = getattr(item, "position", None) or getattr(item, "pos", None)
            residue = getattr(item, "aa", None) or getattr(item, "residue", None)
            if position is not None and residue is not None:
                residues.append((position, normalize_text(residue)))
    return [(position, residue) for position, residue in residues if residue]


def _position_parts(position: Any) -> tuple[object, str]:
    number = getattr(position, "number", None)
    insertion = getattr(position, "letter", None) or getattr(position, "insertion_code", None) or ""
    if number is not None:
        return number, normalize_text(insertion)

    text = normalize_text(position)
    match = re.search(r"(\d+)([A-Za-z]*)$", text)
    if match:
        return int(match.group(1)), match.group(2)
    return text, ""


def _numeric_position(value: object) -> int | None:
    if isinstance(value, int):
        return value
    match = re.search(r"\d+", normalize_text(value))
    return int(match.group(0)) if match else None


def _region_from_imgt_position(position: object) -> str:
    number = _numeric_position(position)
    if number is None:
        return "Unknown"
    if 1 <= number <= 26:
        return "FR1"
    if 27 <= number <= 38:
        return "CDR1"
    if 39 <= number <= 55:
        return "FR2"
    if 56 <= number <= 65:
        return "CDR2"
    if 66 <= number <= 104:
        return "FR3"
    if 105 <= number <= 117:
        return "CDR3"
    if 118 <= number <= 128:
        return "FR4"
    return "Unknown"


def _region_sequences_from_chain(chain_obj: Any) -> dict[str, object]:
    sequences: dict[str, object] = {}
    for region in REGION_ORDER:
        attr = f"{region.lower()}_seq"
        sequences[f"{region}_sequence"] = normalize_text(getattr(chain_obj, attr, ""))
    sequences["CDR1_length"] = len(normalize_text(sequences["CDR1_sequence"]))
    sequences["CDR2_length"] = len(normalize_text(sequences["CDR2_sequence"]))
    sequences["CDR3_length"] = len(normalize_text(sequences["CDR3_sequence"]))
    return sequences


def _labels_from_region_sequences(region_sequences: dict[str, object], residue_count: int) -> list[str]:
    labels: list[str] = []
    for region in REGION_ORDER:
        labels.extend([region] * len(normalize_text(region_sequences.get(f"{region}_sequence"))))
    if len(labels) == residue_count:
        return labels
    return []


def _original_positions(sequence: str, numbered_residues: list[tuple[Any, str]]) -> list[int]:
    numbered_sequence = "".join(residue for _, residue in numbered_residues)
    start = sequence.find(numbered_sequence)
    if numbered_sequence and start >= 0:
        return list(range(start + 1, start + 1 + len(numbered_sequence)))
    return list(range(1, len(numbered_residues) + 1))


def run_imgt_numbering(qc_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run optional AbNumber/ANARCI IMGT numbering for variable-region rows."""
    residue_records: list[dict[str, object]] = []
    summary_records: list[dict[str, object]] = []
    if qc_df.empty:
        return pd.DataFrame(residue_records, columns=NUMBERING_COLUMNS), pd.DataFrame(
            summary_records, columns=REGION_SUMMARY_COLUMNS
        )

    chain_class, import_warning = _load_abnumber_chain()

    for _, row in qc_df.iterrows():
        region_type = normalize_text(row.get("region_type")).lower()
        sequence = normalize_text(row.get("cleaned_sequence"))
        if region_type != "variable_region":
            summary_records.append(
                _summary_record(
                    row,
                    "SKIPPED",
                    "Numbering module currently supports variable_region sequences only",
                    _infer_chain_type(row),
                )
            )
            continue
        if _is_scfv_scope(row):
            summary_records.append(
                _summary_record(
                    row,
                    "SKIPPED",
                    "scFv sequences are not directly numbered by the v0.4 variable-domain numbering module",
                    "unknown",
                )
            )
            continue
        if not sequence:
            summary_records.append(_summary_record(row, "FAIL", "Empty sequence cannot be numbered", _infer_chain_type(row)))
            continue
        if chain_class is None:
            summary_records.append(_summary_record(row, "FAIL", import_warning, _infer_chain_type(row)))
            continue

        warning_parts: list[str] = []
        try:
            chain_obj = chain_class(sequence, scheme="imgt")
            numbered_residues = _iter_numbered_residues(chain_obj)
            if not numbered_residues:
                raise ValueError("AbNumber returned no numbered residues")

            inferred_chain_type = _infer_chain_type(row, chain_obj)
            region_sequences = _region_sequences_from_chain(chain_obj)
            labels = _labels_from_region_sequences(region_sequences, len(numbered_residues))
            if not labels:
                warning_parts.append("Region labels assigned by fallback IMGT position ranges")

            original_positions = _original_positions(sequence, numbered_residues)
            for index, ((position, residue), original_position) in enumerate(zip(numbered_residues, original_positions)):
                imgt_position, insertion_code = _position_parts(position)
                region_label = labels[index] if labels else _region_from_imgt_position(imgt_position)
                residue_record = _base_chain_record(row)
                residue_record.update(
                    {
                        "original_position": original_position,
                        "residue": residue,
                        "numbering_scheme": NUMBERING_SCHEME,
                        "imgt_position": imgt_position,
                        "insertion_code": insertion_code,
                        "region_label": region_label,
                        "numbering_status": "WARNING" if warning_parts else "PASS",
                        "numbering_warning": join_messages(warning_parts),
                    }
                )
                residue_records.append(residue_record)

            summary_records.append(
                _summary_record(
                    row,
                    "WARNING" if warning_parts else "PASS",
                    join_messages(warning_parts),
                    inferred_chain_type,
                    region_sequences,
                    len(numbered_residues),
                )
            )
        except Exception as exc:
            summary_records.append(_summary_record(row, "FAIL", f"IMGT numbering failed: {exc}", _infer_chain_type(row)))

    return pd.DataFrame(residue_records, columns=NUMBERING_COLUMNS), pd.DataFrame(
        summary_records, columns=REGION_SUMMARY_COLUMNS
    )


def _key_columns() -> list[str]:
    return ["antibody_id", "chain_id", "sequence_scope"]


def _region_category(regions: list[str]) -> str:
    clean = {region for region in regions if region and region != "Unknown"}
    if clean and clean.issubset(CDR_REGIONS):
        return "CDR"
    if clean and clean.issubset(FR_REGIONS):
        return "FR"
    if clean & CDR_REGIONS and clean & FR_REGIONS:
        return "Boundary"
    return "Unknown"


def _format_imgt_positions(rows: pd.DataFrame) -> str:
    values: list[str] = []
    for _, row in rows.iterrows():
        pos = normalize_text(row.get("imgt_position"))
        insertion = normalize_text(row.get("insertion_code"))
        if pos:
            values.append(f"{pos}{insertion}")
    return ";".join(values)


def map_liabilities_to_regions(liability_df: pd.DataFrame, numbering_df: pd.DataFrame) -> pd.DataFrame:
    """Map liability motif positions onto numbered IMGT CDR/FR regions."""
    records: list[dict[str, object]] = []
    if liability_df.empty:
        return pd.DataFrame(records, columns=LIABILITY_REGION_MAP_COLUMNS)

    numbering_lookup: dict[tuple[object, object, object], pd.DataFrame] = {}
    if not numbering_df.empty:
        for key, group in numbering_df.groupby(_key_columns(), dropna=False):
            numbering_lookup[key] = group

    for _, site in liability_df.iterrows():
        record = {column: site.get(column, "") for column in LIABILITY_BASE_COLUMNS}
        start = int(site.get("position_start", 0) or 0)
        end = int(site.get("position_end", 0) or 0)
        key = tuple(site.get(column, "") for column in _key_columns())
        chain_numbering = numbering_lookup.get(key, pd.DataFrame())

        if chain_numbering.empty:
            record.update(
                {
                    "mapped_region_start": "",
                    "mapped_region_end": "",
                    "mapped_regions": "Unknown",
                    "cdr_or_fr": "Unknown",
                    "imgt_positions": "",
                    "region_mapping_status": "FAIL",
                    "region_mapping_warning": "No successful numbering residues available for this chain",
                }
            )
            records.append(record)
            continue

        mapped = chain_numbering[
            chain_numbering["original_position"].astype(int).between(start, end, inclusive="both")
        ].copy()
        if mapped.empty:
            record.update(
                {
                    "mapped_region_start": "",
                    "mapped_region_end": "",
                    "mapped_regions": "Unknown",
                    "cdr_or_fr": "Unknown",
                    "imgt_positions": "",
                    "region_mapping_status": "WARNING",
                    "region_mapping_warning": "Liability position was not found in numbered residues",
                }
            )
            records.append(record)
            continue

        regions = [region for region in mapped["region_label"].astype(str).tolist() if region]
        unique_regions = list(dict.fromkeys(regions)) or ["Unknown"]
        record.update(
            {
                "mapped_region_start": unique_regions[0],
                "mapped_region_end": unique_regions[-1],
                "mapped_regions": ";".join(unique_regions),
                "cdr_or_fr": _region_category(unique_regions),
                "imgt_positions": _format_imgt_positions(mapped),
                "region_mapping_status": "PASS" if "Unknown" not in unique_regions else "WARNING",
                "region_mapping_warning": "" if "Unknown" not in unique_regions else "At least one residue mapped to Unknown region",
            }
        )
        records.append(record)

    return pd.DataFrame(records, columns=LIABILITY_REGION_MAP_COLUMNS)


LIABILITY_BASE_COLUMNS = [
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
