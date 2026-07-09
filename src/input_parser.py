"""Input parsing and sequence cleaning."""

from __future__ import annotations

import re
from pathlib import Path
from typing import BinaryIO

import pandas as pd

from .utils import REQUIRED_COLUMNS, VALID_AA, join_messages, normalize_text


REMOVABLE_PATTERN = re.compile(r"[\s\d]+")


def read_input_file(file_obj: BinaryIO | str | Path, filename: str | None = None) -> pd.DataFrame:
    """Read a CSV or XLSX input file into a DataFrame."""
    name = filename or getattr(file_obj, "name", "") or str(file_obj)
    suffix = Path(name).suffix.lower()
    try:
        if suffix == ".csv":
            return pd.read_csv(file_obj)
        if suffix in {".xlsx", ".xls"}:
            return pd.read_excel(file_obj)
    except pd.errors.EmptyDataError as exc:
        raise ValueError("The uploaded file is empty. Please upload a CSV or XLSX file with input rows.") from exc
    except Exception as exc:
        raise ValueError(f"Could not read the uploaded file: {exc}") from exc
    raise ValueError("Unsupported file format. Please upload a CSV or XLSX file.")


def clean_sequence(sequence: object) -> tuple[str, str]:
    """Clean a raw sequence and return cleaned text plus illegal characters."""
    raw = normalize_text(sequence).upper()
    without_spacing_digits = REMOVABLE_PATTERN.sub("", raw)
    illegal = sorted({char for char in without_spacing_digits if char not in VALID_AA})
    cleaned = "".join(char for char in without_spacing_digits if char in VALID_AA)
    return cleaned, "".join(illegal)


def parse_input(file_obj: BinaryIO | str | Path, filename: str | None = None) -> pd.DataFrame:
    """Read, validate, clean, and standardize the input table."""
    df = read_input_file(file_obj, filename)
    if df.empty:
        raise ValueError("The uploaded file has no data rows.")

    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {', '.join(missing)}.")

    df = df.dropna(how="all").copy()
    if df.empty:
        raise ValueError("The uploaded file only contains blank rows.")

    for column in REQUIRED_COLUMNS:
        df[column] = df[column].map(normalize_text)

    df["original_sequence"] = df["sequence"].map(normalize_text)
    cleaned = df["sequence"].map(clean_sequence)
    df["cleaned_sequence"] = cleaned.map(lambda item: item[0])
    df["illegal_characters"] = cleaned.map(lambda item: item[1])

    statuses: list[str] = []
    warnings: list[str] = []
    for _, row in df.iterrows():
        row_warnings: list[str] = []
        if not row["original_sequence"]:
            row_warnings.append("sequence is empty")
        if row["illegal_characters"]:
            row_warnings.append(f"illegal amino-acid character(s): {row['illegal_characters']}")
        if not row["region_type"]:
            row_warnings.append("region_type is missing")
        if not row["sequence_scope"]:
            row_warnings.append("sequence_scope is missing")
        statuses.append("WARNING" if row_warnings else "OK")
        warnings.append(join_messages(row_warnings))

    df["input_status"] = statuses
    df["input_warnings"] = warnings

    return df[
        [
            "antibody_id",
            "molecule_format",
            "chain_id",
            "chain_type",
            "region_type",
            "sequence_scope",
            "original_sequence",
            "cleaned_sequence",
            "input_status",
            "input_warnings",
            "illegal_characters",
        ]
    ]
