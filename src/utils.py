"""Shared constants and helper functions for AbDev-Lite."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


REQUIRED_COLUMNS = [
    "antibody_id",
    "molecule_format",
    "chain_id",
    "chain_type",
    "region_type",
    "sequence_scope",
    "sequence",
]

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
STOP_GAP_CHARS = set("*-.")
HYDROPHOBIC = set("AVILMFWY")
POSITIVE = set("KRH")
NEGATIVE = set("DE")
POLAR = set("STNQCY")

RISK_CLASS_RANK = {"Low": 0, "Medium": 1, "High": 2}


def ensure_output_dir(path: str | Path = "outputs") -> Path:
    """Create and return the output directory."""
    output_dir = Path(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def join_messages(messages: Iterable[str]) -> str:
    """Join warning or status messages for table display."""
    clean = [str(message).strip() for message in messages if str(message).strip()]
    return "; ".join(dict.fromkeys(clean))


def normalize_text(value: object) -> str:
    """Return a stripped string while keeping missing values empty."""
    if value is None:
        return ""
    text = str(value)
    if text.lower() == "nan":
        return ""
    return text.strip()


def local_window(sequence: str, start: int, end: int) -> str:
    """Return a 1-based motif window with the motif highlighted in brackets."""
    left = max(start - 6, 0)
    right = min(end + 5, len(sequence))
    prefix = sequence[left : start - 1]
    motif = sequence[start - 1 : end]
    suffix = sequence[end:right]
    return f"{prefix}[{motif}]{suffix}"


def risk_class_from_score(score: float) -> str:
    """Convert a numeric chain score into a risk class."""
    if score <= 2:
        return "Low"
    if score <= 5:
        return "Medium"
    return "High"


def max_risk_class(classes: Iterable[str]) -> str:
    """Return the highest risk class in an iterable."""
    clean = [risk_class for risk_class in classes if risk_class in RISK_CLASS_RANK]
    if not clean:
        return "Low"
    return max(clean, key=lambda item: RISK_CLASS_RANK[item])
