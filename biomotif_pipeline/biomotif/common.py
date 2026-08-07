from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import math
import shutil
import tempfile
from contextlib import contextmanager
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

import numpy as np


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def safe_float(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    return val if math.isfinite(val) else None


def safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def min_numeric(values: Iterable[Any]) -> float | None:
    nums = [x for x in (safe_float(v) for v in values) if x is not None]
    return min(nums) if nums else None


def json_default(value: Any) -> Any:
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, set):
        return sorted(value)
    raise TypeError(f"Cannot JSON-serialize {type(value).__name__}")


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True, default=json_default)
    path.write_text(text + "\n", encoding="utf-8")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_jsonl(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(dict(row), sort_keys=True, default=json_default) + "\n")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _csv_value(value: Any) -> Any:
    if isinstance(value, (dict, list, tuple, set)):
        return json.dumps(value, sort_keys=True, default=json_default)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    return value


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sequence_digest(sequences: Sequence[str]) -> str:
    clean = sorted("".join(seq.split()).upper() for seq in sequences if seq)
    if not clean:
        return "none"
    return hashlib.sha256("|".join(clean).encode()).hexdigest()[:24]


def stable_id(*parts: Any, prefix: str = "") -> str:
    raw = "|".join(str(part) for part in parts)
    token = hashlib.sha1(raw.encode()).hexdigest()[:12].upper()
    return f"{prefix}{token}" if prefix else token


@contextmanager
def materialized_text_file(path: Path, suffix: str = "") -> Iterator[Path]:
    if not path.name.lower().endswith(".gz"):
        yield path
        return

    temp = tempfile.NamedTemporaryFile(suffix=suffix or path.stem, delete=False)
    temp_path = Path(temp.name)
    temp.close()
    try:
        with gzip.open(path, "rb") as source, temp_path.open("wb") as target:
            shutil.copyfileobj(source, target)
        yield temp_path
    finally:
        temp_path.unlink(missing_ok=True)


def unit(vector: np.ndarray) -> np.ndarray | None:
    size = float(np.linalg.norm(vector))
    if size < 1e-12:
        return None
    return vector / size


def distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def angle_deg(a: np.ndarray, vertex: np.ndarray, c: np.ndarray) -> float | None:
    left = unit(a - vertex)
    right = unit(c - vertex)
    if left is None or right is None:
        return None
    dot = float(np.clip(np.dot(left, right), -1.0, 1.0))
    return math.degrees(math.acos(dot))


def dihedral_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float | None:
    b0 = a - b
    b1 = c - b
    b2 = d - c
    axis = unit(b1)
    if axis is None:
        return None
    v = unit(b0 - np.dot(b0, axis) * axis)
    w = unit(b2 - np.dot(b2, axis) * axis)
    if v is None or w is None:
        return None
    y = float(np.dot(np.cross(axis, v), w))
    x = float(np.dot(v, w))
    return math.degrees(math.atan2(y, x))


def circular_difference(a: float, b: float) -> float:
    return ((a - b + 180.0) % 360.0) - 180.0


def circular_rms(values: Sequence[float], targets: Sequence[float]) -> float:
    diffs = [circular_difference(a, b) for a, b in zip(values, targets)]
    return math.sqrt(sum(val * val for val in diffs) / len(diffs))


def find_executable(requested: str | None, candidates: Sequence[str]) -> str | None:
    if requested and requested.lower() not in {"auto", "none"}:
        found = shutil.which(requested)
        return found or (requested if Path(requested).exists() else None)
    for name in candidates:
        found = shutil.which(name)
        if found:
            return found
    return None


PROTEIN_MOTIF_ALLOWLIST = (
    "alpha_helix",
    "3_10_helix",
    "pi_helix",
    "beta_sheet_parallel",
    "beta_sheet_antiparallel",
    "beta_turn_I",
    "beta_turn_II",
    "beta_turn_I_prime",
    "beta_turn_II_prime",
    "gamma_turn_classic",
    "gamma_turn_inverse",
    "collagen_triple_helix",
)
APPROVED_MOTIFS = frozenset(PROTEIN_MOTIF_ALLOWLIST)
