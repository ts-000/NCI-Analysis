#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import tomllib
import traceback
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Mapping, Sequence

from biomotif.common import (
    PROTEIN_MOTIF_ALLOWLIST,
    find_executable,
    read_csv,
    read_json,
    setup_logging,
    write_csv,
    write_json,
    write_jsonl,
)
from biomotif.stage2 import Stage2Config, process_structure

TABLES = ("structures", "approved_motifs", "protein_residues", "protein_geometry")


def load_config(path: Path) -> Stage2Config:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)["stage2"]
    fields = Stage2Config.__dataclass_fields__
    return Stage2Config(**{name: raw[name] for name in fields if name in raw})


def protein_records(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [dict(row) for row in rows if str(row.get("branch") or "protein") == "protein"]


def build_summary(
    records: Sequence[Mapping[str, Any]],
    structures: Sequence[Mapping[str, Any]],
    sites: Sequence[Mapping[str, Any]],
    errors: Sequence[Mapping[str, Any]],
    dssp: str | None,
) -> dict[str, Any]:
    failed = [str(row.get("pdb_id")) for row in errors]
    counts = Counter(str(site["motif_type"]) for site in sites)
    return {
        "input_proteins": len(records),
        "successfully_parsed_proteins": len(structures),
        "failed_proteins": len(records) - len(structures),
        "failed_protein_ids": failed,
        "approved_motif_counts": {name: counts.get(name, 0) for name in PROTEIN_MOTIF_ALLOWLIST},
        "approved_site_count": len(sites),
        "dssp_executable": dssp,
    }


def update_stage1_summary(stage1: Path, summary: Mapping[str, Any]) -> None:
    path = stage1 / "summary.json"
    current: dict[str, Any] = {}
    if path.exists():
        loaded = read_json(path)
        if isinstance(loaded, dict):
            current = loaded
    current["protein_selected"] = int(summary["input_proteins"])
    current["protein_stage2_successful"] = int(summary["successfully_parsed_proteins"])
    current["protein_stage2_failed"] = int(summary["failed_proteins"])
    current["protein_stage2_failed_ids"] = list(summary["failed_protein_ids"])
    write_json(path, current)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stage 2: annotate approved protein motifs with DSSP and local geometry.")
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("config.toml"))
    parser.add_argument("--stage1", type=Path, default=Path("stage1_gold"))
    parser.add_argument("--output", type=Path, default=Path("stage2_annotations"))
    parser.add_argument("--dssp", default="auto")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--max-structures", type=int, default=0)
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    args = parser.parse_args(argv)
    setup_logging(args.log_level)

    if args.workers <= 0:
        parser.error("--workers must be positive")
    if args.max_structures < 0:
        parser.error("--max-structures cannot be negative")

    manifest = args.stage1 / "manifests" / "selected_structures.csv"
    if not manifest.exists():
        parser.error(f"Stage-1 manifest not found: {manifest}")
    before = manifest.read_bytes()
    rows = protein_records(read_csv(manifest))
    records = rows[:args.max_structures] if args.max_structures else rows
    if not records:
        parser.error("No protein structures were found in the Stage-1 manifest")

    dssp = find_executable(args.dssp, ("mkdssp", "dssp"))
    if dssp is None:
        parser.error("DSSP is required; use --dssp mkdssp or pass its full path")

    cfg = load_config(args.config)
    args.output.mkdir(parents=True, exist_ok=True)
    tables: dict[str, list[dict[str, Any]]] = {name: [] for name in TABLES}
    sites: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    traces: list[str] = []

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = {pool.submit(process_structure, row, args.stage1, cfg, dssp): row for row in records}
        for done, future in enumerate(as_completed(jobs), 1):
            row = jobs[future]
            try:
                result, new_sites = future.result()
                for name, values in result.items():
                    tables[name].extend(values)
                sites.extend(new_sites)
            except Exception as exc:
                errors.append({"pdb_id": row.get("pdb_id"), "branch": "protein", "error": str(exc)})
                logging.error("%s | protein | %s", row.get("pdb_id"), str(exc))
                traces.append(f"\n===== {row.get('pdb_id')} (protein) =====\n{traceback.format_exc()}")
            if done % 10 == 0 or done == len(jobs):
                logging.info("Stage 2: %d/%d structures", done, len(jobs))

    unexpected = sorted({site.get("motif_type") for site in sites} - set(PROTEIN_MOTIF_ALLOWLIST))
    if unexpected:
        raise RuntimeError(f"Unapproved motif categories reached the final output: {unexpected}")

    for name, rows_out in tables.items():
        write_csv(args.output / "tables" / f"{name}.csv", rows_out)
        write_jsonl(args.output / "jsonl" / f"{name}.jsonl", rows_out)
    write_csv(args.output / "tables" / "sites.csv", sites)
    write_jsonl(args.output / "sites.jsonl", sites)
    write_json(args.output / "errors.json", errors)
    (args.output / "debug_tracebacks.log").write_text("".join(traces), encoding="utf-8")

    summary = build_summary(records, tables["structures"], sites, errors, dssp)
    write_json(args.output / "summary.json", summary)
    if not args.max_structures:
        update_stage1_summary(args.stage1, summary)
    if manifest.read_bytes() != before:
        raise RuntimeError("Stage 2 modified the Stage-1 manifest")
    return 0 if tables["structures"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
