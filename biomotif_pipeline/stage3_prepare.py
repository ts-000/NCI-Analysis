#!/usr/bin/env python3
from __future__ import annotations

import argparse
import tomllib
from pathlib import Path
from typing import Sequence

from biomotif.common import setup_logging
from biomotif.stage3 import PreparationConfig, QualityConfig, SelectionConfig, run_stage3


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stage 3: extract approved local motifs, hydrogen-cap artificial cuts, add hydrogens locally, and export Gaussian-ready coordinates."
    )
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("config.toml"))
    parser.add_argument("--stage2", type=Path, default=Path("stage2_annotations"))
    parser.add_argument("--output", type=Path, default=Path("stage3_fragments"))
    parser.add_argument("--max-per-motif", type=int, default=None, help="0 prepares all passing sites")
    parser.add_argument("--include-direct-waters", action="store_true")
    parser.add_argument("--allow-oversize", action="store_true", help="Manual override for fragments above the hard atom limit")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    args = parser.parse_args(argv)
    setup_logging(args.log_level)
    if args.max_per_motif is not None and args.max_per_motif < 0:
        parser.error("--max-per-motif cannot be negative")
    with args.config.open("rb") as handle:
        cfg = tomllib.load(handle)
    selection_cfg = cfg["stage3"]["selection"]
    quality_cfg = cfg["stage3"]["quality"]
    preparation_cfg = cfg["stage3"]["preparation"]
    selection = SelectionConfig(
        max_per_motif=(
            args.max_per_motif
            if args.max_per_motif is not None
            else int(selection_cfg["max_per_motif"])
        )
    )
    quality = QualityConfig(**quality_cfg)
    preparation = PreparationConfig(
        **preparation_cfg,
        include_direct_waters=args.include_direct_waters,
        allow_oversize=args.allow_oversize,
    )
    sites = args.stage2 / "sites.jsonl"
    if not sites.exists():
        parser.error(f"Stage-2 site registry not found: {sites}")
    summary = run_stage3(sites, args.output, selection, quality, preparation)
    return 0 if summary["sites_prepared"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
