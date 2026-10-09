#!/usr/bin/env python3
from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

from biomotif.common import setup_logging, write_json
from biomotif.stage1 import BranchConfig, run_stage1


def main() -> int:
    parser = argparse.ArgumentParser(description="Stage 1: collect proteins deduplicated by entry-level RCSB cluster combinations.")
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("config.toml"))
    parser.add_argument("--output", type=Path, default=Path("stage1_gold"))
    parser.add_argument("--protein-target", type=int)
    parser.add_argument("--candidate-pool", type=int)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--download-structure-factors", action="store_true")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    args = parser.parse_args()
    setup_logging(args.log_level)

    with args.config.open("rb") as handle:
        cfg = tomllib.load(handle)
    project = cfg["project"]
    raw = cfg["stage1"]["protein"]
    protein = BranchConfig(
        name="protein",
        resolution_max=float(raw["resolution_max"]),
        rfree_max=float(raw["rfree_max"]),
        min_length=int(raw["min_length"]),
        target=args.protein_target or int(project["protein_target"]),
    )
    result = run_stage1(
        output=args.output,
        protein=protein,
        candidate_pool=args.candidate_pool or int(project["candidate_pool"]),
        workers=args.workers or int(project["workers"]),
        cluster_identity=int(raw["sequence_cluster_identity"]),
        download_structure_factors=args.download_structure_factors,
    )
    write_json(args.output / "effective_config.json", {"protein": protein, "project": project})
    return 0 if result["selected"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
