#!/usr/bin/env bash
set -euo pipefail

DSSP_BIN="${DSSP_BIN:-mkdssp}"
python verify_install.py --dssp "$DSSP_BIN"
python tests/selftest.py
python stage1_collect.py --output stage1_gold
python stage2_annotate.py --stage1 stage1_gold --output stage2_annotations --dssp "$DSSP_BIN" --workers 2
python stage3_prepare.py --stage2 stage2_annotations --output stage3_fragments --max-per-motif 0
