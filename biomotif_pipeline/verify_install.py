#!/usr/bin/env python3
"""Check the existing protein-pipeline dependencies without modifying them."""
from __future__ import annotations

import argparse
import importlib
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dssp", default="auto", help="auto, an executable path/name, or none for Stage 3 only")
    args = parser.parse_args()
    ok = sys.version_info >= (3, 11)
    print(f"Python: {sys.version.split()[0]} (requires >=3.11)")
    for name in ("requests", "numpy", "Bio", "gemmi", "rdkit"):
        try:
            module = importlib.import_module(name)
            print(f"{name}: {getattr(module, '__version__', 'installed')}")
        except Exception as exc:
            print(f"{name}: MISSING ({exc})")
            ok = False
    if args.dssp.lower() != "none":
        candidates = ("mkdssp", "dssp") if args.dssp == "auto" else (args.dssp,)
        executable = next((shutil.which(item) or (item if Path(item).is_file() else None)
                           for item in candidates if shutil.which(item) or Path(item).is_file()), None)
        try:
            if executable is None:
                raise RuntimeError("executable not found")
            result = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=20, check=True)
            print(f"DSSP: {executable} {(result.stdout or result.stderr).strip()}")
        except Exception as exc:
            print(f"DSSP: MISSING/FAILED ({exc})")
            ok = False
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
