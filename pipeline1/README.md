# Pipeline 1 — Gaussian/NBO Interaction Extraction and Classification

## Overview

Pipeline 1 converts completed Gaussian log files containing an NBO second-order perturbation analysis into compact CSV tables of long-range donor–acceptor interactions.

The normal XYZ workflow performs four steps:

1. Extract one molecular geometry from each Gaussian `.log` file and write an `.xyz` file.
2. Extract the final NBO E(2) table associated with the selected Gaussian job and write a raw interaction `.csv` file.
3. Infer the covalent bond graph from the XYZ coordinates, then retain supported interactions that satisfy the energy and bond-separation filters.
4. Convert the retained NBO orbital labels into broad interaction classes such as `n→σ\\\*`, `n→π\\\*`, `σ→σ\\\*`, `σ→π\\\*`, `π→σ\\\*`, and `π→π\\\*`.

The scripts use only the Python standard library. No third-party Python package is required.

## Folder layout

```text
pipeline\\\_1/
├── chemtools.py
├── molecule\\\_visualizer.py
├── parser.py
├── filter.py
├── classifier.py
├── input/
│   ├── molecule\\\_1.log
│   ├── molecule\\\_2.log
│   └── ...
├── pdb/                     # needed only for pipeline-pdb
│   ├── molecule\\\_1.pdb
│   └── ...
└── output/                  # created automatically if absent
```

File stems must match. For example, `molecule\\\_1.log`, `molecule\\\_1.csv`, `molecule\\\_1.xyz`, and optionally `molecule\\\_1.pdb` must all use the same stem.

## Requirements

* Python 3.12 was used
* Gaussian log files must contain:

  * `Normal termination of Gaussian` for the relevant job.
* The NBO E(2) table must use the conventional Gaussian/NBO fixed-column line format expected by `parser.py`.

## Quick start

Place the `.log` files in `input/`, then start the command interface:

```bash
python chemtools.py
```

At the prompt, run:

```text
pipeline
```

Exit with:

```text
quit
```

The main results are written to `output/`.

## Known limitations and cautions

1. **The folder must contain `input/`.** The scripts do not create it automatically. The `output/` folder is created automatically.
2. **Linked Gaussian jobs require attention.** Both the geometry extractor and NBO parser stop at the first `Normal termination of Gaussian`. This is correct when the intended geometry and NBO table are within that first job, but later Link1 jobs are ignored.
3. **PDB `CONECT` records are not used.** PDB bonds are inferred from Cartesian distances, just as for XYZ. (This part requires improvement)
4. **Bond inference is approximate.** Unusual elements, metal coordination, compressed nonbonded contacts, or atypical bonding may need manual verification or a different topology source.
5. **Only the elements listed in `filter.py` have covalent radii.** An unlisted element causes filtering to fail explicitly.
6. **Reruns replace same-name filtered output.** Keep copies elsewhere when historical outputs must be preserved.

