# Pipeline 2 — Interaction Geometry and Mainframe Assembly

## Overview

Pipeline 2 takes the classified NBO interaction files from Pipeline 1, combines them with matching XYZ coordinates, and adds interaction-specific geometric descriptors.

Two scripts are provided:

- `mainframe_geometry_v9.py` — model-molecule series.
- `mainframe_geometry_motif.py` — biomolecular motifs.

## Input

Each classified CSV must have a matching XYZ file with the same system name.

Model example:

```text
classified_CF3_N1.csv
CF3_N1.xyz
```

Biomotif example:

```text
classified_motif_example.csv
motif_example.xyz
```

Empty filename examples are included under `example/` only to show the expected naming.

## Requirements

```bash
pip install numpy pandas
```

## Run

For the model-molecule mainframe:

```bash
python mainframe_geometry_v9.py
```

Output:

```text
mainframe.csv
```

For biomolecular motifs:

```bash
python mainframe_geometry_motif.py
```

Output:

```text
mainframe_biomotifs.csv
```

## Notes

- The model script uses the fixed substituent/conformer naming used in this project.
- Interaction geometries are calculated only for explicitly supported donor–acceptor subtypes; unsupported geometry fields remain blank.
- Keep model and biomotif input sets separate when running the two workflows.
