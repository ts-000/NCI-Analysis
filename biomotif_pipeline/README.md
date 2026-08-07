# Biomotif Pipeline

A workflow for collecting high-quality crystal structures, annotating selected backbone motifs, and preparing small local fragments for later Gaussian/NBO work.

The code keeps the three stages separate:

1. **Collect** a nonredundant RCSB protein set.
2. **Annotate** approved motifs with DSSP and local geometry.
3. **Prepare** capped, hydrogenated fragments while preserving deposited heavy-atom coordinates.

The pipeline stops at Gaussian-ready coordinates.

## Motifs written to the final registry

- `alpha_helix`
- `3_10_helix`
- `pi_helix`
- `beta_sheet_parallel`
- `beta_sheet_antiparallel`
- `beta_turn_I`
- `beta_turn_II`
- `beta_turn_I_prime`
- `beta_turn_II_prime`
- `gamma_turn_classic`
- `gamma_turn_inverse`
- `collagen_triple_helix`

Internal DSSP and geometry labels do not become output folders unless they map to one of these names.


## Run the pipeline

### Stage 1 — collect structures

```bash
python stage1_collect.py \
  --output stage1_gold \
  --protein-target 100 \
  --candidate-pool 2000 \
  --workers 8
```

Default filters are set in `config.toml`:

- X-ray diffraction
- resolution at or below 2.0 Å
- R-free at or below 0.25
- missing R-free rejected
- minimum protein length of 30 residues
- 30% RCSB sequence-cluster de-redundancy

Stage 1 writes the selected manifest under:

```text
stage1_gold/manifests/selected_structures.csv
```

### Stage 2 — annotate motifs

```bash
python stage2_annotate.py \
  --stage1 stage1_gold \
  --output stage2_annotations \
  --dssp mkdssp \
  --workers 2
```

Important outputs:

```text
stage2_annotations/sites.jsonl
stage2_annotations/tables/sites.csv
stage2_annotations/tables/protein_residues.csv
stage2_annotations/tables/protein_geometry.csv
stage2_annotations/summary.json
stage2_annotations/errors.json
```

A Stage-2 site is a structural/geometric assignment.

### Stage 3 — prepare local fragments

```bash
python stage3_prepare.py \
  --stage2 stage2_annotations \
  --output stage3_fragments \
  --max-per-motif 0
```

`--max-per-motif 0` prepares every passing site. Positive value for a smaller geometry-diverse subset.

For each successful site, Stage 3 writes:

```text
source_fragment.cif
fragment_capped_heavy.pdb
fragment_prepared.pdb
fragment_prepared.xyz
gaussian_coordinates.txt
atom_map.csv
freeze_mask.csv
site_metadata.json
protonation_report.json
charge_multiplicity.json
cap_report.json
environment_report.json
validation_report.json
```

The heavy atoms retain their deposited coordinates. Artificial peptide cuts receive hydrogen link caps. Added hydrogens are movable in `gaussian_coordinates.txt`; heavy atoms are marked frozen.