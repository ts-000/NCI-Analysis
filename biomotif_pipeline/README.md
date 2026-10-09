# Biomotif Pipeline

A workflow for collecting high-quality crystal structures, annotating selected backbone motifs, and preparing small local fragments for later Gaussian/NBO work.

The code keeps the three stages separate:

1. **Collect** an RCSB protein set deduplicated by entry-level cluster combinations.
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
- entry-level deduplication by the combination of 30% RCSB protein-entity clusters

Two selected entries can share an individual entity cluster if their full cluster
combinations differ. Stage 1 does not guarantee global entity-level nonredundancy.

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

Beta turns require nonhelical central residues and canonical torsions: three
angular deviations at most 30 degrees and the fourth at most 40 degrees, in
addition to the existing C-alpha and configured RMS limits. Classic gamma turns
have positive phi/negative psi; inverse gamma turns have negative phi/positive psi.
Sheet windows retain the existing collection scheme, but require an actual retained
cross-strand DSSP hydrogen bond (energy below -0.5 kcal/mol and N...O at most 3.5 A).
Metadata records its direction and energy without assuming a reciprocal bond.
Reannotate old Stage-2 registries to apply these classification corrections.

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

`cap_report.json` and the `boundary_affected`/`boundary_type` columns in
`atom_map.csv` identify the directly altered boundary heavy atoms. H caps can
produce boundary amines and aldehydes; they remain the existing approximation.
Environment reports include nearby neutral protein residues and ligands as well
as waters (including DOD), ions (including CU), and charged residues. Reporting
does not automatically add those neighbours to the fragment.
Minimum contact distances use deposited heavy atoms, excluding H/D. A listed
contact or "strong" geometry label does not establish an NBO energy.

Unverified polymer endpoints, disulfide-linked CYS, and neutral ARG above pH 12.5
are rejected for manual preparation instead of guessing unsupported chemistry.
PDB hydrogen names are unique within each residue; Gaussian atom order, labels,
freeze flags, and the capping/hydrogenation procedure are preserved.

Verify installation and run the protein-only regression tests with:

```bash
python verify_install.py --dssp mkdssp
python tests/selftest.py
```

The optional protein-only helper can rank prepared sites and write Gaussian input
files using the existing settings; it is not part of the three-stage launcher:

```bash
python select_gaussian_candidates_v2.py --stage3 stage3_fragments --max-files-per-motif 1
```

Its shortlist requires explicit successful validation and matching atom maps/freeze
masks, and uses deposited metadata for occupancy/B factors. Review motif identity,
boundary chemistry, and omitted environment before choosing a representative.
