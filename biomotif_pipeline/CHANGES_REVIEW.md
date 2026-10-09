# Six-point protein-only review patch — 2026-10-09

Base: `385d69492bc5267b953d7b736b8a3f017a202b37`.
The patch is limited to the six agreed review points. It preserves the three-stage
workflow and its H-cap model. No DNA implementation or Gaussian calculations were added.

## Exact changes

| Agreed point | Files | Change |
|---|---|---|
| 1. Missing tooling | `verify_install.py` (new), `tests/selftest.py` (new), `select_gaussian_candidates_v2.py` (new), `run_pipeline.ps1` | Restore installation checking and protein regression tests. Restore the existing local Gaussian helper as an optional protein-only script; keep its route, NBO settings, resource defaults, ranking/fallback policy and coordinate rendering. Require explicit validation, required reports, matching atom maps/freeze flags/counts/charge, and the producer's exact `heavy_atom_coordinate_mismatches` field. Read deposited quality metadata instead of prepared-PDB occupancy/B placeholders. Reject fractional integer fields. Add five native-exit checks to the PowerShell launcher. |
| 2. Motif and sheet correctness | `biomotif/stage2.py` | Swap classic/inverse gamma names while preserving the existing target magnitudes and tolerance. Require nonhelical beta-turn centres and three canonical angular errors <=30 degrees, fourth <=40, alongside the existing distance/RMS checks. Record angular errors and central DSSP codes. Preserve all DSSP offsets for window collection, add direction/energy records, and retain only sheet windows with an actual cross-window N...O contact <=3.5 A and DSSP energy <-0.5 kcal/mol. Reject overlapping units. Deduplicate directed bonds; remove the assumed reciprocal interaction. Retain the seed as `dssp_candidate_pair`; `dssp_hbond_pair` and `dssp_hydrogen_bonds` describe actual retained bonds. |
| 3. Hydrogen names and boundary provenance | `biomotif/stage3.py` | Number added PDB hydrogens per residue rather than per heavy parent. Preserve parent-based atom labels, coordinates, order, and freeze flags. Add `affected_heavy_atom_labels` to cap links: N at N-side cuts, C/O at C-side cuts. Append `boundary_affected` and `boundary_type` to `atom_map.csv`, preserving every existing column position and `source_kind` value. |
| 4. Unsupported chemistry guards | `biomotif/stage3.py` | Verify an unconnected endpoint against deposited label-sequence position and entity sequence length; otherwise raise `ambiguous_polymer_endpoint`. Detect selected disulfide CYS from deposited links or S...S <=2.3 A and raise `unsupported_disulfide`. Raise `unsupported_neutral_arginine` above pH12.5; keep the existing protonated-ARG branch. Unsupported cases require manual preparation. |
| 5. Environment reporting | `biomotif/stage3.py` | Recognize DOD water, CU ions, water entities and nonpolymer monatomic ions. Measure contacts using deposited heavy atoms, excluding H/D. Add `nearby_protein_residues` and `nearby_ligands` within the existing 5 A reporting radius. Preserve existing water cutoff and environment-inclusion switches; reporting does not automatically include these neighbours. |
| 6. Redundancy description | `README.md`, `stage1_collect.py` help text | State that deduplication operates on entry-level combinations of 30% entity clusters, not global uniqueness of each entity. The Stage-1 algorithm is unchanged. |

`README.md` also documents points 1–5, the optional helper and reannotation of
older Stage-2 registries. This file is the requested change log.

## File inventory

Existing files changed (added/removed lines):

- `README.md`: 47/3.
- `biomotif/stage2.py`: 69/14.
- `biomotif/stage3.py`: 68/7.
- `run_pipeline.ps1`: 5/0.
- `stage1_collect.py`: 1/1 (help wording only).

Added: `verify_install.py` (42 lines), `tests/selftest.py` (300 lines),
`select_gaussian_candidates_v2.py` (373 lines), and this change log.

Unchanged byte-for-byte: `biomotif/__init__.py`, `biomotif/common.py`,
`biomotif/stage1.py`, `config.toml`, `requirements.txt`, `environment.yml`,
`run_pipeline.sh`, `stage2_annotate.py`, `stage3_prepare.py`.
No files outside `biomotif_pipeline/` are changed in the GitHub patch.

## Verification

- 12 regression tests passed, including fresh DSSP annotation of 10 cached parent
  structures and preparation replay of 15 existing fragments: the original 11,
  three recommended replacements and the optional pure-pi window.
- Eight valid original examples remain available, with corrected gamma names.
  The duplicate/helical beta-I, invalid antiparallel window and helical beta-II-prime
  are excluded from their previously assigned types. The three replacements
  (3ZOJ267–270, 1R6J264–265/200–201, 3QPA74–77) and optional 3ZOJ186–189 are retained.
- Prepared counts, charge/multiplicity, source heavy coordinates, atom labels/order,
  unique PDB names, freeze masks and boundary tags were checked against existing records.
- Independent original-GitHub versus patched preparation comparison: all 15 core
  atom blocks match exactly, including H coordinates; maximum change is 0 A.
- Selector tests reject missing validation, recorded moved atoms, inconsistent counts
  or charge, and wrong freeze flags. Ranking reads actual source B metadata.
- Installation verifier, command-line help, syntax checks and `git diff --check` pass.
  PowerShell failure simulation stops after the first failed command with nonzero exit.
- Main replay environment: Python3.11.15, NumPy2.4.6, Biopython1.87, Gemmi0.7.5,
  RDKit2026.03.4, DSSP4.6.1 with its activated Conda data directory.

The tests read original structures and write separate temporary/test outputs.
Original computer pipeline outputs and Drive files were not edited. Existing
Gaussian calculations were not rerun, and no claim of identical future live-RCSB
collection or Gaussian energies is made. Reannotate old registries to apply the
corrected identities; historical output folders are not silently rewritten.

Classification reference: the beta angular criterion follows the
[PROMOTIF author manual](https://www.uoxray.uoregon.edu/local/manuals/promotif/document_2.html);
gamma names follow [Guruprasad and Rajkumar, 2000](https://doi.org/10.1007/BF03404909).
