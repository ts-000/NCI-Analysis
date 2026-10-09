#!/usr/bin/env python3
"""Rank Stage-3 fragments and write one combined OPT/FREQ/NBO7 input each."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

ROUTE = "#p wb97xd/6-31g(d,p) opt nosymm freq pop=nbo7read"
NBO_BLOCK = "$NBO e2pert=0.1 dipole spnbo nce ARCHIVE $END"
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
AVOID_AA = {"HIS", "ASP", "GLU", "LYS", "ARG"}
IDEAL_AA = {"ALA", "GLY", "VAL", "LEU", "ILE", "PHE"}


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def find(data: Any, *keys: str) -> Any:
    wanted = {key.lower() for key in keys}
    if isinstance(data, dict):
        for key, value in data.items():
            if key.lower() in wanted and value not in (None, ""):
                return value
        for value in data.values():
            result = find(value, *keys)
            if result not in (None, ""):
                return result
    elif isinstance(data, list):
        for value in data:
            result = find(value, *keys)
            if result not in (None, ""):
                return result
    return None


def values_matching(data: Any, *words: str) -> list[Any]:
    found: list[Any] = []
    if isinstance(data, dict):
        for key, value in data.items():
            if any(word.lower() in key.lower() for word in words):
                found.append(value)
            found.extend(values_matching(value, *words))
    elif isinstance(data, list):
        for value in data:
            found.extend(values_matching(value, *words))
    return found


def present(value: Any) -> bool:
    if value in (None, False, 0, ""):
        return False
    if isinstance(value, str):
        return value.strip().lower() not in {"none", "false", "no", "ok", "pass", "passed", "[]", "{}"}
    return bool(value)


def number(value: Any, integer: bool = False) -> float | int | None:
    try:
        result = float(value)
        if not math.isfinite(result):
            return None
        if integer and not result.is_integer():
            return None
        return int(result) if integer else result
    except (TypeError, ValueError):
        return None


def safe_name(value: Any, limit: int = 50) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "unknown")).strip("_")
    return (text or "unknown")[:limit]


def read_coordinates(site: Path) -> list[dict[str, Any]]:
    path = site / "gaussian_coordinates.txt"
    atoms: list[dict[str, Any]] = []
    if not path.exists():
        return atoms
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        body, _, label = line.partition("!")
        fields = body.split()
        if len(fields) < 5:
            continue
        try:
            atoms.append({
                "element": fields[0].capitalize(), "flag": int(fields[1]),
                "x": float(fields[2]), "y": float(fields[3]), "z": float(fields[4]),
                "label": label.strip(),
            })
        except ValueError:
            pass
    return atoms


def pdb_info(path: Path) -> tuple[set[str], bool, list[float], list[float]]:
    residues: set[str] = set()
    altloc = False
    occupancies: list[float] = []
    bfactors: list[float] = []
    if not path.exists():
        return residues, altloc, occupancies, bfactors
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom = line[12:16].strip().upper()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element == "H" or atom.startswith("H"):
            continue
        residues.add(line[17:20].strip().upper())
        altloc |= bool(line[16:17].strip())
        try:
            occupancies.append(float(line[54:60]))
            bfactors.append(float(line[60:66]))
        except ValueError:
            pass
    return residues, altloc, occupancies, bfactors


def selector_residues(metadata: dict[str, Any]) -> set[str]:
    selectors = find(metadata, "residue_selectors")
    if not isinstance(selectors, list):
        return set()
    return {
        str(item.get("resname") or item.get("residue_name") or "").upper()
        for item in selectors if isinstance(item, dict)
    } - {""}


def report_count(report: dict[str, Any], *words: str) -> int:
    total = 0
    for value in values_matching(report, *words):
        if isinstance(value, (list, dict)):
            total += len(value)
        elif present(value):
            total += 1
    return total


def coordinates_match_map(site: Path, atoms: list[dict[str, Any]]) -> bool:
    try:
        with (site / "atom_map.csv").open(newline="", encoding="utf-8") as handle:
            mapping = list(csv.DictReader(handle))
        with (site / "freeze_mask.csv").open(newline="", encoding="utf-8") as handle:
            frozen = list(csv.DictReader(handle))
        if not atoms or len(mapping) != len(atoms) or len(frozen) != len(atoms):
            return False
        for index, (atom, mapped, mask) in enumerate(zip(atoms, mapping, frozen), 1):
            expected = 0 if atom["element"] == "H" else -1
            if int(mapped["fragment_index"]) != index or int(mask["fragment_index"]) != index:
                return False
            if mapped["element"].capitalize() != atom["element"] or mask["element"].capitalize() != atom["element"]:
                return False
            if mapped["label"] != atom["label"] or mask["label"] != atom["label"]:
                return False
            if atom["flag"] != expected or int(mask["freeze_flag"]) != expected:
                return False
            for axis in ("x", "y", "z"):
                value = float(mapped[axis])
                if not math.isfinite(value) or not math.isfinite(atom[axis]) or abs(value - atom[axis]) > 1e-8:
                    return False
        return len({atom["label"] for atom in atoms}) == len(atoms)
    except (OSError, ValueError, KeyError, TypeError):
        return False


def evaluate(site: Path, min_occupancy: float, max_bfactor: float) -> dict[str, Any]:
    metadata = load_json(site / "site_metadata.json")
    charge_json = load_json(site / "charge_multiplicity.json")
    validation = load_json(site / "validation_report.json")
    protonation = load_json(site / "protonation_report.json")
    environment = load_json(site / "environment_report.json")
    caps = load_json(site / "cap_report.json")

    residues, pdb_altloc, _, _ = pdb_info(site / "fragment_prepared.pdb")
    residues |= selector_residues(metadata)
    atoms = read_coordinates(site)

    branch = str(metadata.get("branch") or "protein").lower()
    motif = str(find(metadata, "motif_type", "motif") or site.parent.name)
    pdb_id = str(find(metadata, "pdb_id", "entry_id", "structure_id") or "unknown").upper()
    site_id = str(find(metadata, "site_id") or site.name)
    charge = number(find(charge_json, "formal_charge", "total_formal_charge", "charge"), integer=True)
    multiplicity = number(find(charge_json, "multiplicity", "spin_multiplicity"), integer=True)
    atom_count = int(number(find(validation, "atom_count", "total_atom_count"), integer=True) or len(atoms))

    source_quality = metadata.get("quality") or {}
    if not isinstance(source_quality, dict):
        source_quality = {}
    min_occ = number(source_quality.get("min_occupancy"))
    if min_occ is None:
        preparation_quality = protonation.get("stage2_quality") or {}
        if isinstance(preparation_quality, dict):
            min_occ = number(preparation_quality.get("minimum_occupancy"))
    mean_b = number(source_quality.get("mean_local_bfactor"))

    def flagged(*words: str) -> bool:
        return any(present(value) for value in values_matching(validation, *words))

    altloc = pdb_altloc or any(present(value) for value in values_matching(metadata, "altloc", "alternate_location"))
    waters = report_count(environment, "direct_water", "essential_water", "nearby_water")
    ions = report_count(environment, "nearby_ion", "direct_ion", "ion_count")
    real_terminus = any(present(value) for value in values_matching(caps, "biological_termin", "genuine_termin", "real_termin"))
    fallback_ph = any(present(value) for value in values_matching(protonation, "fallback", "manual_review_required"))

    modified = bool(residues - STANDARD_AA)
    ideal_fraction = len(residues & IDEAL_AA) / len(residues) if residues else 0.0

    validation_state = validation.get("validation_pass")
    hard = {
        "protein_branch": branch == "protein",
        "required_reports_present": all((metadata, charge_json, validation, protonation, environment, caps)),
        "validation_passed": validation_state is True,
        "coordinates_present": bool(atoms),
        "coordinates_match_map": coordinates_match_map(site, atoms),
        "atom_count_matches": number(validation.get("atom_count"), integer=True) == len(atoms),
        "charge_multiplicity_present": charge is not None and multiplicity is not None,
        "charge_multiplicity_consistent": charge == number(validation.get("formal_charge"), integer=True)
            and multiplicity == number(validation.get("multiplicity"), integer=True),
        "no_missing_heavy_atoms": not flagged("missing_required", "missing_heavy", "missing_atom"),
        "no_altloc": not altloc,
        "supported_residues": not modified and not flagged("unsupported_residue", "unsupported_modified"),
        "valid_valence": not flagged("invalid_valence", "valence_error"),
        "no_hydrogen_clash": not flagged("hydrogen_clash", "severe_clash"),
        "heavy_atoms_unchanged": validation.get("heavy_atom_coordinate_mismatches") == []
            and not flagged("heavy_coordinate_mismatch", "heavy_atom_moved"),
        "within_150_atoms": 0 < atom_count <= 150,
    }
    soft = {
        "no_histidine": "HIS" not in residues,
        "no_Asp_Glu_Lys_Arg": not bool(residues & (AVOID_AA - {"HIS"})),
        "neutral_charge": charge == 0,
        "no_nearby_ion": ions == 0,
        "no_direct_water": waters == 0,
        "not_real_terminus": not real_terminus,
        "occupancy_ok": min_occ is not None and min_occ >= min_occupancy,
        "bfactor_ok": mean_b is not None and mean_b <= max_bfactor,
        "experimental_ph": not fallback_ph,
    }
    hard_pass = all(hard.values())
    strict_pass = hard_pass and all(soft.values())
    score = 10 * sum(hard.values()) + 2 * sum(soft.values()) + 3 * ideal_fraction

    return {
        "site": site, "branch": branch, "pdb_id": pdb_id, "motif": motif, "site_id": site_id,
        "residues": ";".join(sorted(residues)), "charge": charge, "multiplicity": multiplicity,
        "atom_count": atom_count, "min_occupancy": min_occ, "mean_bfactor": mean_b,
        "waters": waters, "ions": ions, "fallback_ph": fallback_ph,
        "hard_pass": hard_pass, "strict_pass": strict_pass, "score": round(score, 3),
        "failed_filters": ";".join(key for key, ok in {**hard, **soft}.items() if not ok),
        "atoms": atoms,
    }


def coordinates(atoms: list[dict[str, Any]]) -> str:
    lines = []
    for atom in atoms:
        label = f" ! {atom['label']}" if atom["label"] else ""
        lines.append(
            f" {atom['element']:<2s} {atom['flag']:2d} "
            f"{atom['x']:16.10f} {atom['y']:16.10f} {atom['z']:16.10f}{label}"
        )
    return "\n".join(lines)


def write_input(row: dict[str, Any], folder: Path, rank: int, extension: str) -> str:
    fallback = "_fallback" if not row["strict_pass"] else ""
    base = "_".join((
        safe_name(row["pdb_id"], 12), safe_name(row["motif"], 35),
        safe_name(row["site_id"], 24), f"r{rank:02d}{fallback}",
    ))
    filename = f"{base}.{extension}"
    text = (
        f"%chk={base}.chk\n%nprocshared=18\n%mem=12gb\n{ROUTE}\n\n"
        f"{row['pdb_id']} {row['motif']} {row['site_id']} OPT FREQ NBO7\n\n"
        f"{row['charge']} {row['multiplicity']}\n{coordinates(row['atoms'])}\n\n"
        f"{NBO_BLOCK}\n\n"
    )
    folder.mkdir(parents=True, exist_ok=True)
    (folder / filename).write_text(text, encoding="utf-8")
    return filename


def stage3_roots(project: Path, supplied: list[str], include_pilot: bool) -> list[Path]:
    roots = [Path(value) if Path(value).is_absolute() else project / value for value in supplied]
    if not roots:
        roots = [path for path in (project / "stage3_fragments", project / "stage3_protein") if path.exists()]
    if not roots:
        roots = [path for path in project.glob("stage3*") if path.is_dir() and "dna" not in path.name.lower()]
    return [path for path in roots if include_pilot or "pilot" not in path.name.lower()]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "motif_rank", "generated", "selection", "branch", "pdb_id", "motif", "site_id",
        "score", "strict_pass", "hard_pass", "charge", "multiplicity", "atom_count",
        "residues", "min_occupancy", "mean_bfactor", "waters", "ions", "fallback_ph",
        "failed_filters", "input_file", "site_directory",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                **{key: row.get(key, "") for key in fields},
                "site_directory": row["site"],
            })


def main() -> int:
    parser = argparse.ArgumentParser(description="Rank Stage-3 fragments and write combined OPT/FREQ/NBO7 inputs.")
    parser.add_argument("--stage3", nargs="*", default=[], help="Defaults to stage3_fragments or stage3_protein (protein only).")
    parser.add_argument("--output", default="gaussian_candidates_v2")
    parser.add_argument("--max-files-per-motif", type=int, default=0, help="0 = all strict-pass candidates.")
    parser.add_argument("--extension", choices=("gjf", "com"), default="gjf")
    parser.add_argument("--min-occupancy", type=float, default=0.90)
    parser.add_argument("--max-mean-bfactor", type=float, default=40.0)
    parser.add_argument("--include-pilot", action="store_true")
    args = parser.parse_args()

    project = Path(__file__).resolve().parent
    roots = stage3_roots(project, args.stage3, args.include_pilot)
    sites = sorted({path.parent for root in roots for path in root.rglob("site_metadata.json")})
    if not sites:
        parser.error("No Stage-3 site_metadata.json files found.")

    rows = [evaluate(site, args.min_occupancy, args.max_mean_bfactor) for site in sites]
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[row["motif"]].append(row)

    output = Path(args.output) if Path(args.output).is_absolute() else project / args.output
    output.mkdir(parents=True, exist_ok=True)

    for motif_rows in groups.values():
        motif_rows.sort(key=lambda row: (-row["strict_pass"], -row["hard_pass"], -row["score"], row["pdb_id"], row["site_id"]))
        for rank, row in enumerate(motif_rows, 1):
            row["motif_rank"] = rank
        selected = [row for row in motif_rows if row["strict_pass"]]
        selection = "strict_pass"
        if not selected:
            selected = [row for row in motif_rows if row["hard_pass"]][:1]
            selection = "best_fallback"
        elif args.max_files_per_motif:
            selected = selected[:args.max_files_per_motif]
        for row in selected:
            row["generated"] = True
            row["selection"] = selection
            row["input_file"] = write_input(row, output / "inputs", row["motif_rank"], args.extension)

    rows.sort(key=lambda row: (row["motif"], row["motif_rank"]))
    write_csv(output / "candidate_ranking.csv", rows)
    generated = sum(bool(row.get("generated")) for row in rows)
    print(f"Scanned {len(rows)} fragments; generated {generated} combined Gaussian/NBO inputs.")
    print(f"Ranking: {output / 'candidate_ranking.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
