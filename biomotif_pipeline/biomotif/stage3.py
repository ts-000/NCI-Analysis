from __future__ import annotations

import logging
import math
import statistics
import traceback
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import gemmi
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

from .common import (
    APPROVED_MOTIFS,
    distance,
    materialized_text_file,
    read_jsonl,
    safe_float,
    unit,
    write_csv,
    write_json,
)

STANDARD_PROTEIN = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "HYP",
}
WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
ION_NAMES = {"NA", "K", "LI", "MG", "CA", "ZN", "CL", "CU"}
CHARGED_RESIDUES = {"ASP", "GLU", "LYS", "ARG", "HIS", "CYS", "TYR"}
ATOMIC_NUMBERS = {
    "H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16,
    "CL": 17, "K": 19, "CA": 20, "ZN": 30, "MG": 12, "NA": 11,
    "SE": 34,
}
LINK_BOND_LENGTHS = {"N": 1.01, "C": 1.09}

# Heavy-atom aliases seen in PDB/mmCIF files.
ATOM_ALIASES: dict[str, dict[str, tuple[str, ...]]] = {
    "HYP": {"OD1": ("OD1", "OD", "OG")},
}

# Bond order values: 1.0 single, 2.0 double, 1.5 aromatic.
BACKBONE_BONDS = (("N", "CA", 1.0), ("CA", "C", 1.0), ("C", "O", 2.0))
SIDECHAIN_BONDS: dict[str, tuple[tuple[str, str, float], ...]] = {
    "GLY": (),
    "ALA": (("CA", "CB", 1.0),),
    "VAL": (("CA", "CB", 1.0), ("CB", "CG1", 1.0), ("CB", "CG2", 1.0)),
    "LEU": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD1", 1.0), ("CG", "CD2", 1.0)),
    "ILE": (("CA", "CB", 1.0), ("CB", "CG1", 1.0), ("CB", "CG2", 1.0), ("CG1", "CD1", 1.0)),
    "SER": (("CA", "CB", 1.0), ("CB", "OG", 1.0)),
    "THR": (("CA", "CB", 1.0), ("CB", "OG1", 1.0), ("CB", "CG2", 1.0)),
    "CYS": (("CA", "CB", 1.0), ("CB", "SG", 1.0)),
    "MET": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "SD", 1.0), ("SD", "CE", 1.0)),
    "MSE": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "SE", 1.0), ("SE", "CE", 1.0)),
    "PRO": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD", 1.0), ("CD", "N", 1.0)),
    "HYP": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD", 1.0), ("CD", "N", 1.0), ("CG", "OD1", 1.0)),
    "ASN": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "OD1", 2.0), ("CG", "ND2", 1.0)),
    "ASP": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "OD1", 2.0), ("CG", "OD2", 1.0)),
    "GLN": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD", 1.0), ("CD", "OE1", 2.0), ("CD", "NE2", 1.0)),
    "GLU": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD", 1.0), ("CD", "OE1", 2.0), ("CD", "OE2", 1.0)),
    "LYS": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD", 1.0), ("CD", "CE", 1.0), ("CE", "NZ", 1.0)),
    "ARG": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD", 1.0), ("CD", "NE", 1.0), ("NE", "CZ", 1.0), ("CZ", "NH1", 1.0), ("CZ", "NH2", 1.0)),
    "PHE": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD1", 1.5), ("CD1", "CE1", 1.5), ("CE1", "CZ", 1.5), ("CZ", "CE2", 1.5), ("CE2", "CD2", 1.5), ("CD2", "CG", 1.5)),
    "TYR": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD1", 1.5), ("CD1", "CE1", 1.5), ("CE1", "CZ", 1.5), ("CZ", "CE2", 1.5), ("CE2", "CD2", 1.5), ("CD2", "CG", 1.5), ("CZ", "OH", 1.0)),
    "HIS": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "ND1", 1.5), ("ND1", "CE1", 1.5), ("CE1", "NE2", 1.5), ("NE2", "CD2", 1.5), ("CD2", "CG", 1.5)),
    "TRP": (("CA", "CB", 1.0), ("CB", "CG", 1.0), ("CG", "CD1", 1.5), ("CD1", "NE1", 1.5), ("NE1", "CE2", 1.5), ("CE2", "CD2", 1.5), ("CD2", "CG", 1.5), ("CD2", "CE3", 1.5), ("CE3", "CZ3", 1.5), ("CZ3", "CH2", 1.5), ("CH2", "CZ2", 1.5), ("CZ2", "CE2", 1.5)),
}


@dataclass
class SelectionConfig:
    max_per_motif: int = 0


@dataclass
class QualityConfig:
    min_occupancy: float = 0.90
    preferred_atoms: int = 100
    maximum_atoms: int = 150


@dataclass
class PreparationConfig:
    fallback_ph: float = 7.0
    include_direct_waters: bool = False
    water_oxygen_distance: float = 3.5
    hydrogen_caps: bool = True
    allow_oversize: bool = False


@dataclass
class SourceAtom:
    chain_id: str
    seq_id: int
    insertion_code: str
    resname: str
    atom_name: str
    element: str
    coord: np.ndarray
    occupancy: float
    bfactor: float
    altloc: str
    had_alternates: bool
    serial: int

    @property
    def residue_key(self) -> tuple[str, int, str]:
        return self.chain_id, self.seq_id, self.insertion_code


class PreparationError(RuntimeError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code


class SourceIndex:
    """A coordinate index. It never hydrogenates or repairs the parent structure."""

    def __init__(self, path: Path) -> None:
        self.path = path.resolve()
        with materialized_text_file(self.path, suffix=".cif") as materialized:
            self.structure = gemmi.read_structure(str(materialized))
        if len(self.structure) == 0:
            raise PreparationError("source_parse_failed", "source structure contains no coordinate model")
        self.model = self.structure[0]
        self.residues: dict[tuple[str, int, str], dict[str, Any]] = {}
        self.chain_polymer_keys: dict[str, list[tuple[str, int, str]]] = {}
        serial = 0
        for chain in self.model:
            polymer_keys: list[tuple[str, int, str]] = []
            for residue in chain:
                key = (chain.name, int(residue.seqid.num), str(residue.seqid.icode).strip())
                grouped: dict[str, list[gemmi.Atom]] = {}
                for atom in residue:
                    grouped.setdefault(atom.name.strip(), []).append(atom)
                chosen: dict[str, SourceAtom] = {}
                for atom_name, candidates in grouped.items():
                    selected = max(
                        candidates,
                        key=lambda atom: (
                            float(atom.occ),
                            str(atom.altloc).strip() in {"", "A", "1", "\x00"},
                        ),
                    )
                    serial += 1
                    element = selected.element.name.upper() if selected.element.name else atom_name[0].upper()
                    chosen[atom_name] = SourceAtom(
                        chain_id=chain.name,
                        seq_id=int(residue.seqid.num),
                        insertion_code=str(residue.seqid.icode).strip(),
                        resname=residue.name.strip().upper(),
                        atom_name=atom_name,
                        element=element,
                        coord=np.array([selected.pos.x, selected.pos.y, selected.pos.z], dtype=float),
                        occupancy=float(selected.occ),
                        bfactor=float(selected.b_iso),
                        altloc=str(selected.altloc).strip(),
                        had_alternates=len(candidates) > 1,
                        serial=serial,
                    )
                record = {
                    "resname": residue.name.strip().upper(),
                    "atoms": chosen,
                    "gemmi_residue": residue,
                }
                self.residues[key] = record
                if all(name in chosen for name in ("N", "CA", "C")):
                    polymer_keys.append(key)
            if polymer_keys:
                self.chain_polymer_keys[chain.name] = polymer_keys

    @staticmethod
    def selector_key(selector: Mapping[str, Any]) -> tuple[str, int, str]:
        return (
            str(selector.get("chain_id") or ""),
            int(selector.get("seq_id")),
            str(selector.get("insertion_code") or ""),
        )

    def selected_keys(self, selectors: Sequence[Mapping[str, Any]]) -> set[tuple[str, int, str]]:
        return {self.selector_key(selector) for selector in selectors}

    def selected_atoms(self, keys: Iterable[tuple[str, int, str]]) -> list[SourceAtom]:
        atoms: list[SourceAtom] = []
        for key in keys:
            record = self.residues.get(key)
            if record:
                atoms.extend(record["atoms"].values())
        return atoms

    def atom(self, key: tuple[str, int, str], canonical_name: str) -> SourceAtom | None:
        record = self.residues.get(key)
        if record is None:
            return None
        aliases = ATOM_ALIASES.get(record["resname"], {}).get(canonical_name, (canonical_name,))
        for name in aliases:
            if name in record["atoms"]:
                return record["atoms"][name]
        return None

    def peptide_connected(self, left: tuple[str, int, str], right: tuple[str, int, str]) -> bool:
        carbon = self.atom(left, "C")
        nitrogen = self.atom(right, "N")
        return carbon is not None and nitrogen is not None and distance(carbon.coord, nitrogen.coord) <= 1.90

    def protein_blocks(self, selected: set[tuple[str, int, str]]) -> list[list[tuple[str, int, str]]]:
        blocks: list[list[tuple[str, int, str]]] = []
        for chain_id, ordered in self.chain_polymer_keys.items():
            selected_positions = [index for index, key in enumerate(ordered) if key in selected]
            if not selected_positions:
                continue
            current: list[tuple[str, int, str]] = []
            previous_position: int | None = None
            for position in selected_positions:
                key = ordered[position]
                if (
                    current
                    and previous_position is not None
                    and position == previous_position + 1
                    and self.peptide_connected(ordered[previous_position], key)
                ):
                    current.append(key)
                else:
                    if current:
                        blocks.append(current)
                    current = [key]
                previous_position = position
            if current:
                blocks.append(current)
        return blocks

    def biological_neighbors(self, block: Sequence[tuple[str, int, str]]) -> tuple[tuple[str, int, str] | None, tuple[str, int, str] | None]:
        if not block:
            return None, None
        ordered = self.chain_polymer_keys.get(block[0][0], [])
        positions = {key: index for index, key in enumerate(ordered)}
        start = positions[block[0]]
        end = positions[block[-1]]
        previous = ordered[start - 1] if start > 0 and self.peptide_connected(ordered[start - 1], block[0]) else None
        following = ordered[end + 1] if end + 1 < len(ordered) and self.peptide_connected(block[-1], ordered[end + 1]) else None
        for key, neighbor, side in ((block[0], previous, "N"), (block[-1], following, "C")):
            if neighbor is not None:
                continue
            residue = self.residues[key]["gemmi_residue"]
            entity = next((item for item in self.structure.entities if residue.subchain in item.subchains), None)
            length = len(entity.full_sequence) if entity is not None else 0
            position = residue.label_seq
            if not length or position is None or position != (1 if side == "N" else length):
                raise PreparationError(
                    "ambiguous_polymer_endpoint",
                    f"{side}-side endpoint at {key} is not a verified deposited sequence terminus; manual review required",
                )
        return previous, following

    def has_disulfide(self, key: tuple[str, int, str]) -> bool:
        sulfur = self.atom(key, "SG")
        if sulfur is None:
            return False
        for connection in self.structure.connections:
            if connection.type != gemmi.ConnectionType.Disulf:
                continue
            for partner in (connection.partner1, connection.partner2):
                address = (partner.chain_name, int(partner.res_id.seqid.num), str(partner.res_id.seqid.icode).strip())
                if address == key and partner.atom_name.strip() == sulfur.atom_name:
                    return True
        return any(
            other_key != key and record["resname"] == "CYS"
            and (other := self.atom(other_key, "SG")) is not None
            and distance(sulfur.coord, other.coord) <= 2.3
            for other_key, record in self.residues.items()
        )

    def write_subset_cif(self, keys: set[tuple[str, int, str]], destination: Path) -> None:
        out = gemmi.Structure()
        out.name = self.structure.name or destination.stem
        model = gemmi.Model("1")
        for chain in self.model:
            new_chain = gemmi.Chain(chain.name)
            for residue in chain:
                key = (chain.name, int(residue.seqid.num), str(residue.seqid.icode).strip())
                if key not in keys:
                    continue
                record = self.residues[key]
                new_residue = gemmi.Residue()
                new_residue.name = record["resname"]
                new_residue.seqid = gemmi.SeqId(key[1], key[2] or " ")
                for source_atom in record["atoms"].values():
                    if source_atom.element == "H":
                        continue
                    atom = gemmi.Atom()
                    atom.name = source_atom.atom_name
                    atom.element = gemmi.Element(source_atom.element)
                    atom.pos = gemmi.Position(*map(float, source_atom.coord))
                    atom.occ = source_atom.occupancy
                    atom.b_iso = source_atom.bfactor
                    new_residue.add_atom(atom)
                new_chain.add_residue(new_residue)
            if len(new_chain):
                model.add_chain(new_chain)
        out.add_model(model)
        out.setup_entities()
        destination.parent.mkdir(parents=True, exist_ok=True)
        out.make_mmcif_document().write_file(str(destination))

    def nearby_environment(
        self,
        selected: set[tuple[str, int, str]],
        centers: Sequence[np.ndarray],
        water_cutoff: float,
    ) -> dict[str, Any]:
        waters: list[dict[str, Any]] = []
        ions: list[dict[str, Any]] = []
        charged: list[dict[str, Any]] = []
        protein: list[dict[str, Any]] = []
        ligands: list[dict[str, Any]] = []
        if not centers:
            centers = [atom.coord for atom in self.selected_atoms(selected) if atom.element not in {"H", "D"}]
        for key, record in self.residues.items():
            if key in selected:
                continue
            atoms = [atom for atom in record["atoms"].values() if atom.element not in {"H", "D"}]
            if not atoms or not centers:
                continue
            minimum = min(distance(atom.coord, center) for atom in atoms for center in centers)
            item = {
                "chain_id": key[0], "seq_id": key[1], "insertion_code": key[2],
                "resname": record["resname"], "minimum_distance_A": minimum,
            }
            is_water = record["resname"] in WATER_NAMES or record["gemmi_residue"].entity_type == gemmi.EntityType.Water
            is_ion = record["resname"] in ION_NAMES or (
                record["gemmi_residue"].entity_type != gemmi.EntityType.Polymer
                and record["resname"] not in STANDARD_PROTEIN
                and len(atoms) == 1 and atoms[0].element not in {"C", "N", "O", "S", "P"}
            )
            if is_water and minimum <= water_cutoff:
                waters.append(item)
            elif is_ion and minimum <= 5.0:
                ions.append(item)
            elif record["resname"] in CHARGED_RESIDUES and minimum <= 5.0:
                charged.append(item)
            if not is_water and not is_ion and minimum <= 5.0:
                if record["gemmi_residue"].entity_type == gemmi.EntityType.Polymer or record["resname"] in STANDARD_PROTEIN:
                    protein.append(item)
                else:
                    ligands.append(item)
        return {
            "nearby_waters": sorted(waters, key=lambda item: item["minimum_distance_A"]),
            "nearby_ions": sorted(ions, key=lambda item: item["minimum_distance_A"]),
            "nearby_charged_residues": sorted(charged, key=lambda item: item["minimum_distance_A"]),
            "nearby_protein_residues": sorted(protein, key=lambda item: item["minimum_distance_A"]),
            "nearby_ligands": sorted(ligands, key=lambda item: item["minimum_distance_A"]),
        }


@dataclass
class FragmentAtom:
    element: str
    coord: np.ndarray
    atom_name: str
    chain_id: str
    seq_id: int
    insertion_code: str
    resname: str
    source_serial: int | None
    source_kind: str
    occupancy: float | None = None
    bfactor: float | None = None
    label: str = ""


@dataclass
class BuildResult:
    molecule: Chem.Mol
    initial_atoms: list[FragmentAtom]
    capped_atoms: list[FragmentAtom]
    cap_report: dict[str, Any]
    protonation_report: dict[str, Any]
    environment_report: dict[str, Any]
    source_heavy_coordinates: dict[str, np.ndarray]


def _bond_type(order: float) -> Chem.BondType:
    if order == 1.0:
        return Chem.BondType.SINGLE
    if order == 2.0:
        return Chem.BondType.DOUBLE
    if order == 1.5:
        return Chem.BondType.AROMATIC
    raise ValueError(order)


def _required_atoms(resname: str) -> set[str]:
    if resname in STANDARD_PROTEIN:
        atoms = {"N", "CA", "C", "O"}
        for left, right, _ in SIDECHAIN_BONDS[resname]:
            atoms.update((left, right))
        return atoms
    return set()


def _choose_ph(site: Mapping[str, Any], cfg: PreparationConfig) -> tuple[float, str, bool]:
    exact = safe_float(site.get("crystallization_ph"))
    if exact is not None:
        return exact, "exact_crystallization_pH", False
    minimum = safe_float(site.get("crystallization_ph_min"))
    maximum = safe_float(site.get("crystallization_ph_max"))
    if minimum is not None and maximum is not None:
        return (minimum + maximum) / 2.0, "crystallization_pH_range_midpoint", False
    return cfg.fallback_ph, "fallback_pH_7.0", True


def _centers_for_site(source: SourceIndex, site: Mapping[str, Any]) -> list[np.ndarray]:
    centers: list[np.ndarray] = []
    for selector in site.get("central_interaction_atoms", []):
        key = source.selector_key(selector)
        atom = source.atom(key, str(selector.get("atom_name") or ""))
        if atom is not None:
            centers.append(atom.coord)
    return centers


def validate_local_site(
    source: SourceIndex,
    site: Mapping[str, Any],
    cfg: QualityConfig,
) -> dict[str, Any]:
    if str(site.get("branch") or "protein") != "protein":
        raise PreparationError("unsupported_branch", "Stage 3 accepts protein sites only")
    motif_type = str(site.get("motif_type"))
    if motif_type not in APPROVED_MOTIFS:
        raise PreparationError("unapproved_motif", f"unapproved motif type {motif_type}")
    selectors = site.get("residue_selectors") or []
    if not selectors:
        raise PreparationError("empty_fragment", "site has no residue selectors")
    keys = source.selected_keys(selectors)
    missing_residues = [selector for selector in selectors if source.selector_key(selector) not in source.residues]
    if missing_residues:
        raise PreparationError("missing_residue", f"source residues not found: {missing_residues}")

    missing_atoms: list[str] = []
    minimum_occupancy = 1.0
    modified: list[str] = []
    for key in keys:
        record = source.residues[key]
        resname = record["resname"]
        if resname not in STANDARD_PROTEIN:
            raise PreparationError(
                "unsupported_residue_inside_fragment",
                f"unsupported residue {resname} inside selected fragment at {key[0]}:{key[1]}{key[2]}",
            )
        if resname in {"MSE", "HYP"}:
            modified.append(f"{key[0]}:{resname}{key[1]}{key[2]}")
        if resname == "CYS" and source.has_disulfide(key):
            raise PreparationError("unsupported_disulfide", f"disulfide-linked CYS at {key}; manual preparation required")
        for canonical in _required_atoms(resname):
            atom = source.atom(key, canonical)
            if atom is None:
                missing_atoms.append(f"{key[0]}:{resname}{key[1]}{key[2]}:{canonical}")
            else:
                minimum_occupancy = min(minimum_occupancy, atom.occupancy)
    if missing_atoms:
        raise PreparationError("missing_required_heavy_atom", "missing required heavy atoms: " + ", ".join(missing_atoms[:12]))
    if minimum_occupancy < cfg.min_occupancy:
        raise PreparationError(
            "occupancy_below_threshold",
            f"minimum local occupancy {minimum_occupancy:.2f} is below {cfg.min_occupancy:.2f}",
        )

    central_altlocs: list[str] = []
    for selector in site.get("central_interaction_atoms", []):
        key = source.selector_key(selector)
        atom_name = str(selector.get("atom_name") or "")
        atom = source.atom(key, atom_name)
        if atom is None:
            raise PreparationError(
                "missing_central_atom",
                f"missing central atom {atom_name} at {key[0]}:{key[1]}{key[2]}",
            )
        if atom.had_alternates:
            central_altlocs.append(f"{key[0]}:{key[1]}{key[2]}:{atom_name}")
    if central_altlocs:
        raise PreparationError(
            "unresolved_central_altloc",
            "central motif atoms have alternate conformations: " + ", ".join(central_altlocs),
        )
    return {
        "selected_residue_count": len(keys),
        "minimum_occupancy": minimum_occupancy,
        "modified_residues": modified,
        "local_quality_pass": True,
    }


def _add_atom(
    rw: Chem.RWMol,
    conformer: Chem.Conformer,
    atom_record: FragmentAtom,
    aromatic: bool = False,
) -> int:
    atom = Chem.Atom(atom_record.element.title())
    atom.SetIsAromatic(aromatic)
    index = rw.AddAtom(atom)
    conformer.SetAtomPosition(index, Point3D(*map(float, atom_record.coord)))
    return index


def _set_aromatic_atoms(rw: Chem.RWMol) -> None:
    for bond in rw.GetBonds():
        if bond.GetBondType() == Chem.BondType.AROMATIC:
            bond.SetIsAromatic(True)
            bond.GetBeginAtom().SetIsAromatic(True)
            bond.GetEndAtom().SetIsAromatic(True)


def _add_bond_if_absent(rw: Chem.RWMol, left: int, right: int, order: float) -> None:
    if left == right or rw.GetBondBetweenAtoms(left, right) is not None:
        return
    rw.AddBond(left, right, _bond_type(order))


def _configure_protonation(
    rw: Chem.RWMol,
    atom_index: Mapping[tuple[tuple[str, int, str], str], int],
    source: SourceIndex,
    selected_keys: set[tuple[str, int, str]],
    ph: float,
    terminal_states: Mapping[tuple[str, int, str], dict[str, bool]],
) -> tuple[list[dict[str, Any]], bool]:
    states: list[dict[str, Any]] = []
    manual_review = False

    def idx(key: tuple[str, int, str], name: str) -> int | None:
        atom = source.atom(key, name)
        if atom is None:
            return None
        return atom_index.get((key, atom.atom_name))

    for key in sorted(selected_keys):
        resname = source.residues[key]["resname"]
        state: dict[str, Any] = {
            "chain_id": key[0], "seq_id": key[1], "insertion_code": key[2], "resname": resname,
        }
        if resname == "ASP":
            protonated = ph < 4.0
            oxygen = idx(key, "OD2")
            if oxygen is not None and not protonated:
                rw.GetAtomWithIdx(oxygen).SetFormalCharge(-1)
            state["state"] = "ASH" if protonated else "ASP-"
        elif resname == "GLU":
            protonated = ph < 4.4
            oxygen = idx(key, "OE2")
            if oxygen is not None and not protonated:
                rw.GetAtomWithIdx(oxygen).SetFormalCharge(-1)
            state["state"] = "GLH" if protonated else "GLU-"
        elif resname == "LYS":
            protonated = ph <= 10.5
            nitrogen = idx(key, "NZ")
            if nitrogen is not None and protonated:
                rw.GetAtomWithIdx(nitrogen).SetFormalCharge(1)
            state["state"] = "LYS+" if protonated else "LYN"
        elif resname == "ARG":
            protonated = ph <= 12.5
            if not protonated:
                raise PreparationError("unsupported_neutral_arginine", f"neutral ARG at pH {ph:g} requires manual preparation at {key}")
            carbon = idx(key, "CZ")
            if carbon is not None and protonated:
                rw.GetAtomWithIdx(carbon).SetFormalCharge(1)
            state["state"] = "ARG+" if protonated else "ARG0"
        elif resname == "CYS":
            deprotonated = ph > 8.3
            sulfur = idx(key, "SG")
            if sulfur is not None and deprotonated:
                rw.GetAtomWithIdx(sulfur).SetFormalCharge(-1)
            state["state"] = "CYM" if deprotonated else "CYS"
        elif resname == "TYR":
            deprotonated = ph > 10.1
            oxygen = idx(key, "OH")
            if oxygen is not None and deprotonated:
                rw.GetAtomWithIdx(oxygen).SetFormalCharge(-1)
            state["state"] = "TYM" if deprotonated else "TYR"
        elif resname == "HIS":
            nd1 = idx(key, "ND1")
            ne2 = idx(key, "NE2")
            if ph < 6.0:
                if nd1 is not None:
                    rw.GetAtomWithIdx(nd1).SetNumExplicitHs(1)
                    rw.GetAtomWithIdx(nd1).SetNoImplicit(True)
                if ne2 is not None:
                    rw.GetAtomWithIdx(ne2).SetNumExplicitHs(1)
                    rw.GetAtomWithIdx(ne2).SetFormalCharge(1)
                    rw.GetAtomWithIdx(ne2).SetNoImplicit(True)
                state["state"] = "HIP"
            else:
                if ne2 is not None:
                    rw.GetAtomWithIdx(ne2).SetNumExplicitHs(1)
                    rw.GetAtomWithIdx(ne2).SetNoImplicit(True)
                if nd1 is not None:
                    rw.GetAtomWithIdx(nd1).SetNoImplicit(True)
                state["state"] = "HIE_default"
                state["manual_review_required"] = True
                manual_review = True
        elif resname == "TRP":
            ne1 = idx(key, "NE1")
            if ne1 is not None:
                rw.GetAtomWithIdx(ne1).SetNumExplicitHs(1)
                rw.GetAtomWithIdx(ne1).SetNoImplicit(True)
            state["state"] = "TRP"
        else:
            state["state"] = "standard"

        termini = terminal_states.get(key, {})
        if termini.get("genuine_n_terminus"):
            nitrogen = idx(key, "N")
            if nitrogen is not None and ph < 8.0:
                rw.GetAtomWithIdx(nitrogen).SetFormalCharge(1)
            state["genuine_n_terminus"] = True
        if termini.get("genuine_c_terminus"):
            oxt = idx(key, "OXT")
            if oxt is None:
                raise PreparationError(
                    "unsupported_genuine_c_terminus",
                    f"genuine C terminus lacks OXT at {key[0]}:{key[1]}{key[2]}",
                )
            if ph > 3.1:
                rw.GetAtomWithIdx(oxt).SetFormalCharge(-1)
            state["genuine_c_terminus"] = True
        states.append(state)
    return states, manual_review



def build_local_fragment(
    source: SourceIndex,
    site: Mapping[str, Any],
    quality: QualityConfig,
    preparation: PreparationConfig,
) -> BuildResult:
    quality_report = validate_local_site(source, site, quality)
    selected_keys = source.selected_keys(site["residue_selectors"])
    ph, ph_source, ph_manual_review = _choose_ph(site, preparation)
    centers = _centers_for_site(source, site)
    environment = source.nearby_environment(
        selected_keys,
        centers,
        preparation.water_oxygen_distance,
    )
    included_water_keys: set[tuple[str, int, str]] = set()
    if preparation.include_direct_waters:
        for item in environment["nearby_waters"]:
            key = (item["chain_id"], int(item["seq_id"]), item.get("insertion_code", ""))
            oxygen = next((atom for atom in source.residues[key]["atoms"].values() if atom.element == "O"), None)
            if oxygen is not None and oxygen.occupancy >= quality.min_occupancy:
                included_water_keys.add(key)

    rw = Chem.RWMol()
    initial_records: list[FragmentAtom] = []
    conformer = Chem.Conformer()
    atom_index: dict[tuple[tuple[str, int, str], str], int] = {}
    source_heavy_coordinates: dict[str, np.ndarray] = {}

    # Create selected heavy atoms only. Deposited hydrogens are intentionally ignored.
    for key in sorted(selected_keys):
        record = source.residues[key]
        for source_atom in sorted(record["atoms"].values(), key=lambda atom: atom.serial):
            if source_atom.element == "H":
                continue
            fragment_atom = FragmentAtom(
                element=source_atom.element,
                coord=source_atom.coord.copy(),
                atom_name=source_atom.atom_name,
                chain_id=key[0],
                seq_id=key[1],
                insertion_code=key[2],
                resname=record["resname"],
                source_serial=source_atom.serial,
                source_kind="source_heavy",
                occupancy=source_atom.occupancy,
                bfactor=source_atom.bfactor,
                label=f"{key[0]}:{record['resname']}{key[1]}{key[2]}:{source_atom.atom_name}",
            )
            index = _add_atom(rw, conformer, fragment_atom)
            atom_index[(key, source_atom.atom_name)] = index
            initial_records.append(fragment_atom)
            source_heavy_coordinates[fragment_atom.label] = fragment_atom.coord.copy()

    # Optional direct waters are the only environment atoms that can be included.
    for key in sorted(included_water_keys):
        oxygen = next((atom for atom in source.residues[key]["atoms"].values() if atom.element == "O"), None)
        if oxygen is None:
            continue
        fragment_atom = FragmentAtom(
            element="O",
            coord=oxygen.coord.copy(),
            atom_name="O",
            chain_id=key[0],
            seq_id=key[1],
            insertion_code=key[2],
            resname=source.residues[key]["resname"],
            source_serial=oxygen.serial,
            source_kind="source_heavy_water",
            occupancy=oxygen.occupancy,
            bfactor=oxygen.bfactor,
            label=f"{key[0]}:{source.residues[key]['resname']}{key[1]}{key[2]}:O",
        )
        index = _add_atom(rw, conformer, fragment_atom)
        atom_index[(key, oxygen.atom_name)] = index
        initial_records.append(fragment_atom)
        source_heavy_coordinates[fragment_atom.label] = fragment_atom.coord.copy()

    rw.AddConformer(conformer)

    def mapped_index(key: tuple[str, int, str], canonical: str) -> int:
        source_atom = source.atom(key, canonical)
        if source_atom is None:
            raise PreparationError("missing_required_heavy_atom", f"missing {canonical} at {key}")
        index = atom_index.get((key, source_atom.atom_name))
        if index is None:
            raise PreparationError("atom_mapping_failed", f"atom mapping failed for {key}:{canonical}")
        return index

    # Intra-residue bonds.
    for key in sorted(selected_keys):
        resname = source.residues[key]["resname"]
        bonds = BACKBONE_BONDS + SIDECHAIN_BONDS[resname]
        for left_name, right_name, order in bonds:
            _add_bond_if_absent(rw, mapped_index(key, left_name), mapped_index(key, right_name), order)
        oxt = source.atom(key, "OXT")
        if oxt is not None and (key, oxt.atom_name) in atom_index:
            _add_bond_if_absent(rw, mapped_index(key, "C"), atom_index[(key, oxt.atom_name)], 1.0)

    # Peptide bonds only inside the extracted local fragment.
    terminal_states: dict[tuple[str, int, str], dict[str, bool]] = {}
    cap_links: list[dict[str, Any]] = []
    capped_records = list(initial_records)
    blocks = source.protein_blocks(selected_keys)
    for block in blocks:
        for left, right in zip(block, block[1:]):
            _add_bond_if_absent(rw, mapped_index(left, "C"), mapped_index(right, "N"), 1.0)
        previous, following = source.biological_neighbors(block)
        terminal_states.setdefault(block[0], {})["genuine_n_terminus"] = previous is None
        terminal_states.setdefault(block[-1], {})["genuine_c_terminus"] = following is None

        if previous is not None:
            if not preparation.hydrogen_caps:
                raise PreparationError("uncapped_artificial_cut", f"artificial N-side cut at {block[0]}")
            selected_atom = source.atom(block[0], "N")
            removed_atom = source.atom(previous, "C")
            if selected_atom is None or removed_atom is None:
                raise PreparationError("artificial_cut_cannot_be_capped", f"cannot cap N-side cut at {block[0]}")
            direction = unit(removed_atom.coord - selected_atom.coord)
            if direction is None:
                raise PreparationError("artificial_cut_cannot_be_capped", f"zero-length N-side cut at {block[0]}")
            record = FragmentAtom(
                element="H",
                coord=selected_atom.coord + LINK_BOND_LENGTHS["N"] * direction,
                atom_name="HLN",
                chain_id=block[0][0],
                seq_id=block[0][1],
                insertion_code=block[0][2],
                resname=source.residues[block[0]]["resname"],
                source_serial=None,
                source_kind="link_hydrogen",
                label=f"LINK:N:{block[0][0]}:{block[0][1]}{block[0][2]}",
            )
            index = _add_atom(rw, rw.GetConformer(), record)
            _add_bond_if_absent(rw, mapped_index(block[0], "N"), index, 1.0)
            capped_records.append(record)
            cap_links.append({
                "boundary_type": "artificial_N_side_peptide_cut",
                "affected_heavy_atom_labels": [initial_records[mapped_index(block[0], "N")].label],
                "selected_atom": f"{block[0]}:N",
                "removed_external_atom": f"{previous}:C",
                "link_hydrogen": record.label,
                "bond_length_A": LINK_BOND_LENGTHS["N"],
                "treated_as_biological_terminus": False,
            })
        if following is not None:
            if not preparation.hydrogen_caps:
                raise PreparationError("uncapped_artificial_cut", f"artificial C-side cut at {block[-1]}")
            selected_atom = source.atom(block[-1], "C")
            removed_atom = source.atom(following, "N")
            if selected_atom is None or removed_atom is None:
                raise PreparationError("artificial_cut_cannot_be_capped", f"cannot cap C-side cut at {block[-1]}")
            direction = unit(removed_atom.coord - selected_atom.coord)
            if direction is None:
                raise PreparationError("artificial_cut_cannot_be_capped", f"zero-length C-side cut at {block[-1]}")
            record = FragmentAtom(
                element="H",
                coord=selected_atom.coord + LINK_BOND_LENGTHS["C"] * direction,
                atom_name="HLC",
                chain_id=block[-1][0],
                seq_id=block[-1][1],
                insertion_code=block[-1][2],
                resname=source.residues[block[-1]]["resname"],
                source_serial=None,
                source_kind="link_hydrogen",
                label=f"LINK:C:{block[-1][0]}:{block[-1][1]}{block[-1][2]}",
            )
            index = _add_atom(rw, rw.GetConformer(), record)
            _add_bond_if_absent(rw, mapped_index(block[-1], "C"), index, 1.0)
            capped_records.append(record)
            cap_links.append({
                "boundary_type": "artificial_C_side_peptide_cut",
                "affected_heavy_atom_labels": [initial_records[mapped_index(block[-1], name)].label for name in ("C", "O")],
                "selected_atom": f"{block[-1]}:C",
                "removed_external_atom": f"{following}:N",
                "link_hydrogen": record.label,
                "bond_length_A": LINK_BOND_LENGTHS["C"],
                "treated_as_biological_terminus": False,
            })

    _set_aromatic_atoms(rw)
    protonation_states, histidine_manual = _configure_protonation(
        rw,
        atom_index,
        source,
        selected_keys,
        ph,
        terminal_states,
    )

    try:
        molecule = rw.GetMol()
        Chem.SanitizeMol(molecule)
        molecule = Chem.AddHs(molecule, addCoords=True)
        Chem.SanitizeMol(molecule)
    except Exception as exc:
        raise PreparationError("invalid_valence", f"local fragment valence/hydrogenation failed: {exc}") from exc

    environment_report = {
        **environment,
        "included_direct_waters": [list(key) for key in sorted(included_water_keys)],
        "automatic_ions_included": False,
        "automatic_charged_residues_included": False,
    }
    cap_report = {
        "policy": "hydrogen link caps at artificial peptide cuts",
        "peptide_link_caps": cap_links,
        "artificial_cuts_are_charged_termini": False,
    }
    protonation_report = {
        "pH": ph,
        "pH_source": ph_source,
        "method": "local RDKit valence hydrogenation after fragment extraction",
        "residue_states": protonation_states,
        "manual_review_required": ph_manual_review or histidine_manual,
        "manual_review_reasons": [
            reason for condition, reason in (
                (ph_manual_review, "crystallization pH missing; fallback pH used"),
                (histidine_manual, "neutral histidine tautomer assigned as HIE by default"),
            ) if condition
        ],
    }
    protonation_report["stage2_quality"] = quality_report
    return BuildResult(
        molecule=molecule,
        initial_atoms=initial_records,
        capped_atoms=capped_records,
        cap_report=cap_report,
        protonation_report=protonation_report,
        environment_report=environment_report,
        source_heavy_coordinates=source_heavy_coordinates,
    )


def _prepared_atom_records(molecule: Chem.Mol, capped_atoms: Sequence[FragmentAtom]) -> list[FragmentAtom]:
    conformer = molecule.GetConformer()
    records = list(capped_atoms)
    hydrogen_counts: Counter[str] = Counter()
    residue_hydrogen_counts: Counter[tuple[str, int, str]] = Counter()
    for index in range(len(capped_atoms), molecule.GetNumAtoms()):
        atom = molecule.GetAtomWithIdx(index)
        position = conformer.GetAtomPosition(index)
        neighbors = list(atom.GetNeighbors())
        if atom.GetSymbol().upper() != "H" or not neighbors:
            raise PreparationError("hydrogen_mapping_failed", f"unexpected added atom at RDKit index {index}")
        neighbor_index = neighbors[0].GetIdx()
        neighbor = records[neighbor_index]
        base = neighbor.label or f"{neighbor.chain_id}:{neighbor.resname}{neighbor.seq_id}:{neighbor.atom_name}"
        hydrogen_counts[base] += 1
        residue_key = (neighbor.chain_id, neighbor.seq_id, neighbor.insertion_code)
        residue_hydrogen_counts[residue_key] += 1
        records.append(FragmentAtom(
            element="H",
            coord=np.array([position.x, position.y, position.z], dtype=float),
            atom_name=f"H{residue_hydrogen_counts[residue_key]}",
            chain_id=neighbor.chain_id,
            seq_id=neighbor.seq_id,
            insertion_code=neighbor.insertion_code,
            resname=neighbor.resname,
            source_serial=None,
            source_kind="added_hydrogen",
            label=f"{base}:H{hydrogen_counts[base]}",
        ))
    # Update coordinates of capped atoms from the final molecule without changing heavy atoms.
    for index, record in enumerate(records[:len(capped_atoms)]):
        position = conformer.GetAtomPosition(index)
        record.coord = np.array([position.x, position.y, position.z], dtype=float)
    return records


def _write_pdb(path: Path, atoms: Sequence[FragmentAtom]) -> None:
    lines: list[str] = []
    for index, atom in enumerate(atoms, 1):
        record_name = "HETATM" if atom.resname in WATER_NAMES else "ATOM  "
        chain = (atom.chain_id or "A")[:1]
        insertion = (atom.insertion_code or " ")[:1]
        atom_name = atom.atom_name[:4]
        x, y, z = map(float, atom.coord)
        lines.append(
            f"{record_name}{index:5d} {atom_name:^4s} {atom.resname[:3]:>3s} {chain}{atom.seq_id:4d}{insertion}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {atom.element:>2s}"
        )
    path.write_text("\n".join(lines + ["END"]) + "\n", encoding="utf-8")


def validate_prepared_fragment(
    molecule: Chem.Mol,
    atoms: Sequence[FragmentAtom],
    source_coordinates: Mapping[str, np.ndarray],
    quality: QualityConfig,
    allow_oversize: bool,
) -> dict[str, Any]:
    charge = int(Chem.GetFormalCharge(molecule))
    electrons = sum(ATOMIC_NUMBERS.get(atom.element.upper(), 0) for atom in atoms) - charge
    if electrons % 2:
        raise PreparationError(
            "charge_electron_inconsistency",
            f"formal charge {charge} gives odd electron count {electrons} for an expected closed-shell fragment",
        )
    multiplicity = 1
    mismatches: list[dict[str, Any]] = []
    for atom in atoms:
        if atom.source_kind not in {"source_heavy", "source_heavy_water"}:
            continue
        original = source_coordinates.get(atom.label)
        if original is None:
            continue
        displacement = distance(atom.coord, original)
        if displacement > 1.0e-6:
            mismatches.append({"label": atom.label, "displacement_A": displacement})
    if mismatches:
        raise PreparationError("heavy_atom_coordinate_changed", f"{len(mismatches)} source heavy atoms moved")

    adjacency = {index: {neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors()} for index in range(molecule.GetNumAtoms())}
    coords = np.vstack([atom.coord for atom in atoms])
    clashes: list[dict[str, Any]] = []
    for index, atom in enumerate(atoms):
        if atom.element != "H":
            continue
        for other_index, other in enumerate(atoms):
            if other_index == index or other_index in adjacency[index]:
                continue
            if distance(atom.coord, other.coord) < 0.55:
                clashes.append({"hydrogen": atom.label, "other": other.label, "distance_A": distance(atom.coord, other.coord)})
                break
    if clashes:
        raise PreparationError("hydrogen_clash", f"{len(clashes)} severe nonbonded hydrogen clashes detected")

    atom_count = len(atoms)
    if atom_count > quality.maximum_atoms and not allow_oversize:
        raise PreparationError(
            "atom_count_above_hard_limit",
            f"prepared fragment has {atom_count} atoms; hard maximum is {quality.maximum_atoms}",
        )
    return {
        "atom_count": atom_count,
        "heavy_atom_count": sum(atom.element != "H" for atom in atoms),
        "hydrogen_count": sum(atom.element == "H" for atom in atoms),
        "preferred_atom_limit": quality.preferred_atoms,
        "hard_atom_limit": quality.maximum_atoms,
        "large_fragment": atom_count > quality.preferred_atoms,
        "oversize_override_used": atom_count > quality.maximum_atoms and allow_oversize,
        "formal_charge": charge,
        "multiplicity": multiplicity,
        "electron_count": electrons,
        "electron_parity": "even",
        "heavy_atom_coordinate_mismatches": mismatches,
        "hydrogen_clashes": clashes,
        "validation_pass": True,
    }


def write_fragment_outputs(
    folder: Path,
    site: Mapping[str, Any],
    capped_atoms: Sequence[FragmentAtom],
    prepared_atoms: Sequence[FragmentAtom],
    validation: Mapping[str, Any],
    reports: Mapping[str, Any],
) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    _write_pdb(folder / "fragment_capped_heavy.pdb", capped_atoms)
    _write_pdb(folder / "fragment_prepared.pdb", prepared_atoms)

    xyz_lines = [str(len(prepared_atoms)), str(site["site_id"])]
    gaussian_lines: list[str] = []
    atom_map: list[dict[str, Any]] = []
    freeze_mask: list[dict[str, Any]] = []
    boundary_types = {
        label: link["boundary_type"]
        for link in reports.get("cap_report", {}).get("peptide_link_caps", [])
        for label in link.get("affected_heavy_atom_labels", [])
    }
    for index, atom in enumerate(prepared_atoms, 1):
        x, y, z = map(float, atom.coord)
        freeze_flag = 0 if atom.element == "H" else -1
        label = atom.label or f"{atom.chain_id}:{atom.resname}{atom.seq_id}:{atom.atom_name}"
        xyz_lines.append(f"{atom.element:<2s} {x: .10f} {y: .10f} {z: .10f}")
        gaussian_lines.append(f"{atom.element:<2s} {freeze_flag:2d} {x: .10f} {y: .10f} {z: .10f} ! {label}")
        atom_map.append({
            "fragment_index": index,
            "element": atom.element,
            "atom_name": atom.atom_name,
            "label": label,
            "chain_id": atom.chain_id,
            "seq_id": atom.seq_id,
            "insertion_code": atom.insertion_code,
            "resname": atom.resname,
            "source_serial": atom.source_serial,
            "source_kind": atom.source_kind,
            "x": x,
            "y": y,
            "z": z,
            "boundary_affected": atom.label in boundary_types,
            "boundary_type": boundary_types.get(atom.label, ""),
        })
        freeze_mask.append({
            "fragment_index": index,
            "label": label,
            "element": atom.element,
            "freeze_flag": freeze_flag,
            "movable": atom.element == "H",
        })
    (folder / "fragment_prepared.xyz").write_text("\n".join(xyz_lines) + "\n", encoding="utf-8")
    (folder / "gaussian_coordinates.txt").write_text("\n".join(gaussian_lines) + "\n", encoding="utf-8")
    write_csv(folder / "atom_map.csv", atom_map)
    write_csv(folder / "freeze_mask.csv", freeze_mask)

    site_metadata = dict(site)
    site_metadata["large_fragment"] = bool(validation["large_fragment"])
    write_json(folder / "site_metadata.json", site_metadata)
    for name, report in reports.items():
        write_json(folder / f"{name}.json", report)


def prepare_site(
    site: Mapping[str, Any],
    source: SourceIndex,
    output_root: Path,
    quality: QualityConfig,
    preparation: PreparationConfig,
) -> dict[str, Any]:
    selected_keys = source.selected_keys(site.get("residue_selectors") or [])
    folder = output_root / str(site["branch"]) / str(site["motif_type"]) / str(site["site_id"])
    build = build_local_fragment(source, site, quality, preparation)
    prepared_atoms = _prepared_atom_records(build.molecule, build.capped_atoms)
    validation = validate_prepared_fragment(
        build.molecule,
        prepared_atoms,
        build.source_heavy_coordinates,
        quality,
        preparation.allow_oversize,
    )
    source.write_subset_cif(selected_keys, folder / "source_fragment.cif")
    charge_report = {
        "formal_charge": validation["formal_charge"],
        "multiplicity": validation["multiplicity"],
        "electron_count": validation["electron_count"],
        "neutrality_forced": False,
    }
    reports = {
        "protonation_report": build.protonation_report,
        "charge_multiplicity": charge_report,
        "cap_report": build.cap_report,
        "environment_report": build.environment_report,
        "validation_report": validation,
    }
    write_fragment_outputs(folder, site, build.capped_atoms, prepared_atoms, validation, reports)
    return {
        "site_id": site["site_id"],
        "pdb_id": site["pdb_id"],
        "branch": site["branch"],
        "motif_type": site["motif_type"],
        "folder": str(folder),
        "atom_count": validation["atom_count"],
        "large_fragment": validation["large_fragment"],
        "formal_charge": validation["formal_charge"],
        "multiplicity": validation["multiplicity"],
    }


def _geometry_vector(site: Mapping[str, Any]) -> tuple[float, ...]:
    geometry = site.get("geometry") or {}
    primary = geometry.get("primary") or {}
    values = (
        safe_float(primary.get("o_donor_c_acceptor_distance_A")),
        safe_float(primary.get("o_donor_c_acceptor_o_acceptor_angle_deg")),
        safe_float(primary.get("donor_phi_deg")),
        safe_float(primary.get("donor_psi_deg")),
        safe_float(geometry.get("strand_axis_dot")),
    )
    return tuple(0.0 if value is None else value for value in values)


def select_sites(sites: Sequence[Mapping[str, Any]], max_per_motif: int) -> list[dict[str, Any]]:
    if max_per_motif == 0:
        return [dict(site) for site in sites]
    selected: list[dict[str, Any]] = []
    for motif_type in sorted({str(site["motif_type"]) for site in sites}):
        group = sorted(
            (dict(site) for site in sites if site["motif_type"] == motif_type),
            key=lambda site: (str(site.get("pdb_id")), str(site.get("site_id"))),
        )
        if len(group) <= max_per_motif:
            selected.extend(group)
            continue
        # Deterministic diversity: first prefer unique parent structures, then farthest geometry.
        chosen: list[dict[str, Any]] = []
        used_pdb: set[str] = set()
        for site in group:
            pdb_id = str(site.get("pdb_id"))
            if pdb_id not in used_pdb:
                chosen.append(site)
                used_pdb.add(pdb_id)
            if len(chosen) == max_per_motif:
                break
        remaining = [site for site in group if site not in chosen]
        while len(chosen) < max_per_motif and remaining:
            def min_distance(candidate: Mapping[str, Any]) -> float:
                vector = np.asarray(_geometry_vector(candidate), dtype=float)
                return min(
                    float(np.linalg.norm(vector - np.asarray(_geometry_vector(existing), dtype=float)))
                    for existing in chosen
                ) if chosen else 0.0
            next_site = max(remaining, key=lambda site: (min_distance(site), str(site["site_id"])))
            chosen.append(next_site)
            remaining.remove(next_site)
        selected.extend(chosen[:max_per_motif])
    return selected


def run_stage3(
    sites_path: Path,
    output: Path,
    selection: SelectionConfig,
    quality: QualityConfig,
    preparation: PreparationConfig,
) -> dict[str, Any]:
    all_sites = read_jsonl(sites_path)
    sites = [site for site in all_sites if site.get("branch") == "protein" and site.get("motif_type") in APPROVED_MOTIFS]
    selected = select_sites(sites, selection.max_per_motif)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "selected_sites.csv", selected)

    source_cache: dict[str, SourceIndex] = {}
    prepared: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    traceback_blocks: list[str] = []
    logged_signatures: set[tuple[str, str]] = set()

    for index, site in enumerate(selected, 1):
        try:
            source_path = str(site["source_mmcif"])
            if source_path not in source_cache:
                source_cache[source_path] = SourceIndex(Path(source_path))
            prepared.append(
                prepare_site(
                    site,
                    source_cache[source_path],
                    output / "fragments",
                    quality,
                    preparation,
                )
            )
        except Exception as exc:
            code = exc.code if isinstance(exc, PreparationError) else "unexpected_error"
            error = {
                "site_id": site.get("site_id"),
                "pdb_id": site.get("pdb_id"),
                "motif_type": site.get("motif_type"),
                "reason": code,
                "error": str(exc),
            }
            errors.append(error)
            signature = (str(site.get("pdb_id")), code)
            if signature not in logged_signatures:
                logging.error(
                    "%s | %s | %s",
                    site.get("site_id"),
                    site.get("motif_type"),
                    str(exc),
                )
                logged_signatures.add(signature)
            traceback_blocks.append(
                f"\n===== {site.get('site_id')} | {site.get('pdb_id')} | {site.get('motif_type')} =====\n{traceback.format_exc()}"
            )
        if index % 10 == 0 or index == len(selected):
            logging.info("Stage 3: %d/%d sites", index, len(selected))

    write_csv(output / "prepared_fragments.csv", prepared)
    write_json(output / "errors.json", errors)
    (output / "debug_tracebacks.log").write_text("".join(traceback_blocks), encoding="utf-8")
    failure_counts = Counter(error["reason"] for error in errors)
    prepared_counts = Counter(item["motif_type"] for item in prepared)
    summary = {
        "sites_received": len(all_sites),
        "approved_sites_received": len(sites),
        "sites_selected": len(selected),
        "sites_prepared": len(prepared),
        "sites_failed": len(errors),
        "failures_by_reason": dict(sorted(failure_counts.items())),
        "prepared_by_motif": dict(sorted(prepared_counts.items())),
        "source_structures_loaded": len(source_cache),
    }
    write_json(output / "summary.json", summary)
    return summary
