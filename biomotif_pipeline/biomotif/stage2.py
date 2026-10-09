from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import is_aa

from .common import (
    PROTEIN_MOTIF_ALLOWLIST,
    angle_deg,
    circular_rms,
    dihedral_deg,
    distance,
    materialized_text_file,
    safe_float,
    stable_id,
    unit,
)

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HYP": "P", "SEC": "U",
}
HELIX_CODES = {"H": "alpha_helix", "G": "3_10_helix", "I": "pi_helix"}
HELIX_MIN_LENGTH = {"H": 4, "G": 3, "I": 5}
TURN_TARGETS = {
    "beta_turn_I": (-60.0, -30.0, -90.0, 0.0),
    "beta_turn_II": (-60.0, 120.0, 80.0, 0.0),
    "beta_turn_I_prime": (60.0, 30.0, 90.0, 0.0),
    "beta_turn_II_prime": (60.0, -120.0, -80.0, 0.0),
}
GAMMA_TARGETS = {
    "gamma_turn_classic": (70.0, -60.0),
    "gamma_turn_inverse": (-70.0, 60.0),
}


@dataclass
class Stage2Config:
    peptide_bond_max: float = 1.90
    npi_distance_max: float = 3.20
    npi_angle_min: float = 99.0
    npi_angle_max: float = 119.0
    npi_border_distance: float = 0.15
    npi_border_angle: float = 10.0
    beta_turn_ca_max: float = 7.0
    beta_turn_tolerance: float = 35.0
    gamma_turn_ca_max: float = 6.0
    gamma_turn_on_max: float = 3.5
    gamma_turn_tolerance: float = 35.0
    collagen_min_length: int = 6
    collagen_gly_fraction: float = 0.60
    local_min_occupancy: float = 0.90


@dataclass
class AtomChoice:
    name: str
    element: str
    coord: np.ndarray
    occupancy: float | None
    bfactor: float | None
    altloc: str
    had_alternates: bool


@dataclass
class ProteinResidue:
    pdb_id: str
    chain_id: str
    seq_id: int
    insertion_code: str
    resname: str
    one_letter: str
    chain_index: int
    residue_obj: Any = field(repr=False)
    atoms: dict[str, AtomChoice] = field(default_factory=dict, repr=False)
    phi: float | None = None
    psi: float | None = None
    omega: float | None = None
    dssp_code: str = "-"
    dssp_index: int | None = None
    dssp_hbond_offsets: list[int] = field(default_factory=list)
    dssp_hbonds: list[dict[str, Any]] = field(default_factory=list)
    relative_asa: float | None = None
    connected_prev: bool = False
    connected_next: bool = False
    min_occupancy: float | None = None
    mean_bfactor: float | None = None
    max_bfactor: float | None = None
    has_altloc: bool = False

    @property
    def key(self) -> tuple[str, int, str]:
        return self.chain_id, self.seq_id, self.insertion_code

    @property
    def selector(self) -> dict[str, Any]:
        return {
            "chain_id": self.chain_id,
            "seq_id": self.seq_id,
            "insertion_code": self.insertion_code,
            "resname": self.resname,
        }

    @property
    def label(self) -> str:
        return f"{self.chain_id}:{self.resname}{self.seq_id}{self.insertion_code}"


def _atom_candidates(atom: Any) -> list[Any]:
    if getattr(atom, "is_disordered", lambda: False)():
        getter = getattr(atom, "disordered_get_list", None)
        if getter:
            values = list(getter())
            if values:
                return values
    return [atom]


def choose_atom(residue: Any, name: str) -> AtomChoice | None:
    if name not in residue:
        return None
    candidates = _atom_candidates(residue[name])
    selected = max(
        candidates,
        key=lambda atom: (
            float(atom.get_occupancy() or 0.0),
            str(atom.get_altloc()).strip() in {"", "A", "1"},
        ),
    )
    element = str(getattr(selected, "element", "") or name[0]).strip().upper()
    return AtomChoice(
        name=name,
        element=element,
        coord=np.asarray(selected.get_coord(), dtype=float),
        occupancy=safe_float(selected.get_occupancy()),
        bfactor=safe_float(selected.get_bfactor()),
        altloc=str(selected.get_altloc()).strip(),
        had_alternates=len(candidates) > 1,
    )


def residue_quality(atoms: Iterable[AtomChoice]) -> tuple[float | None, float | None, float | None, bool]:
    values = list(atoms)
    occupancies = [atom.occupancy for atom in values if atom.occupancy is not None]
    bfactors = [atom.bfactor for atom in values if atom.bfactor is not None]
    return (
        min(occupancies) if occupancies else None,
        statistics.fmean(bfactors) if bfactors else None,
        max(bfactors) if bfactors else None,
        any(atom.had_alternates or bool(atom.altloc) for atom in values),
    )


def collect_protein_residues(
    pdb_id: str,
    model: Any,
    cfg: Stage2Config,
) -> dict[str, list[ProteinResidue]]:
    chains: dict[str, list[ProteinResidue]] = {}
    for chain in model:
        residues: list[ProteinResidue] = []
        for residue in chain:
            if not is_aa(residue, standard=False):
                continue
            resname = str(residue.get_resname()).strip().upper()
            atoms: dict[str, AtomChoice] = {}
            for atom_wrapper in residue:
                atom_name = str(atom_wrapper.get_name()).strip()
                atom = choose_atom(residue, atom_name)
                if atom is not None and atom.element != "H":
                    atoms[atom_name] = atom
            min_occ, mean_b, max_b, altloc = residue_quality(atoms.values())
            residues.append(
                ProteinResidue(
                    pdb_id=pdb_id,
                    chain_id=str(chain.id),
                    seq_id=int(residue.id[1]),
                    insertion_code=str(residue.id[2]).strip(),
                    resname=resname,
                    one_letter=AA3_TO_1.get(resname, "X"),
                    chain_index=len(residues),
                    residue_obj=residue,
                    atoms=atoms,
                    min_occupancy=min_occ,
                    mean_bfactor=mean_b,
                    max_bfactor=max_b,
                    has_altloc=altloc,
                )
            )
        for index, residue in enumerate(residues):
            if index > 0 and "C" in residues[index - 1].atoms and "N" in residue.atoms:
                residue.connected_prev = (
                    distance(residues[index - 1].atoms["C"].coord, residue.atoms["N"].coord)
                    <= cfg.peptide_bond_max
                )
            if index + 1 < len(residues) and "C" in residue.atoms and "N" in residues[index + 1].atoms:
                residue.connected_next = (
                    distance(residue.atoms["C"].coord, residues[index + 1].atoms["N"].coord)
                    <= cfg.peptide_bond_max
                )
        calculate_torsions(residues)
        if residues:
            chains[str(chain.id)] = residues
    return chains


def calculate_torsions(residues: Sequence[ProteinResidue]) -> None:
    for index, residue in enumerate(residues):
        atoms = residue.atoms
        if index > 0 and residue.connected_prev:
            previous = residues[index - 1].atoms
            if "C" in previous and all(name in atoms for name in ("N", "CA", "C")):
                residue.phi = dihedral_deg(
                    previous["C"].coord,
                    atoms["N"].coord,
                    atoms["CA"].coord,
                    atoms["C"].coord,
                )
            if all(name in previous for name in ("CA", "C")) and all(name in atoms for name in ("N", "CA")):
                residue.omega = dihedral_deg(
                    previous["CA"].coord,
                    previous["C"].coord,
                    atoms["N"].coord,
                    atoms["CA"].coord,
                )
        if index + 1 < len(residues) and residue.connected_next:
            following = residues[index + 1].atoms
            if all(name in atoms for name in ("N", "CA", "C")) and "N" in following:
                residue.psi = dihedral_deg(
                    atoms["N"].coord,
                    atoms["CA"].coord,
                    atoms["C"].coord,
                    following["N"].coord,
                )


def assign_dssp(
    model: Any,
    cif_path: Path,
    chains: dict[str, list[ProteinResidue]],
    executable: str,
) -> None:
    dssp = DSSP(model, str(cif_path), dssp=executable, file_type="MMCIF")
    lookup = {residue.key: residue for residues in chains.values() for residue in residues}
    for key in dssp.keys():
        chain_id, residue_id = key
        _, seq_id, insertion = residue_id
        residue = lookup.get((str(chain_id), int(seq_id), str(insertion).strip()))
        if residue is None:
            continue
        data = dssp[key]
        residue.dssp_index = int(data[0])
        residue.dssp_code = str(data[2]).strip() or "-"
        residue.relative_asa = safe_float(data[3])
        offsets: list[int] = []
        hbonds: list[dict[str, Any]] = []
        for position in (6, 8, 10, 12):
            try:
                offset = int(data[position])
            except (TypeError, ValueError, IndexError):
                continue
            if offset:
                offsets.append(offset)
                energy = safe_float(data[position + 1]) if position + 1 < len(data) else None
                hbonds.append({
                    "partner_offset": offset,
                    "role": "donor" if position in (6, 10) else "acceptor",
                    "energy_kcal_mol": energy,
                })
        residue.dssp_hbond_offsets = offsets
        residue.dssp_hbonds = hbonds


def contiguous(residues: Sequence[ProteinResidue], code: str) -> list[list[ProteinResidue]]:
    segments: list[list[ProteinResidue]] = []
    index = 0
    while index < len(residues):
        if residues[index].dssp_code != code:
            index += 1
            continue
        end = index + 1
        while end < len(residues) and residues[end].dssp_code == code and residues[end].connected_prev:
            end += 1
        segments.append(list(residues[index:end]))
        index = end
    return segments


def geometry_label(distance_value: float | None, angle_value: float | None, cfg: Stage2Config) -> str:
    if distance_value is None or angle_value is None:
        return "outside_window"
    if (
        distance_value <= cfg.npi_distance_max
        and cfg.npi_angle_min <= angle_value <= cfg.npi_angle_max
    ):
        return "strong"
    if (
        distance_value <= cfg.npi_distance_max + cfg.npi_border_distance
        and cfg.npi_angle_min - cfg.npi_border_angle
        <= angle_value
        <= cfg.npi_angle_max + cfg.npi_border_angle
    ):
        return "borderline"
    return "outside_window"


# Compatibility name retained for existing notebooks/tests.
def geometry_bucket(distance_value: float | None, angle_value: float | None, cfg: Stage2Config) -> str:
    return geometry_label(distance_value, angle_value, cfg)


def carbonyl_pair_geometry(
    donor: ProteinResidue,
    acceptor: ProteinResidue,
    cfg: Stage2Config,
) -> dict[str, Any] | None:
    if not donor.connected_next:
        return None
    if not all(name in donor.atoms for name in ("O", "C")):
        return None
    if not all(name in acceptor.atoms for name in ("C", "O")):
        return None
    distance_value = distance(donor.atoms["O"].coord, acceptor.atoms["C"].coord)
    angle_value = angle_deg(
        donor.atoms["O"].coord,
        acceptor.atoms["C"].coord,
        acceptor.atoms["O"].coord,
    )
    relevant_atoms = [donor.atoms["O"], acceptor.atoms["C"], acceptor.atoms["O"]]
    occupancies = [atom.occupancy for atom in relevant_atoms if atom.occupancy is not None]
    bfactors = [atom.bfactor for atom in relevant_atoms if atom.bfactor is not None]
    return {
        "donor": donor.selector,
        "acceptor": acceptor.selector,
        "o_donor_c_acceptor_distance_A": distance_value,
        "o_donor_c_acceptor_o_acceptor_angle_deg": angle_value,
        "donor_phi_deg": donor.phi,
        "donor_psi_deg": donor.psi,
        "acceptor_phi_deg": acceptor.phi,
        "acceptor_psi_deg": acceptor.psi,
        "min_occupancy": min(occupancies) if occupancies else None,
        "mean_bfactor": statistics.fmean(bfactors) if bfactors else None,
        "has_altloc": any(atom.had_alternates or bool(atom.altloc) for atom in relevant_atoms),
        "geometry_label": geometry_label(distance_value, angle_value, cfg),
    }


def carbonyl_geometries(
    residues: Sequence[ProteinResidue],
    cfg: Stage2Config,
) -> list[dict[str, Any]]:
    geometries: list[dict[str, Any]] = []
    by_chain: dict[str, list[ProteinResidue]] = {}
    for residue in residues:
        by_chain.setdefault(residue.chain_id, []).append(residue)
    for chain_residues in by_chain.values():
        ordered = sorted(chain_residues, key=lambda residue: residue.chain_index)
        for left, right in zip(ordered, ordered[1:]):
            if right.chain_index != left.chain_index + 1:
                continue
            result = carbonyl_pair_geometry(left, right, cfg)
            if result is not None:
                geometries.append(result)
    return geometries


def quality_summary(residues: Sequence[ProteinResidue], cfg: Stage2Config) -> dict[str, Any]:
    occupancies = [residue.min_occupancy for residue in residues if residue.min_occupancy is not None]
    bfactors = [residue.mean_bfactor for residue in residues if residue.mean_bfactor is not None]
    min_occupancy = min(occupancies) if occupancies else None
    has_altloc = any(residue.has_altloc for residue in residues)
    return {
        "min_occupancy": min_occupancy,
        "mean_local_bfactor": statistics.fmean(bfactors) if bfactors else None,
        "has_altloc": has_altloc,
        "stage2_quality_pass": (
            min_occupancy is not None
            and min_occupancy >= cfg.local_min_occupancy
            and not has_altloc
        ),
    }


def _central_atoms_from_geometry(geometry: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    if not geometry:
        return []
    donor = geometry.get("donor") or {}
    acceptor = geometry.get("acceptor") or {}
    return [
        {
            "chain_id": donor.get("chain_id"),
            "seq_id": donor.get("seq_id"),
            "insertion_code": donor.get("insertion_code", ""),
            "atom_name": "O",
        },
        {
            "chain_id": acceptor.get("chain_id"),
            "seq_id": acceptor.get("seq_id"),
            "insertion_code": acceptor.get("insertion_code", ""),
            "atom_name": "C",
        },
        {
            "chain_id": acceptor.get("chain_id"),
            "seq_id": acceptor.get("seq_id"),
            "insertion_code": acceptor.get("insertion_code", ""),
            "atom_name": "O",
        },
    ]


def make_site(
    pdb_id: str,
    motif_type: str,
    residues: Sequence[ProteinResidue],
    cfg: Stage2Config,
    assignment_method: str,
    geometry: Mapping[str, Any] | None = None,
    central_atoms: Sequence[Mapping[str, Any]] | None = None,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if motif_type not in PROTEIN_MOTIF_ALLOWLIST:
        raise ValueError(f"Unapproved motif type: {motif_type}")
    selectors = [residue.selector for residue in residues]
    geometry_dict = dict(geometry or {})
    if "carbonyl_pairs" not in geometry_dict:
        geometry_dict["carbonyl_pairs"] = carbonyl_geometries(residues, cfg)
    if "primary" not in geometry_dict:
        pairs = geometry_dict["carbonyl_pairs"]
        geometry_dict["primary"] = min(
            pairs,
            key=lambda item: item["o_donor_c_acceptor_distance_A"],
        ) if pairs else None
    primary = geometry_dict.get("primary")
    geometry_dict["geometry_label"] = (
        primary.get("geometry_label") if isinstance(primary, Mapping) else "outside_window"
    )
    quality = quality_summary(residues, cfg)
    site_id = stable_id(
        pdb_id,
        motif_type,
        [(selector["chain_id"], selector["seq_id"], selector.get("insertion_code", "")) for selector in selectors],
        geometry_dict.get("primary"),
        prefix="SITE_",
    )
    record: dict[str, Any] = {
        "site_id": site_id,
        "pdb_id": pdb_id,
        "branch": "protein",
        "motif_type": motif_type,
        "chain_ids": sorted({residue.chain_id for residue in residues}),
        "residue_selectors": selectors,
        "central_interaction_atoms": list(central_atoms or _central_atoms_from_geometry(primary)),
        "assignment_method": assignment_method,
        "geometry": geometry_dict,
        "quality": quality,
        "manual_review_required": False,
    }
    if extra:
        record.update(extra)
    return record


def _four_residue_window(
    chain: Sequence[ProteinResidue],
    pair_start_index: int,
) -> list[ProteinResidue] | None:
    if len(chain) < 4:
        return None
    candidates: list[tuple[int, list[ProteinResidue]]] = []
    for start in range(max(0, pair_start_index - 2), min(pair_start_index + 1, len(chain) - 4) + 1):
        window = list(chain[start:start + 4])
        if not (start <= pair_start_index and pair_start_index + 1 <= start + 3):
            continue
        if all(window[index].connected_next for index in range(3)):
            score = abs((start + 1.5) - (pair_start_index + 0.5))
            candidates.append((int(score * 1000), window))
    if not candidates:
        return None
    return sorted(candidates, key=lambda item: (item[0], item[1][0].chain_index))[0][1]


def annotate_helix_sites(
    pdb_id: str,
    chains: Mapping[str, Sequence[ProteinResidue]],
    cfg: Stage2Config,
) -> list[dict[str, Any]]:
    sites: list[dict[str, Any]] = []
    seen: set[str] = set()
    for chain in chains.values():
        for code, motif_type in HELIX_CODES.items():
            for segment in contiguous(chain, code):
                if len(segment) < HELIX_MIN_LENGTH[code]:
                    continue
                for donor, acceptor in zip(segment, segment[1:]):
                    if acceptor.chain_index != donor.chain_index + 1:
                        continue
                    window = _four_residue_window(chain, donor.chain_index)
                    if window is None:
                        continue
                    primary = carbonyl_pair_geometry(donor, acceptor, cfg)
                    geometry = {
                        "primary": primary,
                        "carbonyl_pairs": carbonyl_geometries(window, cfg),
                        "dssp_code": code,
                    }
                    site = make_site(
                        pdb_id,
                        motif_type,
                        window,
                        cfg,
                        assignment_method="DSSP_four_residue_window",
                        geometry=geometry,
                    )
                    if site["site_id"] not in seen:
                        sites.append(site)
                        seen.add(site["site_id"])
    return sites


def annotate_turn_sites(
    pdb_id: str,
    chains: Mapping[str, Sequence[ProteinResidue]],
    cfg: Stage2Config,
) -> list[dict[str, Any]]:
    sites: list[dict[str, Any]] = []
    for chain in chains.values():
        for index in range(len(chain) - 3):
            window = list(chain[index:index + 4])
            if not all(window[position].connected_next for position in range(3)):
                continue
            if not all("CA" in residue.atoms for residue in window):
                continue
            if any(residue.dssp_code in HELIX_CODES for residue in window[1:3]):
                continue
            ca_distance = distance(window[0].atoms["CA"].coord, window[3].atoms["CA"].coord)
            if ca_distance > cfg.beta_turn_ca_max:
                continue
            torsions = (window[1].phi, window[1].psi, window[2].phi, window[2].psi)
            if any(value is None for value in torsions):
                continue
            numeric = tuple(float(value) for value in torsions if value is not None)
            matches = []
            for name, target in TURN_TARGETS.items():
                deviations = [abs((value - ideal + 180.0) % 360.0 - 180.0) for value, ideal in zip(numeric, target)]
                ordered = sorted(deviations)
                if ordered[-1] <= 40.0 and ordered[-2] <= 30.0:
                    matches.append((circular_rms(numeric, target), name, deviations))
            if not matches:
                continue
            rms, motif_type, deviations = min(matches)
            if rms > cfg.beta_turn_tolerance:
                continue
            sites.append(
                make_site(
                    pdb_id,
                    motif_type,
                    window,
                    cfg,
                    assignment_method="torsion_geometry",
                    extra={
                        "turn_parameters": {
                            "ca_i_i3_distance_A": ca_distance,
                            "torsion_rms_deg": rms,
                            "canonical_angle_deviations_deg": deviations,
                            "central_dssp_codes": [residue.dssp_code for residue in window[1:3]],
                            "phi2_deg": torsions[0],
                            "psi2_deg": torsions[1],
                            "phi3_deg": torsions[2],
                            "psi3_deg": torsions[3],
                        }
                    },
                )
            )
        for index in range(len(chain) - 2):
            window = list(chain[index:index + 3])
            if not window[0].connected_next or not window[1].connected_next:
                continue
            if not all("CA" in residue.atoms for residue in window):
                continue
            central = window[1]
            if central.phi is None or central.psi is None:
                continue
            ca_distance = distance(window[0].atoms["CA"].coord, window[2].atoms["CA"].coord)
            on_distance = (
                distance(window[0].atoms["O"].coord, window[2].atoms["N"].coord)
                if "O" in window[0].atoms and "N" in window[2].atoms
                else None
            )
            rms, motif_type = min(
                (
                    circular_rms((central.phi, central.psi), target),
                    name,
                )
                for name, target in GAMMA_TARGETS.items()
            )
            if ca_distance > cfg.gamma_turn_ca_max or rms > cfg.gamma_turn_tolerance:
                continue
            if on_distance is not None and on_distance > cfg.gamma_turn_on_max:
                continue
            sites.append(
                make_site(
                    pdb_id,
                    motif_type,
                    window,
                    cfg,
                    assignment_method="torsion_geometry",
                    extra={
                        "turn_parameters": {
                            "ca_i_i2_distance_A": ca_distance,
                            "o_i_n_i2_distance_A": on_distance,
                            "torsion_rms_deg": rms,
                        }
                    },
                )
            )
    return sites


def _local_direction(residue: ProteinResidue, chain: Sequence[ProteinResidue]) -> np.ndarray | None:
    index = residue.chain_index
    if index > 0 and index + 1 < len(chain) and residue.connected_prev and residue.connected_next:
        previous = chain[index - 1]
        following = chain[index + 1]
        if "CA" in previous.atoms and "CA" in following.atoms:
            return unit(following.atoms["CA"].coord - previous.atoms["CA"].coord)
    if "N" in residue.atoms and "C" in residue.atoms:
        return unit(residue.atoms["C"].coord - residue.atoms["N"].coord)
    return None


def _local_peptide_unit(
    residue: ProteinResidue,
    chain: Sequence[ProteinResidue],
) -> list[ProteinResidue]:
    index = residue.chain_index
    if index + 1 < len(chain) and residue.connected_next:
        return [residue, chain[index + 1]]
    if index > 0 and residue.connected_prev:
        return [chain[index - 1], residue]
    return [residue]


def annotate_beta_sheet_sites(
    pdb_id: str,
    chains: Mapping[str, Sequence[ProteinResidue]],
    cfg: Stage2Config,
) -> list[dict[str, Any]]:
    by_dssp = {
        residue.dssp_index: residue
        for chain in chains.values()
        for residue in chain
        if residue.dssp_index is not None
    }
    unique_pairs: dict[tuple[tuple[str, int, str], tuple[str, int, str]], tuple[ProteinResidue, ProteinResidue]] = {}
    for chain in chains.values():
        for residue in chain:
            if residue.dssp_code != "E" or residue.dssp_index is None:
                continue
            for offset in residue.dssp_hbond_offsets:
                partner = by_dssp.get(residue.dssp_index + offset)
                if partner is None or partner.dssp_code != "E" or partner.key == residue.key:
                    continue
                pair_key = tuple(sorted((residue.key, partner.key)))
                unique_pairs[pair_key] = (residue, partner)

    sites: list[dict[str, Any]] = []
    for pair_key in sorted(unique_pairs):
        first, second = unique_pairs[pair_key]
        chain_first = chains[first.chain_id]
        chain_second = chains[second.chain_id]
        first_direction = _local_direction(first, chain_first)
        second_direction = _local_direction(second, chain_second)
        if first_direction is None or second_direction is None:
            continue
        dot = float(np.dot(first_direction, second_direction))
        motif_type = "beta_sheet_parallel" if dot >= 0 else "beta_sheet_antiparallel"
        first_unit = _local_peptide_unit(first, chain_first)
        second_unit = _local_peptide_unit(second, chain_second)
        first_keys = {residue.key for residue in first_unit}
        second_keys = {residue.key for residue in second_unit}
        if first_keys & second_keys:
            continue
        residues: list[ProteinResidue] = []
        for residue in first_unit + second_unit:
            if residue.key not in {item.key for item in residues}:
                residues.append(residue)
        retained_keys = first_keys | second_keys
        directed: dict[tuple[Any, Any], dict[str, Any]] = {}
        for residue in residues:
            for bond in residue.dssp_hbonds:
                partner = by_dssp.get(residue.dssp_index + bond["partner_offset"])
                energy = bond.get("energy_kcal_mol")
                if partner is None or partner.key not in retained_keys or energy is None or energy >= -0.5:
                    continue
                if (residue.key in first_keys) == (partner.key in first_keys):
                    continue
                donor, acceptor = (residue, partner) if bond["role"] == "donor" else (partner, residue)
                if donor.resname == "PRO" or "N" not in donor.atoms or "O" not in acceptor.atoms:
                    continue
                separation = distance(donor.atoms["N"].coord, acceptor.atoms["O"].coord)
                if separation > 3.5:
                    continue
                pair = (donor.key, acceptor.key)
                item = {
                    "donor": {**donor.selector, "atom_name": "N"},
                    "acceptor": {**acceptor.selector, "atom_name": "O"},
                    "energy_kcal_mol": energy,
                    "n_o_distance_A": separation,
                }
                if pair not in directed or energy < directed[pair]["energy_kcal_mol"]:
                    directed[pair] = item
        hbonds = [directed[key] for key in sorted(directed)]
        if not hbonds:
            continue
        central_atoms = []
        seen_atoms: set[tuple[Any, ...]] = set()
        for bond in hbonds:
            for selector in (bond["donor"], bond["acceptor"]):
                key = (selector["chain_id"], selector["seq_id"], selector["insertion_code"], selector["atom_name"])
                if key not in seen_atoms:
                    central_atoms.append(selector)
                    seen_atoms.add(key)
        geometry = {
            "strand_axis_dot": dot,
            "orientation": "parallel" if dot >= 0 else "antiparallel",
            "dssp_candidate_pair": [first.selector, second.selector],
            "dssp_hbond_pair": [hbonds[0]["donor"], hbonds[0]["acceptor"]],
            "dssp_hydrogen_bonds": hbonds,
            "ca_distance_A": (
                distance(first.atoms["CA"].coord, second.atoms["CA"].coord)
                if "CA" in first.atoms and "CA" in second.atoms
                else None
            ),
            "carbonyl_pairs": carbonyl_geometries(residues, cfg),
            "primary": None,
        }
        site = make_site(
            pdb_id,
            motif_type,
            residues,
            cfg,
            assignment_method="DSSP_hydrogen_bond_relationship+strand_direction",
            geometry=geometry,
            central_atoms=central_atoms,
        )
        sites.append(site)
    return sites


def _gly_phase_fraction(segment: Sequence[ProteinResidue]) -> float:
    values: list[float] = []
    for phase in range(3):
        subset = segment[phase::3]
        values.append(
            sum(residue.one_letter == "G" for residue in subset) / len(subset)
            if subset else 0.0
        )
    return max(values, default=0.0)


def annotate_collagen_sites(
    pdb_id: str,
    chains: Mapping[str, Sequence[ProteinResidue]],
    cfg: Stage2Config,
) -> list[dict[str, Any]]:
    segments: list[tuple[list[ProteinResidue], np.ndarray, np.ndarray, float]] = []
    for chain in chains.values():
        compatible = [
            residue.phi is not None
            and residue.psi is not None
            and -105.0 <= residue.phi <= -45.0
            and 90.0 <= residue.psi <= 180.0
            for residue in chain
        ]
        index = 0
        while index < len(chain):
            if not compatible[index]:
                index += 1
                continue
            end = index + 1
            while end < len(chain) and compatible[end] and chain[end].connected_prev:
                end += 1
            segment = list(chain[index:end])
            if len(segment) >= cfg.collagen_min_length:
                gly_fraction = _gly_phase_fraction(segment)
                ca_coords = [residue.atoms["CA"].coord for residue in segment if "CA" in residue.atoms]
                if gly_fraction >= cfg.collagen_gly_fraction and len(ca_coords) >= 2:
                    axis = unit(ca_coords[-1] - ca_coords[0])
                    if axis is not None:
                        segments.append((segment, np.vstack(ca_coords).mean(axis=0), axis, gly_fraction))
            index = end

    sites: list[dict[str, Any]] = []
    for group in combinations(segments, 3):
        chain_ids = {segment[0].chain_id for segment, _, _, _ in group}
        if len(chain_ids) != 3:
            continue
        centroid_distances = [
            distance(group[left][1], group[right][1])
            for left, right in ((0, 1), (0, 2), (1, 2))
        ]
        axis_dots = [
            abs(float(np.dot(group[left][2], group[right][2])))
            for left, right in ((0, 1), (0, 2), (1, 2))
        ]
        if max(centroid_distances) > 12.0 or min(axis_dots) < 0.70:
            continue
        local_residues: list[ProteinResidue] = []
        for segment, _, _, _ in group:
            center = len(segment) // 2
            start = max(0, min(center - 1, len(segment) - 3))
            local_residues.extend(segment[start:start + 3])
        site = make_site(
            pdb_id,
            "collagen_triple_helix",
            local_residues,
            cfg,
            assignment_method="three_chain_Gly-X-Y+PPII_geometry",
            geometry={
                "max_chain_centroid_distance_A": max(centroid_distances),
                "min_abs_axis_dot": min(axis_dots),
                "gly_repeat_fractions": [item[3] for item in group],
                "carbonyl_pairs": carbonyl_geometries(local_residues, cfg),
                "primary": None,
            },
            extra={
                "manual_review_required": True,
                "review_reason": "Conservative automatic collagen triple-helix candidate; verify three-chain register manually.",
            },
        )
        sites.append(site)
    return sites


def protein_residue_rows(chains: Mapping[str, Sequence[ProteinResidue]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for chain in chains.values():
        for residue in chain:
            rows.append({
                "pdb_id": residue.pdb_id,
                "chain_id": residue.chain_id,
                "seq_id": residue.seq_id,
                "insertion_code": residue.insertion_code,
                "resname": residue.resname,
                "phi_deg": residue.phi,
                "psi_deg": residue.psi,
                "omega_deg": residue.omega,
                "dssp_code_internal": residue.dssp_code,
                "relative_asa": residue.relative_asa,
                "min_occupancy": residue.min_occupancy,
                "mean_bfactor": residue.mean_bfactor,
                "has_altloc": residue.has_altloc,
            })
    return rows


def protein_geometry_rows(sites: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for site in sites:
        geometry = site.get("geometry") or {}
        for index, pair in enumerate(geometry.get("carbonyl_pairs") or []):
            rows.append({
                "site_id": site["site_id"],
                "pdb_id": site["pdb_id"],
                "motif_type": site["motif_type"],
                "pair_index": index,
                **dict(pair),
            })
    return rows



def process_structure(
    record: Mapping[str, Any],
    stage1_root: Path,
    cfg: Stage2Config,
    dssp_executable: str | None,
) -> tuple[dict[str, list[dict[str, Any]]], list[dict[str, Any]]]:
    pdb_id = str(record["pdb_id"]).upper()
    if str(record.get("branch") or "protein") != "protein":
        raise RuntimeError(f"Stage 2 accepts protein records only: {pdb_id}")

    source = Path(str(record["mmcif_path"]))
    cif_path = source if source.is_absolute() else stage1_root / source
    cif_path = cif_path.resolve()
    tables: dict[str, list[dict[str, Any]]] = {
        "structures": [],
        "approved_motifs": [],
        "protein_residues": [],
        "protein_geometry": [],
    }

    if dssp_executable is None:
        raise RuntimeError("DSSP executable not found")

    with materialized_text_file(cif_path, suffix=".cif") as materialized:
        parser = MMCIFParser(auth_chains=True, auth_residues=True, QUIET=True)
        structure = parser.get_structure(pdb_id, str(materialized))
        model = next(structure.get_models())
        chains = collect_protein_residues(pdb_id, model, cfg)
        assign_dssp(model, materialized, chains, dssp_executable)
        sites = annotate_helix_sites(pdb_id, chains, cfg)
        sites += annotate_beta_sheet_sites(pdb_id, chains, cfg)
        sites += annotate_turn_sites(pdb_id, chains, cfg)
        sites += annotate_collagen_sites(pdb_id, chains, cfg)

    sites = [site for site in sites if site["motif_type"] in PROTEIN_MOTIF_ALLOWLIST]
    for site in sites:
        site["source_mmcif"] = str(cif_path)
        site["resolution_angstrom"] = safe_float(record.get("resolution_angstrom"))
        site["r_free"] = safe_float(record.get("r_free"))
        site["crystallization_ph"] = safe_float(record.get("crystallization_ph"))
        site["crystallization_ph_min"] = safe_float(record.get("crystallization_ph_min"))
        site["crystallization_ph_max"] = safe_float(record.get("crystallization_ph_max"))

    tables["approved_motifs"] = list(sites)
    tables["protein_residues"] = protein_residue_rows(chains)
    tables["protein_geometry"] = protein_geometry_rows(sites)
    tables["structures"] = [{
        "pdb_id": pdb_id,
        "branch": "protein",
        "source_mmcif": str(cif_path),
        "resolution_angstrom": safe_float(record.get("resolution_angstrom")),
        "r_free": safe_float(record.get("r_free")),
        "crystallization_ph": safe_float(record.get("crystallization_ph")),
        "approved_motif_count": len(sites),
    }]
    return tables, sites
