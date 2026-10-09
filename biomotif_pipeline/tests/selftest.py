#!/usr/bin/env python3
"""Protein-only regression checks; optional read-only replay of existing sites."""
from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import gemmi
import numpy as np
from rdkit import Chem

import biomotif.stage2 as stage2
import biomotif.stage3 as stage3
import select_gaussian_candidates_v2 as selector


def residues(count: int, chain: str = "A") -> list[stage2.ProteinResidue]:
    result = []
    for index in range(count):
        atoms = {}
        for name, xyz in {"CA": (index, 0, 0), "N": (index - .2, 0, 0),
                          "C": (index + .2, 0, 0), "O": (index + .2, 1, 0)}.items():
            atoms[name] = stage2.AtomChoice(name, name[0], np.array(xyz, dtype=float), 1., 5., "", False)
        result.append(stage2.ProteinResidue(
            "TEST", chain, index + 1, "", "GLY", "G", index, None, atoms=atoms,
            connected_prev=index > 0, connected_next=index + 1 < count,
            min_occupancy=1., mean_bfactor=5., phi=-60., psi=120.,
        ))
    return result


def source_record(name: str, number: int, atoms: dict[str, tuple[str, tuple[float, float, float]]],
                  entity_type: gemmi.EntityType = gemmi.EntityType.NonPolymer) -> dict:
    residue = gemmi.Residue()
    residue.name = name
    residue.entity_type = entity_type
    residue.subchain = "A"
    residue.label_seq = number
    return {"resname": name, "gemmi_residue": residue, "atoms": {
        atom: stage3.SourceAtom("A", number, "", name, atom, element, np.array(xyz), 1., 5., "", False, number)
        for atom, (element, xyz) in atoms.items()
    }}


class ProteinTests(unittest.TestCase):
    def test_beta_turn_identity_and_helix_exclusion(self):
        chain = residues(4)
        chain[1].phi, chain[1].psi, chain[2].phi, chain[2].psi = (-60., -30., -90., 0.)
        sites = stage2.annotate_turn_sites("TEST", {"A": chain}, stage2.Stage2Config())
        self.assertIn("beta_turn_I", [site["motif_type"] for site in sites])
        chain[1].dssp_code = "H"
        sites = stage2.annotate_turn_sites("TEST", {"A": chain}, stage2.Stage2Config())
        self.assertFalse(any(site["motif_type"].startswith("beta_turn") for site in sites))

    def test_beta_turn_rejects_two_large_angle_errors(self):
        chain = residues(4)
        chain[1].phi, chain[1].psi, chain[2].phi, chain[2].psi = (-60., -65., -70., 40.)
        sites = stage2.annotate_turn_sites("TEST", {"A": chain}, stage2.Stage2Config())
        self.assertFalse(any(site["motif_type"].startswith("beta_turn") for site in sites))

    def test_gamma_names_from_reviewed_examples(self):
        for torsions, expected in (((-88.911, 72.192), "gamma_turn_inverse"),
                                   ((84.541, -58.081), "gamma_turn_classic")):
            with self.subTest(expected=expected):
                chain = residues(3)
                chain[1].phi, chain[1].psi = torsions
                sites = stage2.annotate_turn_sites("TEST", {"A": chain}, stage2.Stage2Config())
                self.assertEqual([site["motif_type"] for site in sites], [expected])

    def test_dssp_direction_energy_and_retained_sheet_bonds(self):
        chains = {"A": residues(2), "B": residues(2, "B")}
        chains["A"][0].atoms["N"].coord = np.array((0., 0., 0.))
        chains["A"][0].atoms["C"].coord = np.array((1., 0., 0.))
        chains["B"][0].atoms["N"].coord = np.array((10., 0., 0.))
        chains["B"][0].atoms["C"].coord = np.array((9., 0., 0.))
        chains["A"][1].atoms["N"].coord = np.array((0., 0., 0.))
        chains["B"][1].atoms["O"].coord = np.array((2.9, 0., 0.))
        chains["B"][1].atoms["N"].coord = np.array((0., 5., 0.))
        chains["A"][1].atoms["O"].coord = np.array((2.9, 5., 0.))
        data = {}
        for index, residue in enumerate([*chains["A"], *chains["B"]], 1):
            row = [index, "G", "E", .1, -120., 120., 0, 0., 0, 0., 0, 0., 0, 0.]
            data[(residue.chain_id, (" ", residue.seq_id, " "))] = row
        data[("A", (" ", 1, " "))][8:10] = [2, 0.]
        data[("A", (" ", 2, " "))][6:10] = [2, -1.8, 2, -1.9]
        data[("B", (" ", 2, " "))][6:10] = [-2, -1.9, -2, -1.8]
        with patch.object(stage2, "DSSP", return_value=data):
            stage2.assign_dssp(None, Path("unused.cif"), chains, "unused")
        self.assertEqual(chains["A"][1].dssp_hbonds[0]["role"], "donor")
        self.assertEqual(chains["A"][1].dssp_hbonds[1]["role"], "acceptor")
        sites = stage2.annotate_beta_sheet_sites("TEST", chains, stage2.Stage2Config())
        window = next(site for site in sites if len(site["residue_selectors"]) == 4)
        self.assertEqual(window["motif_type"], "beta_sheet_antiparallel")
        bonds = window["geometry"]["dssp_hydrogen_bonds"]
        self.assertEqual(len(bonds), 2)
        self.assertEqual({bond["energy_kcal_mol"] for bond in bonds}, {-1.8, -1.9})
        self.assertTrue(all(bond["donor"]["atom_name"] == "N" and bond["acceptor"]["atom_name"] == "O" for bond in bonds))
        chains["A"][1].atoms["N"].coord += 50.
        chains["B"][1].atoms["N"].coord += 50.
        self.assertEqual(stage2.annotate_beta_sheet_sites("TEST", chains, stage2.Stage2Config()), [])

    def test_unobserved_polymer_end_is_not_a_genuine_terminus(self):
        source = stage3.SourceIndex.__new__(stage3.SourceIndex)
        key = ("A", 3, "")
        source.chain_polymer_keys = {"A": [key]}
        source.residues = {key: source_record("GLY", 3, {"N": ("N", (0., 0., 0.)), "C": ("C", (1., 0., 0.))})}
        source.structure = SimpleNamespace(entities=[SimpleNamespace(subchains=["A"], full_sequence=["GLY"] * 8)])
        with self.assertRaises(stage3.PreparationError) as error:
            source.biological_neighbors([key])
        self.assertEqual(error.exception.code, "ambiguous_polymer_endpoint")

    def test_disulfide_and_neutral_arg_guards(self):
        source = stage3.SourceIndex.__new__(stage3.SourceIndex)
        source.structure = SimpleNamespace(connections=[])
        source.residues = {
            ("A", 1, ""): source_record("CYS", 1, {"SG": ("S", (0., 0., 0.))}),
            ("A", 2, ""): source_record("CYS", 2, {"SG": ("S", (2.05, 0., 0.))}),
        }
        self.assertTrue(source.has_disulfide(("A", 1, "")))
        source.residues[("A", 3, "")] = source_record("ARG", 3, {})
        with self.assertRaises(stage3.PreparationError) as error:
            stage3._configure_protonation(Chem.RWMol(), {}, source, {("A", 3, "")}, 13., {})
        self.assertEqual(error.exception.code, "unsupported_neutral_arginine")

    def test_environment_discloses_dod_copper_neutral_protein_and_ligands(self):
        source = stage3.SourceIndex.__new__(stage3.SourceIndex)
        source.residues = {
            ("A", 1, ""): source_record("DOD", 1, {"O": ("O", (2.5, 0., 0.)), "D1": ("D", (1.7, 0., 0.))}),
            ("A", 2, ""): source_record("CU", 2, {"CU": ("CU", (3.1, 0., 0.))}),
            ("A", 3, ""): source_record("ALA", 3, {"N": ("N", (3.05, 0., 0.))}),
            ("A", 4, ""): source_record("SO4", 4, {"S": ("S", (4., 0., 0.)), "O": ("O", (4.1, 0., 0.))}),
            ("A", 5, ""): source_record("MSE", 5, {"SE": ("SE", (4.2, 0., 0.))}, gemmi.EntityType.Polymer),
        }
        report = source.nearby_environment(set(), [np.zeros(3)], 3.5)
        self.assertEqual(report["nearby_waters"][0]["minimum_distance_A"], 2.5)
        self.assertEqual(report["nearby_ions"][0]["resname"], "CU")
        self.assertEqual(len(report["nearby_ions"]), 1)
        self.assertEqual(report["nearby_protein_residues"][0]["resname"], "ALA")
        self.assertEqual(report["nearby_ligands"][0]["resname"], "SO4")

    def test_unique_pdb_hydrogen_names_preserve_parent_labels(self):
        molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
        molecule.AddConformer(Chem.Conformer(molecule.GetNumAtoms()))
        heavy = [stage3.FragmentAtom("C", np.zeros(3), f"C{i}", "A", 1, "", "GLY", i, "source_heavy", label=f"A:GLY1:C{i}") for i in (1, 2)]
        records = stage3._prepared_atom_records(molecule, heavy)
        self.assertEqual(len({atom.atom_name for atom in records}), len(records))
        self.assertIn("A:GLY1:C1:H1", {atom.label for atom in records})
        self.assertIn("A:GLY1:C2:H1", {atom.label for atom in records})

    def test_selector_rejects_missing_reports(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = selector.evaluate(Path(temporary), .9, 40.)
        self.assertFalse(result["hard_pass"])
        self.assertIn("validation_passed", result["failed_filters"])

    def test_selector_requires_explicit_validation_matching_maps_and_source_quality(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            reports = {
                "site_metadata": {"branch": "protein", "motif_type": "alpha_helix", "site_id": "TEST",
                                  "quality": {"min_occupancy": 1., "mean_local_bfactor": 12.},
                                  "residue_selectors": [{"resname": "GLY"}]},
                "charge_multiplicity": {"formal_charge": 0, "multiplicity": 1},
                "validation_report": {"validation_pass": True, "atom_count": 2, "formal_charge": 0,
                                      "multiplicity": 1, "heavy_atom_coordinate_mismatches": []},
                "protonation_report": {"pH": 7., "manual_review_required": False},
                "environment_report": {"nearby_waters": [], "nearby_ions": []},
                "cap_report": {"peptide_link_caps": []},
            }
            for name, report in reports.items():
                (folder / f"{name}.json").write_text(json.dumps(report), encoding="utf-8")
            (folder / "gaussian_coordinates.txt").write_text("H 0 0 0 0 ! H1\nH 0 .74 0 0 ! H2\n")
            for filename in ("atom_map.csv", "freeze_mask.csv"):
                with (folder / filename).open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=["fragment_index", "label", "element", "x", "y", "z", "freeze_flag"])
                    writer.writeheader()
                    for index, x in ((1, 0.), (2, .74)):
                        writer.writerow({"fragment_index": index, "label": f"H{index}", "element": "H",
                                         "x": x, "y": 0., "z": 0., "freeze_flag": 0})
            result = selector.evaluate(folder, .9, 40.)
            self.assertTrue(result["hard_pass"], result["failed_filters"])
            self.assertEqual(result["mean_bfactor"], 12.)
            validation_path = folder / "validation_report.json"
            good = reports["validation_report"]
            for bad in ({key: value for key, value in good.items() if key != "validation_pass"},
                        {**good, "heavy_atom_coordinate_mismatches": [{"label": "C", "displacement_A": .1}]},
                        {**good, "atom_count": 3}, {**good, "formal_charge": 1}):
                validation_path.write_text(json.dumps(bad), encoding="utf-8")
                self.assertFalse(selector.evaluate(folder, .9, 40.)["hard_pass"])
            validation_path.write_text(json.dumps(good), encoding="utf-8")
            (folder / "gaussian_coordinates.txt").write_text("H -1 0 0 0 ! H1\nH 0 .74 0 0 ! H2\n")
            self.assertFalse(selector.evaluate(folder, .9, 40.)["hard_pass"])


REFERENCE_CASES = [
    ("alpha_helix", "SITE_4558F38E28E1", "alpha_helix"),
    ("3_10_helix", "SITE_63B5A70D4D4A", "3_10_helix"),
    ("pi_helix", "SITE_295B31E12372", "pi_helix"),
    ("beta_sheet_parallel", "SITE_2550F1BF5D14", "beta_sheet_parallel"),
    ("beta_sheet_antiparallel", "SITE_3351046338C0", None),
    ("beta_turn_I", "SITE_1E27477F21BB", None),
    ("beta_turn_II", "SITE_2DBE5B71C1BA", "beta_turn_II"),
    ("beta_turn_I_prime", "SITE_F9A11D83265F", "beta_turn_I_prime"),
    ("beta_turn_II_prime", "SITE_8B0D3E7F52AC", None),
    ("gamma_turn_classic", "SITE_3F6573CC9241", "gamma_turn_inverse"),
    ("gamma_turn_inverse", "SITE_FFBFE213EF18", "gamma_turn_classic"),
    ("beta_turn_I", "SITE_E7C0D8D87D6B", "beta_turn_I"),
    ("beta_sheet_antiparallel", "SITE_57436A5897F0", "beta_sheet_antiparallel"),
    ("beta_turn_II_prime", "SITE_9E7ED009E3C6", "beta_turn_II_prime"),
    ("pi_helix", "SITE_62D8C06400ED", "pi_helix"),
]


def selector_key(site: dict) -> tuple:
    return tuple(sorted((item["chain_id"], int(item["seq_id"]), item.get("insertion_code", ""))
                        for item in site["residue_selectors"]))


class ReferenceTests(unittest.TestCase):
    root: Path
    dssp: str

    @classmethod
    def setUpClass(cls):
        cls.sources = {}
        cls.annotations = {}
        cls.sites = []
        for motif, site_id, expected in REFERENCE_CASES:
            folder = cls.root / "stage3_protein/fragments/protein" / motif / site_id
            site = json.loads((folder / "site_metadata.json").read_text(encoding="utf-8"))
            pdb_id = site["pdb_id"]
            parent = cls.root / "stage1_gold/protein/mmcif" / f"{pdb_id.lower()}.cif.gz"
            site["source_mmcif"] = str(parent)
            cls.sites.append((folder, site, expected))
            if pdb_id not in cls.sources:
                cls.sources[pdb_id] = stage3.SourceIndex(parent)
                record = {**site, "mmcif_path": str(parent)}
                _, cls.annotations[pdb_id] = stage2.process_structure(record, cls.root / "stage1_gold", stage2.Stage2Config(), cls.dssp)

    def test_fresh_dssp_keeps_valid_examples_and_rejects_reviewed_misassignments(self):
        for folder, site, expected in self.sites:
            with self.subTest(site=site["site_id"]):
                matches = [item for item in self.annotations[site["pdb_id"]] if selector_key(item) == selector_key(site)]
                if expected is None:
                    self.assertFalse(any(item["motif_type"] == site["motif_type"] for item in matches))
                else:
                    self.assertIn(expected, [item["motif_type"] for item in matches])

    def test_preparation_replays_existing_counts_coordinates_freezes_and_charge(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            for folder, original, expected in self.sites:
                with self.subTest(site=original["site_id"]):
                    site = dict(original)
                    if expected is not None:
                        site["motif_type"] = expected
                    prepared = stage3.prepare_site(site, self.sources[site["pdb_id"]], output, stage3.QualityConfig(), stage3.PreparationConfig())
                    reference = json.loads((folder / "validation_report.json").read_text(encoding="utf-8"))
                    self.assertEqual(prepared["atom_count"], reference["atom_count"])
                    self.assertEqual(prepared["formal_charge"], reference["formal_charge"])
                    self.assertEqual(prepared["multiplicity"], reference["multiplicity"])
                    destination = Path(prepared["folder"])
                    with (folder / "atom_map.csv").open(newline="", encoding="utf-8") as handle:
                        before = list(csv.DictReader(handle))
                    with (destination / "atom_map.csv").open(newline="", encoding="utf-8") as handle:
                        after = list(csv.DictReader(handle))
                    self.assertEqual([row["label"] for row in after], [row["label"] for row in before])
                    self.assertEqual([row["element"] for row in after], [row["element"] for row in before])
                    for old, new in zip(before, after):
                        if old["element"] != "H":
                            self.assertTrue(all(abs(float(old[axis]) - float(new[axis])) < 1e-6 for axis in ("x", "y", "z")))
                    selected = selector.evaluate(destination, .9, 40.)
                    self.assertTrue(selected["hard_pass"], selected["failed_filters"])
                    self.assertGreater(selected["mean_bfactor"], 0.)
                    cap_report = json.loads((destination / "cap_report.json").read_text())
                    affected = {label for link in cap_report["peptide_link_caps"] for label in link["affected_heavy_atom_labels"]}
                    self.assertEqual({row["label"] for row in after if row["boundary_affected"] == "True"}, affected)
                    names = [(row["chain_id"], row["seq_id"], row["insertion_code"], row["atom_name"]) for row in after]
                    self.assertEqual(len(set(names)), len(names))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, help="Existing local pipeline output folder; read only")
    parser.add_argument("--dssp", default="mkdssp")
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ProteinTests)
    if args.reference_root:
        ReferenceTests.root, ReferenceTests.dssp = args.reference_root.resolve(), args.dssp
        suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(ReferenceTests))
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
