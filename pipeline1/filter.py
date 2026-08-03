import argparse
import csv
import math
import os
import re
from collections import deque

COVALENT_RADII = {
    "H": 0.37, "C": 0.77, "N": 0.75, "O": 0.73,
    "F": 0.71, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Br": 1.20, "I": 1.39, "Si": 1.11, "B": 0.85,
    "Na": 1.66, "K": 2.03, "Ca": 1.74, "Fe": 1.26,
    "Zn": 1.22, "Cu": 1.32, "Mg": 1.41, "Al": 1.18,
}


def parse_xyz(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    natoms = int(lines[0].strip())
    atom_lines = lines[2:2 + natoms]  # line 2 is the XYZ comment, even when blank
    if len(atom_lines) != natoms:
        raise ValueError(f"XYZ atom count mismatch in {path}")
    return [
        (parts[0], *map(float, parts[1:4]))
        for parts in (line.split() for line in atom_lines)
    ]


def parse_pdb(path):
    coords = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip() or re.sub(r"[^A-Za-z]", "", line[12:16])[0]
                coords.append((
                    element.capitalize(),
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ))
    return coords


def parse_structure(path):
    extension = os.path.splitext(path)[1].lower()
    if extension == ".xyz":
        return parse_xyz(path)
    if extension == ".pdb":
        return parse_pdb(path)
    raise ValueError(f"Unsupported file format: {extension}")


def label_atoms(coords):
    return {i: f"{element}{i}" for i, (element, _, _, _) in enumerate(coords, 1)}


def build_bond_graph(coords, add_tol=0.45, verbose=False, labels=None):
    graph = {i: set() for i in range(1, len(coords) + 1)}
    for i, (element_i, xi, yi, zi) in enumerate(coords):
        if element_i not in COVALENT_RADII:
            raise ValueError(f"Unsupported element: {element_i}")
        for j, (element_j, xj, yj, zj) in enumerate(coords[i + 1:], i + 1):
            if element_j not in COVALENT_RADII:
                raise ValueError(f"Unsupported element: {element_j}")
            distance = math.dist((xi, yi, zi), (xj, yj, zj))
            if distance <= COVALENT_RADII[element_i] + COVALENT_RADII[element_j] + add_tol:
                graph[i + 1].add(j + 1)
                graph[j + 1].add(i + 1)
                if verbose and labels:
                    print(f"Bond: {labels[i + 1]}-{labels[j + 1]} ({distance:.2f} Å)")
    return graph


def parse_indices(label):
    numbers = re.findall(r"\d+", label or "")
    if not numbers:
        return []
    return [int(value) for value in numbers[-2:]] if "-" in (label or "") else [int(numbers[-1])]


def orbital(label):
    match = re.search(r"(LP|BD\*?|RY)\s*\(\s*(\d+)\s*\)", label.upper())
    return (match.group(1), int(match.group(2))) if match else (None, None)


def supported_interaction(donor, acceptor):
    donor_type, donor_number = orbital(donor)
    acceptor_type, acceptor_number = orbital(acceptor)
    donor_ok = donor_type == "LP" or (donor_type == "BD" and donor_number in (1, 2))
    return donor_ok and acceptor_type == "BD*" and acceptor_number in (1, 2)


def bond_distance(graph, atoms_a, atoms_b):
    if not atoms_a or not atoms_b:
        return None
    targets = set(atoms_b)
    queue = deque((atom, 0) for atom in atoms_a if atom in graph)
    visited = {atom for atom in atoms_a if atom in graph}
    while queue:
        atom, distance = queue.popleft()
        if atom in targets:
            return distance
        for neighbour in graph[atom]:
            if neighbour not in visited:
                visited.add(neighbour)
                queue.append((neighbour, distance + 1))
    return None


def main():
    parser = argparse.ArgumentParser(description="Filter NBO interactions")
    parser.add_argument("files", nargs="+")
    parser.add_argument("--min-e2", type=float, default=0.1)
    parser.add_argument("--min-bond-distance", type=int, default=3)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    csv_files = [path for path in args.files if path.lower().endswith(".csv")]
    structure_files = {
        os.path.splitext(os.path.basename(path))[0]: path
        for path in args.files
        if path.lower().endswith((".xyz", ".pdb"))
    }
    headers = [
        "donor", "acceptor", "E2_kcal/mol",
        "donor_labels", "acceptor_labels", "bond_distance",
    ]
    for interactions_csv in csv_files:
        base = os.path.splitext(os.path.basename(interactions_csv))[0]
        structure_file = structure_files.get(base)
        if structure_file is None:
            raise ValueError(f"No matching structure found for {interactions_csv}")
        coords = parse_structure(structure_file)
        labels = label_atoms(coords)
        graph = build_bond_graph(coords, verbose=args.verbose, labels=labels)
        with open(interactions_csv, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        out_rows = []
        for row in rows:
            if float(row["E2_kcal/mol"]) < args.min_e2:
                continue
            if not supported_interaction(row["donor"], row["acceptor"]):
                continue
            donor = parse_indices(row["donor"])
            acceptor = parse_indices(row["acceptor"])
            distance = bond_distance(graph, donor, acceptor)
            #Disconnected fragments have no covalent path. The original pipeline correctly kept such interactions as long-range.
            if distance is not None and distance < args.min_bond_distance:
                continue
            out_rows.append({
                "donor": row["donor"],
                "acceptor": row["acceptor"],
                "E2_kcal/mol": row["E2_kcal/mol"],
                "donor_labels": "|".join(labels[i] for i in donor),
                "acceptor_labels": "|".join(labels[i] for i in acceptor),
                "bond_distance": "" if distance is None else distance,
            })
        out = f"filtered_{base}.csv"
        with open(out, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader()
            writer.writerows(out_rows)
        print(f"Wrote {len(out_rows)} rows to {out}")


if __name__ == "__main__":
    main()
