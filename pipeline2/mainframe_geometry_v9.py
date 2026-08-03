"""Build the model-molecule mainframe using explicit interaction geometry rules."""

import re

import numpy as np
import pandas as pd


SUBSTITUENTS = ["CF3", "Cl", "CN", "F", "H", "NH2", "NMe2", "NO2", "OH"]
CONFORMERS = ["N1", "N2", "P1", "P2"]

atom_sets = {}
coords_dict = {}


def distance(a, b):
    return float(np.linalg.norm(a - b))


def angle_between(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return np.nan
    cosine = np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0)
    return float(np.degrees(np.arccos(cosine)))


def dihedral(p0, p1, p2, p3):
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    norm_b1 = np.linalg.norm(b1)
    if norm_b1 == 0:
        return np.nan

    b1 = b1 / norm_b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if np.linalg.norm(v) == 0 or np.linalg.norm(w) == 0:
        return np.nan

    return float(
        np.degrees(
            np.arctan2(
                np.dot(np.cross(b1, v), w),
                np.dot(v, w),
            )
        )
    )


def absolute_dihedral(p0, p1, p2, p3):
    value = dihedral(p0, p1, p2, p3)
    return abs(value) if not pd.isna(value) else np.nan


def read_xyz(file):
    coords = {}
    with open(file, encoding="utf-8") as handle:
        for index, line in enumerate(handle.readlines()[2:], start=1):
            parts = line.split()
            if len(parts) < 4:
                continue
            atom, x, y, z = parts[:4]
            coords[f"{atom}{index}"] = np.array([float(x), float(y), float(z)])
    return coords


def parse_label(label):
    return re.findall(r"[A-Z][a-z]?\d+", str(label))


def normalize_label_format(label):
    atoms = sorted(
        parse_label(label),
        key=lambda atom: int(re.search(r"\d+", atom).group()),
    )
    return "|".join(atoms)


def strip_serial(label):
    return re.sub(r"^\s*\d+\.\s*", "", str(label)) if pd.notna(label) else label


def normalize_H_atom(atom, number):
    number = int(number)
    if atom == "H" and number == 8:
        number = 28
    elif number > 8:
        number -= 1
    return f"{atom}{number}"


def normalize_H(label):
    return "".join(
        normalize_H_atom(atom, number)
        for atom, number in re.findall(r"([A-Z][a-z]?)(\d+)", str(label))
    )


def clean_orbital_text(text, is_H_substituent=False):
    value = re.sub(r"\s+", " ", str(strip_serial(text))).strip()
    value = re.sub(r"\(\s*(\d+)\s*\)", r"(\1)", value)

    def replace_atom(match):
        atom, number = match.group(1), int(match.group(2))
        return (
            normalize_H_atom(atom, number)
            if is_H_substituent
            else f"{atom}{number}"
        )

    value = re.sub(r"\b([A-Z][a-z]?)\s*(\d+)\b", replace_atom, value)
    return re.sub(r"\s*-\s*", "-", value)


def normalize_interaction_type(value):
    text = str(value).strip().lower()
    text = text.replace("σ", "sig").replace("sigma", "sig")
    text = text.replace("π", "pi")
    text = text.replace("→", "-").replace("->", "-")
    text = text.replace("–", "-").replace("—", "-")
    text = text.replace("∗", "").replace("*", "")
    text = re.sub(r"[\s_]+", "", text)
    return re.sub(r"-+", "-", text)


def element_of(label):
    match = re.match(r"[A-Z][a-z]?", str(label))
    if not match:
        raise ValueError(f"Cannot determine element from atom label: {label}")
    return match.group(0)


def matches_elements(labels, expected):
    return sorted(element_of(label) for label in labels) == sorted(expected)


def atom_of_element(labels, element):
    for label in labels:
        if element_of(label) == element:
            return label
    raise ValueError(f"Expected element {element} in labels: {'|'.join(labels)}")


def two_atoms_of_element(labels, element):
    matches = [label for label in labels if element_of(label) == element]
    if len(matches) != 2:
        raise ValueError(f"Expected two {element} atoms in labels: {'|'.join(labels)}")
    return matches


def atom_distance(label_a, label_b, coords):
    return distance(coords[label_a], coords[label_b])


def nearer_and_other(reference, atom_1, atom_2, coords):
    if atom_distance(reference, atom_1, coords) <= atom_distance(reference, atom_2, coords):
        return atom_1, atom_2
    return atom_2, atom_1


def bonded_carbon_toward_acceptor(donor_atom, acceptor_c, coords):
    """Select a carbon directly bonded to donor N/O and oriented toward acceptor C."""
    donor_element = element_of(donor_atom)
    upper_limit = 1.65 if donor_element == "O" else 1.75
    candidates = [
        atom
        for atom in coords
        if atom != acceptor_c
        and element_of(atom) == "C"
        and 1.10 < atom_distance(donor_atom, atom, coords) < upper_limit
    ]
    if not candidates:
        raise ValueError(f"No carbon bonded directly to donor atom {donor_atom}")
    if len(candidates) == 1:
        return candidates[0]
    return min(candidates, key=lambda atom: atom_distance(atom, acceptor_c, coords))


def geometry_result(distance_value=np.nan, angle_value=np.nan, dihedral_value=np.nan):
    return pd.Series([distance_value, angle_value, dihedral_value])


def geometry_n_pi(donor, acceptor, coords):
    donor_atom = donor[0]
    acceptor_c = atom_of_element(acceptor, "C")
    acceptor_o = atom_of_element(acceptor, "O")
    return geometry_result(
        atom_distance(donor_atom, acceptor_c, coords),
        angle_between(
            coords[donor_atom] - coords[acceptor_c],
            coords[acceptor_o] - coords[acceptor_c],
        ),
    )


def geometry_n_sigma_ch(donor, acceptor, coords):
    donor_atom = donor[0]
    acceptor_c = atom_of_element(acceptor, "C")
    acceptor_h = atom_of_element(acceptor, "H")
    donor_c = bonded_carbon_toward_acceptor(
        donor_atom,
        acceptor_c,
        coords,
    )
    return geometry_result(
        atom_distance(donor_atom, acceptor_h, coords),
        dihedral_value=absolute_dihedral(
            coords[donor_atom],
            coords[donor_c],
            coords[acceptor_c],
            coords[acceptor_h],
        ),
    )


def geometry_n_sigma_cc(donor, acceptor, coords):
    donor_atom = donor[0]
    carbon_1, carbon_2 = two_atoms_of_element(acceptor, "C")
    c_alpha, c_beta = nearer_and_other(donor_atom, carbon_1, carbon_2, coords)
    return geometry_result(
        angle_value=angle_between(
            coords[donor_atom] - coords[c_alpha],
            coords[c_beta] - coords[c_alpha],
        )
    )


def geometry_pi_pi(donor, acceptor, coords):
    donor_c = atom_of_element(donor, "C")
    donor_o = atom_of_element(donor, "O")
    acceptor_c = atom_of_element(acceptor, "C")
    acceptor_o = atom_of_element(acceptor, "O")
    return geometry_result(
        atom_distance(donor_o, acceptor_o, coords),
        dihedral_value=absolute_dihedral(
            coords[donor_o],
            coords[donor_c],
            coords[acceptor_c],
            coords[acceptor_o],
        ),
    )


def geometry_sigma_pi_nc_co(donor, acceptor, coords):
    donor_c = atom_of_element(donor, "C")
    donor_n = atom_of_element(donor, "N")
    acceptor_c = atom_of_element(acceptor, "C")
    acceptor_o = atom_of_element(acceptor, "O")
    return geometry_result(
        atom_distance(donor_c, acceptor_o, coords),
        dihedral_value=dihedral(
            coords[donor_c], coords[donor_n], coords[acceptor_c], coords[acceptor_o]
        ),
    )


def geometry_sigma_pi_ch_co(donor, acceptor, coords):
    donor_c = atom_of_element(donor, "C")
    donor_h = atom_of_element(donor, "H")
    acceptor_c = atom_of_element(acceptor, "C")
    acceptor_o = atom_of_element(acceptor, "O")
    return geometry_result(
        atom_distance(donor_h, acceptor_o, coords),
        dihedral_value=dihedral(
            coords[donor_h], coords[donor_c], coords[acceptor_c], coords[acceptor_o]
        ),
    )


def geometry_sigma_pi_ch_cc(donor, acceptor, coords):
    donor_c = atom_of_element(donor, "C")
    donor_h = atom_of_element(donor, "H")
    carbon_1, carbon_2 = two_atoms_of_element(acceptor, "C")
    contacted_c, other_c = nearer_and_other(donor_h, carbon_1, carbon_2, coords)
    return geometry_result(
        atom_distance(donor_h, contacted_c, coords),
        angle_between(
            coords[donor_c] - coords[contacted_c],
            coords[other_c] - coords[contacted_c],
        ),
    )


def geometry_pi_sigma(donor, acceptor, coords):
    donor_o = atom_of_element(donor, "O")
    acceptor_c = atom_of_element(acceptor, "C")
    acceptor_h = atom_of_element(acceptor, "H")
    return geometry_result(
        angle_value=angle_between(
            coords[acceptor_h] - coords[acceptor_c],
            coords[donor_o] - coords[acceptor_c],
        )
    )


def geometry_sigma_sigma_ch_cc(donor, acceptor, coords):
    donor_c = atom_of_element(donor, "C")
    donor_h = atom_of_element(donor, "H")
    carbon_1, carbon_2 = two_atoms_of_element(acceptor, "C")
    acceptor_near_c, acceptor_far_c = nearer_and_other(
        donor_c, carbon_1, carbon_2, coords
    )
    return geometry_result(
        atom_distance(donor_c, acceptor_far_c, coords),
        dihedral_value=dihedral(
            coords[donor_h],
            coords[donor_c],
            coords[acceptor_near_c],
            coords[acceptor_far_c],
        ),
    )


def geometry_sigma_sigma_ch_ch(donor, acceptor, coords):
    donor_c = atom_of_element(donor, "C")
    donor_h = atom_of_element(donor, "H")
    acceptor_c = atom_of_element(acceptor, "C")
    acceptor_h = atom_of_element(acceptor, "H")
    return geometry_result(
        atom_distance(donor_h, acceptor_h, coords),
        dihedral_value=dihedral(
            coords[donor_h], coords[donor_c], coords[acceptor_c], coords[acceptor_h]
        ),
    )


def geometry_sigma_sigma_ch_cn(donor, acceptor, coords):
    donor_c = atom_of_element(donor, "C")
    donor_h = atom_of_element(donor, "H")
    acceptor_c = atom_of_element(acceptor, "C")
    acceptor_n = atom_of_element(acceptor, "N")
    return geometry_result(
        atom_distance(donor_h, acceptor_c, coords),
        dihedral_value=dihedral(
            coords[donor_h], coords[donor_c], coords[acceptor_n], coords[acceptor_c]
        ),
    )


def compute_geometry(row, coords):
    interaction_type = normalize_interaction_type(row["interaction_type"])
    donor = parse_label(row["donor_labels"])
    acceptor = parse_label(row["acceptor_labels"])
    blank_geometry = geometry_result()

    if any(atom not in coords for atom in donor + acceptor):
        return blank_geometry

    try:
        if interaction_type == "n-pi":
            if (
                len(donor) == 1
                and element_of(donor[0]) in {"N", "O"}
                and matches_elements(acceptor, ("C", "O"))
            ):
                return geometry_n_pi(donor, acceptor, coords)
            return blank_geometry

        if interaction_type in {"n-sig", "n-sigch", "n-sigcc"}:
            valid_donor = len(donor) == 1 and element_of(donor[0]) in {"N", "O"}
            if not valid_donor:
                return blank_geometry

            if interaction_type != "n-sigcc" and matches_elements(
                acceptor, ("C", "H")
            ):
                return geometry_n_sigma_ch(donor, acceptor, coords)

            if interaction_type != "n-sigch" and matches_elements(
                acceptor, ("C", "C")
            ):
                return geometry_n_sigma_cc(donor, acceptor, coords)

            return blank_geometry

        if interaction_type == "pi-pi":
            if matches_elements(donor, ("C", "O")) and matches_elements(
                acceptor, ("C", "O")
            ):
                return geometry_pi_pi(donor, acceptor, coords)
            return blank_geometry

        if interaction_type == "sig-pi":
            if matches_elements(donor, ("C", "N")) and matches_elements(
                acceptor, ("C", "O")
            ):
                return geometry_sigma_pi_nc_co(donor, acceptor, coords)
            if matches_elements(donor, ("C", "H")) and matches_elements(
                acceptor, ("C", "O")
            ):
                return geometry_sigma_pi_ch_co(donor, acceptor, coords)
            if matches_elements(donor, ("C", "H")) and matches_elements(
                acceptor, ("C", "C")
            ):
                return geometry_sigma_pi_ch_cc(donor, acceptor, coords)
            return blank_geometry

        if interaction_type == "pi-sig":
            if matches_elements(donor, ("C", "O")) and matches_elements(
                acceptor, ("C", "H")
            ):
                return geometry_pi_sigma(donor, acceptor, coords)
            return blank_geometry

        if interaction_type == "sig-sig":
            if not matches_elements(donor, ("C", "H")):
                return blank_geometry
            if matches_elements(acceptor, ("C", "C")):
                return geometry_sigma_sigma_ch_cc(donor, acceptor, coords)
            if matches_elements(acceptor, ("C", "H")):
                return geometry_sigma_sigma_ch_ch(donor, acceptor, coords)
            if matches_elements(acceptor, ("C", "N")):
                return geometry_sigma_sigma_ch_cn(donor, acceptor, coords)

        return blank_geometry

    except (KeyError, ValueError):
        # Keep the row blank when a required numbered point is unavailable;
        # do not discard the full model system or invent another descriptor.
        return blank_geometry


def classify_strength(E2):
    if pd.isna(E2):
        return "DE"
    if E2 > 5:
        return "strong"
    if E2 > 1:
        return "moderate"
    return "weak"


def check_possible(row):
    key = (row["substituent"], row["conformer"])
    atoms = parse_label(row.get("donor_labels", "")) + parse_label(
        row.get("acceptor_labels", "")
    )
    return int(all(atom in atom_sets.get(key, set()) for atom in atoms))


def make_interaction_id(row):
    return "__".join(
        str(row[column])
        for column in (
            "donor_orbital",
            "acceptor_orbital",
            "donor_labels",
            "acceptor_labels",
            "interaction_type",
        )
    )


def process(substituent, conformer):
    df = pd.read_csv(f"classified_{substituent}_{conformer}.csv")
    coords = read_xyz(f"{substituent}_{conformer}.xyz")
    is_H = substituent == "H"

    if is_H:
        coords = {normalize_H(label): xyz for label, xyz in coords.items()}

    df["donor_orbital"] = df["donor"].apply(
        lambda value: clean_orbital_text(value, is_H)
    )
    df["acceptor_orbital"] = df["acceptor"].apply(
        lambda value: clean_orbital_text(value, is_H)
    )

    if is_H:
        df["donor_labels"] = df["donor_labels"].apply(normalize_H)
        df["acceptor_labels"] = df["acceptor_labels"].apply(normalize_H)

    df["donor_labels"] = df["donor_labels"].apply(normalize_label_format)
    df["acceptor_labels"] = df["acceptor_labels"].apply(normalize_label_format)
    df["donor"] = df["donor_labels"]
    df["acceptor"] = df["acceptor_labels"]
    df[["distance", "angle", "dihedral"]] = df.apply(
        lambda row: compute_geometry(row, coords), axis=1
    )
    df["substituent"], df["conformer"] = substituent, conformer
    return df


def fill_geometry(row):
    if row["possibility_flag"] and pd.isna(row.get("distance")):
        coords = coords_dict.get((row["substituent"], row["conformer"]))
        if coords:
            return compute_geometry(row, coords)
    return pd.Series([row.get("distance"), row.get("angle"), row.get("dihedral")])


def main():
    all_data = []

    for substituent in SUBSTITUENTS:
        for conformer in CONFORMERS:
            try:
                print(f"Processing: {substituent}_{conformer}")
                coords = read_xyz(f"{substituent}_{conformer}.xyz")
                if substituent == "H":
                    coords = {normalize_H(label): xyz for label, xyz in coords.items()}

                atom_sets[(substituent, conformer)] = set(coords)
                coords_dict[(substituent, conformer)] = coords
                all_data.append(process(substituent, conformer))
            except Exception as error:
                print(f"Skipped {substituent}_{conformer}: {error}")

    if not all_data:
        raise ValueError("No model-molecule systems were processed.")

    raw_df = pd.concat(all_data, ignore_index=True)
    raw_df["E2_kcal/mol"] = pd.to_numeric(raw_df["E2_kcal/mol"], errors="coerce")

    interaction_cols = [
        "donor_orbital",
        "acceptor_orbital",
        "donor_labels",
        "acceptor_labels",
        "interaction_type",
    ]
    sample_cols = ["substituent", "conformer"]
    merge_cols = sample_cols + interaction_cols
    raw_df["interaction_id"] = raw_df.apply(make_interaction_id, axis=1)

    exact_dup_cols = merge_cols + ["E2_kcal/mol"]
    exact_duplicates = raw_df[
        raw_df.duplicated(subset=exact_dup_cols, keep=False)
    ].copy()
    if not exact_duplicates.empty:
        exact_duplicates.to_csv("exact_duplicate_entries_removed.csv", index=False)

    raw_df = raw_df.drop_duplicates(subset=exact_dup_cols, keep="first").copy()
    conflicts = raw_df[raw_df.duplicated(subset=merge_cols, keep=False)].copy()
    if not conflicts.empty:
        conflicts.sort_values(merge_cols + ["E2_kcal/mol"]).to_csv(
            "unresolved_duplicate_conflicts.csv", index=False
        )
        raise ValueError(
            "Unresolved duplicate conflicts found. "
            "Check unresolved_duplicate_conflicts.csv."
        )

    print("Duplicate entries after cleanup:", raw_df.duplicated(subset=merge_cols).sum())

    all_interactions = raw_df[
        interaction_cols + ["interaction_id", "bond_distance"]
    ].drop_duplicates(subset=interaction_cols)
    full_grid = raw_df[sample_cols].drop_duplicates().merge(all_interactions, how="cross")
    merged = full_grid.merge(raw_df, on=merge_cols, how="left", suffixes=("_ref", ""))

    merged["interaction_id"] = merged["interaction_id_ref"].combine_first(
        merged.get("interaction_id")
    )
    if "bond_distance_ref" in merged:
        merged["bond_distance"] = merged["bond_distance"].combine_first(
            merged["bond_distance_ref"]
        )

    merged["exists_flag"] = merged["E2_kcal/mol"].notna().astype(int)
    merged["possibility_flag"] = merged.apply(check_possible, axis=1)
    merged[["distance", "angle", "dihedral"]] = merged.apply(fill_geometry, axis=1)
    merged["interaction_strength"] = merged.apply(
        lambda row: (
            "NA"
            if row["possibility_flag"] == 0
            else classify_strength(row["E2_kcal/mol"])
        ),
        axis=1,
    )

    final_dup_count = merged.duplicated(subset=merge_cols).sum()
    print("Duplicate entries in final mainframe:", final_dup_count)
    if final_dup_count:
        merged[merged.duplicated(subset=merge_cols, keep=False)].to_csv(
            "final_duplicate_check_failed.csv", index=False
        )
        raise ValueError(
            "Final mainframe still has duplicate rows. "
            "Check final_duplicate_check_failed.csv."
        )

    print("\nRows per substituent/conformer:")
    print(merged.groupby(sample_cols).size().unstack())
    print("\nUnique interactions in universe:", merged["interaction_id"].nunique())

    final_columns = [
        "interaction_id",
        "E2_kcal/mol",
        "donor_orbital",
        "acceptor_orbital",
        "donor_labels",
        "acceptor_labels",
        "bond_distance",
        "interaction_type",
        "distance",
        "angle",
        "dihedral",
        "substituent",
        "conformer",
        "exists_flag",
        "possibility_flag",
        "interaction_strength",
    ]
    merged[final_columns].to_csv("mainframe.csv", index=False)
    print("DONE: mainframe.csv created")


if __name__ == "__main__":
    main()
